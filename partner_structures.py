"""Experimental human protein–protein interfaces using assembly coordinates and SIFTS."""

from __future__ import annotations

import gzip
import io
import json
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
from biotite.structure.io import pdbx
from Bio.Data.PDBData import protein_letters_3to1
from scipy.spatial import cKDTree

from binding_sites import AMINO_ACIDS


PDBE = "https://www.ebi.ac.uk/pdbe/"


def mapped_structures(cache, accession, positions=()):
    url = PDBE + "api/mappings/best_structures/" + accession

    def parse(data):
        obj = json.loads(data)
        if not isinstance(obj, dict) or not isinstance(obj.get(accession), list):
            raise ValueError("Unexpected PDBe best_structures response")
        records = obj[accession]
        for row in records:
            if not isinstance(row, dict) or not isinstance(row.get("pdb_id"), str):
                raise ValueError("Missing PDB ID in PDBe response")
        return records

    records = cache.read(
        f"discovery/{accession}.json", url, parse, missing=json.dumps({accession: []}).encode()
    )

    def overlap(row):
        lo, hi = row.get("unp_start") or 0, row.get("unp_end") or 0
        return sum(lo <= pos + 1 <= hi for pos in positions)

    return list(dict.fromkeys(row["pdb_id"].lower() for row in sorted(records, key=overlap, reverse=True)))


def load_entry(cache, pdb_id):
    def parse(data):
        with gzip.open(io.BytesIO(data), "rt", encoding="utf-8") as handle:
            block = pdbx.CIFFile.read(handle).block
        if block["entry"]["id"].as_item().lower() != pdb_id or "entity_poly" not in block:
            raise ValueError(f"Unexpected PDBe mmCIF for {pdb_id}")
        return block

    return cache.read(
        "pdbe/" + pdb_id + ".cif.gz", PDBE + f"entry-files/download/{pdb_id}_updated.cif.gz", parse
    )


def columns(block, category, *keys):
    if category not in block:
        return ()
    table = block[category]
    return zip(*(table[key].as_array(str) for key in keys))


@dataclass
class Chain:
    accession: str
    sequence: str
    name: str
    polymer_form: str
    intervals: list = field(default_factory=list)
    positions: dict = field(default_factory=dict)
    modified_residues: dict = field(default_factory=dict)


def human_chains(block, proteins, parent):
    poly = {
        e: ("".join(seq.split()).upper(), "".join(form.split()))
        for e, t, seq, form in columns(
            block,
            "entity_poly",
            "entity_id",
            "type",
            "pdbx_seq_one_letter_code_can",
            "pdbx_seq_one_letter_code",
        )
        if t == "polypeptide(L)"
    }
    names = dict(columns(block, "entity", "id", "pdbx_description"))
    modified = defaultdict(dict)
    for entity, num, monomer in columns(block, "entity_poly_seq", "entity_id", "num", "mon_id"):
        if num.isdigit() and monomer not in protein_letters_3to1:
            modified[entity][int(num)] = monomer
    taxa = defaultdict(set)
    for cat, key in [
        ("entity_src_gen", "pdbx_gene_src_ncbi_taxonomy_id"),
        ("entity_src_nat", "pdbx_ncbi_taxonomy_id"),
        ("pdbx_entity_src_syn", "ncbi_taxonomy_id"),
    ]:
        for entity, taxon in columns(block, cat, "entity_id", key):
            taxa[entity].add(taxon)
    segments = defaultdict(list)
    for ent, asym, acc, start, end, best in columns(
        block,
        "pdbx_sifts_unp_segments",
        "entity_id",
        "asym_id",
        "unp_acc",
        "unp_start",
        "unp_end",
        "best_mapping",
    ):
        if ent in poly and taxa[ent] == {"9606"} and start.isdigit() and end.isdigit():
            segments[asym].append((acc, ent, int(start), int(end), best == "y"))
    chains = {}
    for asym, records in segments.items():
        # A fusion may map to two different proteins; do not attribute it to either.
        if len({r[0].split("-")[0] for r in records}) != 1:
            continue
        choices = [r for r in records if r[0] in proteins]
        if not choices:
            continue
        acc, ent, *_ = min(choices, key=lambda r: (r[0] != parent.accession, not r[4], "-" in r[0], r[0]))
        seq, form = poly[ent]
        chains[asym] = Chain(
            acc, seq, names.get(ent, ""), form, [(r[2], r[3]) for r in choices if r[0] == acc]
        )
        chains[asym].modified_residues = modified[ent]
    for asym, label, acc, pos, residue in columns(
        block, "pdbx_sifts_xref_db", "asym_id", "seq_id", "unp_acc", "unp_num", "unp_res"
    ):
        chain = chains.get(asym)
        if not chain or acc != chain.accession or not label.isdigit() or not pos.isdigit():
            continue
        pos, label = int(pos), int(label)
        reference = proteins[acc].sequence
        if 1 <= pos <= len(reference) and 1 <= label <= len(chain.sequence) and reference[pos - 1] == residue:
            previous = chain.positions.get(label, pos)
            chain.positions[label] = pos if pos == previous else 0
    return chains


def site_result(contacts, chains, query, parent, match, max_changes, min_identity):
    """Count union of observed contact sites across all contacting copies and target segments."""
    sites = defaultdict(set)
    missing = set()
    unknown = set()
    modifications = set()
    for asym, labels in contacts.items():
        chain = chains[asym]
        for label in labels:
            if label in chain.modified_residues:
                modifications.add(f"{asym}:{label}:{chain.modified_residues[label]}")
            pos = chain.positions.get(label)
            # An exact deposited query also supplies its own unambiguous label coordinates.
            if chain.sequence == query:
                qpos = label - 1
            elif chain.accession == parent.accession and pos:
                qpos = match.mapping.get(pos - 1)
            else:
                qpos = None
            if qpos is None or not 0 <= qpos < len(query):
                missing.add((chain.accession, pos) if pos else (asym, label))
                continue
            aa = chain.sequence[label - 1]
            sites[qpos].add(aa)
            if aa not in AMINO_ACIDS or query[qpos] not in AMINO_ACIDS:
                unknown.add(qpos)
    changes = sum(any(aa != query[q] for aa in vals) for q, vals in sites.items())
    total = len(sites) + len(missing)
    metrics = dict(
        site_size=total,
        site_covered=len(sites),
        site_coverage=len(sites) / total,
        site_substitutions=changes,
        site_identity=(len(sites) - changes) / len(sites) if sites else "",
        query_contact_positions="|".join(str(q + 1) for q in sorted(sites)),
        target_modifications="|".join(sorted(modifications)),
    )
    seqs = {chains[c].sequence for c in contacts}
    if missing:
        reason = "additional_target_chain_outside_query" if len(seqs) > 1 else "site_not_fully_mapped"
    elif len(seqs) > 1:
        reason = "additional_target_chain_required"
    elif unknown:
        reason = "unknown_contact_residue"
    elif modifications:
        reason = "modified_target_contact"
    elif changes > max_changes or metrics["site_identity"] < min_identity:
        reason = "site_difference_limit"
    else:
        tier = "exact_structure" if seqs == {query} else ("site_similar" if changes else "site_conserved")
        return dict(**metrics, status=tier, reason="")
    return dict(**metrics, status="structural_context", reason=reason)


def chain_instances(atoms):
    return "|".join(
        sorted({f"{c}[{a}]#{s}" for c, a, s in zip(atoms.chain_id, atoms.auth_asym_id, atoms.sym_id)})
    )


def partner_components(atoms, partner_ids, indices, cutoff):
    """Keep contacting subunits of the same precursor together, including a non-contacting peptide chain."""
    parts = [
        idx[atoms.sym_id[idx] == sym]
        for chain in partner_ids
        for idx in [indices[chain]]
        for sym in np.unique(atoms.sym_id[idx])
    ]
    if len(parts) <= 1:
        return parts
    trees = [cKDTree(atoms.coord[idx]) for idx in parts]
    neighbours = defaultdict(set)
    for i, tree in enumerate(trees):
        for j in range(i):
            if tree.count_neighbors(trees[j], cutoff):
                neighbours[i].add(j)
                neighbours[j].add(i)
    remaining = set(range(len(parts)))
    groups = []
    while remaining:
        group, frontier = set(), {min(remaining)}
        while frontier:
            group.update(frontier)
            frontier = set().union(*(neighbours[i] for i in frontier)) - group
        remaining -= group
        groups.append(np.concatenate([parts[i] for i in sorted(group)]))
    return groups


def assembly_rows(atoms, chains, targets, query, parent, match, common, cutoff, max_changes, min_identity):
    atoms = atoms[(atoms.element != "H") & (atoms.element != "D") & (atoms.res_id > 0)]
    indices = {chain: np.flatnonzero(atoms.chain_id == chain) for chain in chains}
    target_ids = [chain for chain in targets if indices[chain].size]
    partners = defaultdict(list)
    for chain, item in chains.items():
        if chain not in targets and indices[chain].size:
            partners[item.accession].append(chain)
    if not target_ids:
        return []
    rows = []
    target_atoms = atoms[np.concatenate([indices[c] for c in target_ids])]
    components = [
        (acc, idx)
        for acc, ids in sorted(partners.items())
        for idx in partner_components(atoms, ids, indices, cutoff)
    ]
    assembly_accessions = "|".join(sorted({chains[c].accession for c in indices if indices[c].size}))
    for accession, pidx in components:
        partner_atoms = atoms[pidx]
        ptree = cKDTree(partner_atoms.coord)
        close = ptree.query(target_atoms.coord, k=1)[0] <= cutoff
        if not close.any():
            continue
        hit = target_atoms[close]
        target_tree = cKDTree(hit.coord)
        p_hit = partner_atoms[target_tree.query(partner_atoms.coord, k=1)[0] <= cutoff]
        contacts = {c: set(hit.res_id[hit.chain_id == c].tolist()) for c in sorted(set(hit.chain_id))}
        binding_ids = sorted(set(partner_atoms.chain_id))
        seqs = sorted({chains[c].sequence for c in binding_ids})
        context = "|".join(
            f"{c}:{chains[c].accession}:" + ",".join(f"{start}-{end}" for start, end in chains[c].intervals)
            for c in contacts
        )
        modifications = "|".join(
            sorted(
                {chains[c].polymer_form for c in binding_ids if chains[c].polymer_form != chains[c].sequence}
            )
        )
        binding_type = "peptide" if all(len(s) <= 50 for s in seqs) else "protein"
        stoichiometry = [
            sum(
                len(set(partner_atoms.sym_id[partner_atoms.chain_id == c]))
                for c in binding_ids
                if chains[c].sequence == seq
            )
            for seq in seqs
        ]
        row = dict(
            **common,
            partner_id=accession,
            partner_name=" | ".join(dict.fromkeys(chains[c].name for c in binding_ids)),
            partner_type=binding_type,
            partner_sequences="|".join(seqs),
            sequence_basis="deposited_polypeptides",
            target_chains=chain_instances(hit),
            partner_chains=chain_instances(partner_atoms),
            complex_context=context,
            assembly_components=assembly_accessions,
            biological_role="observed_protein_partner",
            modifications=modifications,
            partner_stoichiometry="|".join(map(str, stoichiometry)),
            partner_chain_sequences=json.dumps({c: chains[c].sequence for c in binding_ids}, sort_keys=True),
        )
        row.update(site_result(contacts, chains, query, parent, match, max_changes, min_identity))
        if any(not chains[c].positions.get(label) for c, label in set(zip(p_hit.chain_id, p_hit.res_id))):
            row.update(status="structural_context", reason="partner_contact_unmapped")
        for other, other_idx in components:
            if other == accession:
                continue
            other_tree = cKDTree(atoms.coord[other_idx])
            if (
                row["status"] != "structural_context"
                and other_tree.count_neighbors(ptree, cutoff)
                and other_tree.count_neighbors(target_tree, cutoff)
            ):
                row.update(status="structural_context", reason="multicomponent_interface")
                break
        if any(not seq or not set(seq) <= AMINO_ACIDS for seq in seqs):
            row.update(status="structural_context", reason="unresolved_partner_sequence")
        rows.append(row)
    return rows


def structure_rows(block, pdb_id, proteins, parent, query, match, cutoff, max_changes, min_identity):
    methods = {method for (method,) in columns(block, "exptl", "method")}
    experimental = {
        "X-RAY DIFFRACTION",
        "NEUTRON DIFFRACTION",
        "ELECTRON MICROSCOPY",
        "ELECTRON CRYSTALLOGRAPHY",
        "SOLUTION NMR",
        "SOLID-STATE NMR",
        "FIBER DIFFRACTION",
    }
    base = dict(
        source="PDBe",
        record_id=pdb_id,
        citation=f"https://www.rcsb.org/structure/{pdb_id.upper()}",
        method="|".join(sorted(methods)),
    )
    if not {method.upper() for method in methods} & experimental:
        return [dict(**base, status="structure_review", reason="experimental_atomic_coordinates_unavailable")]
    if "citation" in block and "pdbx_database_id_PubMed" in block["citation"]:
        for key, pmid in columns(block, "citation", "id", "pdbx_database_id_PubMed"):
            if key.lower() == "primary" and pmid.isdigit():
                base["reference"] = "pubmed:" + pmid
    chains = human_chains(block, proteins, parent)
    targets = [c for c, item in chains.items() if item.accession.split("-")[0] == parent.canonical]
    if not targets:
        return [dict(**base, status="structure_review", reason="no_human_target_mapping")]
    if all(chains[c].accession == parent.accession for c in targets):
        mapped = set(match.mapping)
        if not any(
            any(p - 1 in mapped for p in range(lo, hi + 1)) for c in targets for lo, hi in chains[c].intervals
        ):
            spans = "|".join(f"{c}:{lo}-{hi}" for c in targets for lo, hi in chains[c].intervals)
            return [
                dict(
                    **base,
                    status="outside_query",
                    reason="target_domain_outside_query",
                    complex_context=spans,
                )
            ]
    assemblies = pdbx.list_assemblies(block)
    if not assemblies:
        return [dict(**base, status="structure_review", reason="no_biological_assembly")]
    rows = []
    for assembly in assemblies:
        common = dict(base, assembly_id=str(assembly))
        atoms = pdbx.get_assembly(
            block,
            assembly_id=assembly,
            model=1,
            altloc="occupancy",
            use_author_fields=False,
            extra_fields=["auth_asym_id"],
        )
        found = assembly_rows(
            atoms, chains, targets, query, parent, match, common, cutoff, max_changes, min_identity
        )
        rows.extend(
            found or [dict(**common, status="structure_review", reason="no_other_human_polypeptide_contact")]
        )
    return rows
