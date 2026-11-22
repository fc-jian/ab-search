"""Public reference data and published human protein partner records."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

from binding_sites import AMINO_ACIDS, Match, sequence_match


PER_PAGE = 200
UNIPROT = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
    {
        "format": "fasta",
        "compressed": "true",
        "includeIsoform": "true",
        "query": "(proteome:UP000005640) AND (reviewed:true)",
    }
)
CPDB = "https://raw.githubusercontent.com/ventolab/cellphonedb-data/v5.0.0/data/"
SYNBIP = "https://synbip.idrblab.net/sites/default/files/synbip_download/"
INTACT = "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/"
GTOP = "https://www.guidetopharmacology.org/"
AFFINITY_FIELDS = {"kd": "kd_nM", "ki": "ki_nM", "ic50": "ic50_nM"}
CONCENTRATION_NM = {"pm": 0.001, "nm": 1, "um": 1000, "µm": 1000, "μm": 1000, "mm": 1e6, "m": 1e9}


def measurement(parameter, value, unit, relation):
    field = AFFINITY_FIELDS.get(parameter.casefold())
    scale = CONCENTRATION_NM.get(unit.lower())
    if not field or scale is None:
        return {}
    try:
        nm = float(value) * scale
    except (ValueError, TypeError):
        return {}
    if not math.isfinite(nm) or nm <= 0:
        return {}
    return {field: nm, "measurement_relation": relation or "="}


class Cache:
    def __init__(self, root, offline=False):
        self.root = Path(root)
        self.offline = offline
        self.used = {}

    def read(self, name, url, parse, *, headers=None, missing=None):
        path = self.root / name
        if path.is_file():
            data = path.read_bytes()
            try:
                result = parse(data)
            except (ValueError, OSError, EOFError, KeyError, csv.Error):
                if self.offline:
                    raise
            else:
                self.record(name, url, data)
                return result
        if self.offline:
            raise OSError(f"offline cache miss: {name}")
        request = urllib.request.Request(url, headers={"User-Agent": "ab-search/1.0", **(headers or {})})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request, timeout=35) as response:
                    data = response.read()
                break
            except urllib.error.HTTPError as exc:
                if exc.code == 404 and missing is not None:
                    data = missing
                    break
                if exc.code not in (404, 429, 500, 502, 503, 504) or attempt == 2:
                    raise
            except (urllib.error.URLError, TimeoutError):
                if attempt == 2:
                    raise
            time.sleep(1 + attempt)
        result = parse(data)
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=path.parent) as tmp:
            staging = Path(tmp) / "download"
            staging.write_bytes(data)
            staging.replace(path)
        self.record(name, url, data)
        return result

    def record(self, name, url, data):
        self.used[name] = {"url": url, "sha256": hashlib.sha256(data).hexdigest()}


def csv_records(data, required):
    lines = io.StringIO(data.decode("utf-8-sig"), newline="")
    first = lines.readline()
    if not first.lstrip('"').startswith("#"):
        lines.seek(0)
    reader = csv.DictReader(lines)
    if not set(required) <= set(reader.fieldnames or ()):
        raise ValueError(f"Unexpected source columns, expected {', '.join(required)}")
    return list(reader)


@dataclass(frozen=True)
class Protein:
    accession: str
    gene: str
    name: str
    sequence: str

    @property
    def canonical(self):
        return self.accession.split("-")[0]


def read_proteome(data):
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    proteins = {}
    for head, seq in SimpleFastaParser(io.StringIO(data.decode("utf-8"))):
        if not head.startswith("sp|") or not re.search(r"\bOX=9606\b", head):
            continue
        acc = head.split("|", 2)[1]
        gene = re.search(r"\bGN=(\S+)", head)
        name = head.partition(" ")[2].partition(" OS=")[0]
        proteins[acc] = Protein(acc, gene[1] if gene else acc.split("-")[0], name, seq.upper())
    if not proteins:
        raise ValueError("Expected a reviewed human UniProt FASTA with sp| headers and OX=9606")
    return proteins


def resolve_query(query, proteins):
    """Return accession, parent-to-query placement and all equivalent accessions."""
    if not set(query) <= AMINO_ACIDS:
        return None, None, [], "unresolved_query_residue"
    exact = [p for p in proteins.values() if query in p.sequence]
    if exact:
        options = sorted(exact, key=lambda p: (p.sequence != query, "-" in p.accession, p.accession))
    else:
        # Require two independent exact anchors and near end-to-end identity for variant assignment.
        if len(query) < 24:
            return None, None, [], "no_reviewed_human_sequence"
        seeds = [(i, query[i : i + 8]) for i in range(0, len(query) - 7, max(8, len(query) // 8))]
        words = {s for _, s in seeds}
        scored = []
        for p in proteins.values():
            if sum(seed in p.sequence for seed in words) < 2:
                continue
            # Include terminal substitutions when the full query fits a unique reference window.
            offsets = set()
            for i, seed in seeds:
                hit = p.sequence.find(seed)
                while hit >= 0:
                    if 0 <= hit - i <= len(p.sequence) - len(query):
                        offsets.add(hit - i)
                    hit = p.sequence.find(seed, hit + 1)
            windows = [
                (sum(a == b for a, b in zip(query, p.sequence[j : j + len(query)])), j) for j in offsets
            ]
            score = max((count for count, _ in windows), default=0) / len(query)
            if score >= 0.95:
                starts = [j for count, j in windows if count / len(query) == score]
                m = Match("variant", {starts[0] + i: i for i in range(len(query))}, score, len(starts) != 1)
            else:
                m = sequence_match(query, p.sequence, allow_local=True, min_identity=0.95)
            if m and len(m.mapping) == len(query):
                scored.append((m.identity, p, m))
        if not scored:
            return None, None, [], "no_reviewed_human_sequence"
        best = max(identity for identity, _, _ in scored)
        options = sorted(
            (p for identity, p, _ in scored if identity == best),
            key=lambda p: ("-" in p.accession, p.accession),
        )
    ids = [p.accession for p in options]
    if len({p.canonical for p in options}) != 1:
        return None, None, ids, "ambiguous_human_parent"
    parent = options[0]
    match = (
        sequence_match(query, parent.sequence, allow_local=False, min_identity=0.95)
        if exact
        else next(m for _, p, m in scored if p == parent)
    )
    if not match or match.ambiguous:
        return None, None, ids, "ambiguous_query_position"
    return parent, match, ids, ""


def reference_partner(accessions, proteins):
    ids = sorted(set(accessions))
    if not ids or any(acc not in proteins for acc in ids):
        return None
    items = [proteins[acc] for acc in ids]
    return dict(
        partner_id="|".join(ids),
        partner_name="|".join(p.name for p in items),
        partner_type="protein_complex" if len(ids) > 1 else "protein",
        partner_sequences="|".join(p.sequence for p in items),
        sequence_basis="UniProt_reference",
    )


def cellphonedb(cache, proteins):
    complexes = cache.read(
        "cellphonedb/v5.0.0/complex_input.csv",
        CPDB + "complex_input.csv",
        lambda d: csv_records(d, ["complex_name", "uniprot_1", "uniprot_2"]),
    )
    subunits = {
        c["complex_name"]: [v for k, v in c.items() if k.startswith("uniprot_") and v]
        for c in complexes
        if c.get("other_desc", "").casefold() != "biosynthesis_enzyme"
    }
    records = cache.read(
        "cellphonedb/v5.0.0/interaction_input.csv",
        CPDB + "interaction_input.csv",
        lambda d: csv_records(d, ["partner_a", "partner_b", "is_ppi", "source", "directionality"]),
    )
    out = defaultdict(list)
    for line, r in enumerate(records, 2):
        if r["is_ppi"].lower() != "true":
            continue
        a, b = r["partner_a"], r["partner_b"]
        for left, right, side in [(a, b, "a"), (b, a, "b")] if a != b else [(a, b, "a")]:
            partner = reference_partner(subunits.get(right, [right]), proteins)
            if not partner:
                continue
            if right in subunits:
                partner.update(partner_name=right, partner_type="protein_complex")
            role = r["directionality"].split("-")
            query_role = role[0 if side == "a" else 1] if len(role) == 2 else ""
            for acc in subunits.get(left, [left]):
                out[acc].append(
                    dict(
                        **partner,
                        source="CellPhoneDB",
                        record_id=f"v5.0.0:{line}",
                        status="parent_only",
                        reason="query_complex_subunit" if left in subunits else "construct_unspecified",
                        query_role=query_role,
                        biological_role=r.get("classification", ""),
                        complex_context="|".join(x for x in (left, right) if x in subunits),
                        method=r["annotation_strategy"],
                        reference=r["source"],
                        citation="https://github.com/ventolab/cellphonedb-data/blob/v5.0.0/data/interaction_input.csv",
                    )
                )
    return out


def synbip(cache, proteins):
    catalog_name = "General_Information_of_SBP.csv"
    catalog = cache.read(
        "synbip/" + catalog_name,
        SYNBIP + catalog_name,
        lambda d: csv_records(d, ["SBPID", "Protein Scaffold Name"]),
    )
    scaffold_names = {r["SBPID"]: r["Protein Scaffold Name"] for r in catalog}
    supported = {
        "monobody",
        "adnectin",
        "affibody",
        "darpin",
        "anticalin",
        "affimer",
        "affilin",
        "affitin",
        "centyrin",
        "fynomer",
        "repebody",
        "designed ankyrin repeat protein",
        "peptide aptamer",
        "macrocyclic peptide",
        "abd-derived affinity protein",
    }
    mapping_name = "Target_to_SBP_Mapping_with_Mechanism_of_Action.csv"
    records = cache.read(
        "synbip/" + mapping_name,
        SYNBIP + mapping_name,
        lambda d: csv_records(d, ["SBP ID", "SBP Name", "Affinity", "Target Uniprot ID"]),
    )
    filename = "Sequence_Information_of_SBPs.txt"

    def parse_seqs(data):
        seqs = {
            head.split("|")[0]: "".join(seq.split()).upper()
            for head, seq in SimpleFastaParser(io.StringIO(data.decode("utf-8-sig")))
        }
        if not seqs:
            raise ValueError("Empty SYNBIP sequence file")
        return seqs

    seqs = cache.read("synbip/" + filename, SYNBIP + filename, parse_seqs)
    out = defaultdict(list)
    for line, r in enumerate(records, 2):
        acc = r["Target Uniprot ID"].strip()
        if acc not in proteins:
            continue
        key = r["SBP ID"].strip()
        scaffold = scaffold_names.get(key, "")
        if scaffold.casefold() not in supported:
            continue
        seq = seqs.get(key, "")
        if not set(seq) <= AMINO_ACIDS:
            seq = ""
        reported = re.fullmatch(
            r"\s*(Kd|Ki|IC50)\s*[:=]\s*([<>=~≤≥≈]*)\s*(\d*\.?\d+(?:e[+-]?\d+)?)\s*([pnuµμm]?M)\s*",
            r["Affinity"],
            flags=re.I,
        )
        quantity = measurement(reported[1], reported[3], reported[4], reported[2]) if reported else {}
        out[acc].append(
            dict(
                source="SYNBIP",
                record_id=f"{key}:{r['TargetID']}:{line}",
                partner_id=key,
                partner_name=r["SBP Name"],
                partner_type="engineered_protein",
                partner_sequences=seq,
                sequence_basis="reported_SBP" if seq else "sequence_unavailable",
                status="parent_only",
                reason="construct_unspecified",
                query_role="target",
                biological_role=scaffold,
                affinity=r["Affinity"],
                citation="https://synbip.idrblab.net/data/sbp/details/" + key.lower(),
                **quantity,
            )
        )
    return out


def intact(cache, accession, proteins):
    expr = urllib.parse.quote(f"identifier:{accession} AND taxidA:9606 AND taxidB:9606", safe="")
    first, out = 0, []
    while True:
        url = f"{INTACT}{expr}?format=tab27&firstResult={first}&maxResults={PER_PAGE}"

        def parse(data):
            rows = [
                line.split("\t") for line in data.decode().splitlines() if line and not line.startswith("#")
            ]
            if any(len(r) < 38 for r in rows):
                raise ValueError("Expected PSI-MITAB 2.7 interaction records")
            return rows

        rows = cache.read(f"intact/{accession}.{first}.tsv", url, parse)
        for r in rows:
            ids = [r[i].removeprefix("uniprotkb:") for i in (0, 1)]
            if accession not in ids or not all("taxid:9606(" in r[i] for i in (9, 10)):
                continue
            partner = reference_partner([ids[1 - ids.index(accession)]], proteins)
            if not partner:
                continue
            side = ids.index(accession)
            direct = "MI:0407" in r[11] and r[15] in ("-", "")
            negative = r[35].lower() == "true"
            out.append(
                dict(
                    **partner,
                    source="IntAct",
                    record_id=r[13],
                    status="reported_negative" if negative else "parent_only",
                    reason="negative_assay"
                    if negative
                    else ("direct_interaction" if direct else "association_or_complex"),
                    method=r[6],
                    reference=r[8],
                    assay_parameters=r[29] if r[29] != "-" else "",
                    query_feature=r[36 + side] if r[36 + side] != "-" else "",
                    partner_feature=r[37 - side] if r[37 - side] != "-" else "",
                    citation="https://www.ebi.ac.uk/intact/search?query=" + accession,
                )
            )
        if len(rows) < PER_PAGE:
            return out
        first += len(rows)


def gtopdb(cache, proteins):
    key = os.environ.get("GTP_API_KEY")
    versions = set()

    def read(filename, fields):
        def parse(data):
            first = data.split(b"\n", 1)[0]
            banner = first.decode("utf-8-sig").rstrip("\r").strip('"')
            if not banner.startswith("# GtoPdb Version:"):
                raise ValueError(
                    "Expected a versioned GtoPdb export; authenticated downloads may need GTP_API_KEY"
                )
            records = csv_records(data, fields)
            versions.add(banner)
            return records

        url = GTOP + (f"gtopdb-downloads?file={filename}" if key else f"DATA/{filename}")
        return cache.read("gtopdb/" + filename, url, parse, headers={"GTP-API-Key": key} if key else {})

    peptides = {
        r["Ligand id"]: r
        for r in read(
            "peptides.csv",
            ["Ligand id", "Single letter amino acid sequence", "Species", "Subunit ids", "UniProt id"],
        )
    }
    pairs = read(
        "endogenous_ligands.csv", ["Ligand ID", "Ligand UniProt ID", "Target UniProt ID", "Target ID"]
    )
    measured = read(
        "interactions.csv", ["Ligand ID", "Ligand Type", "Target Species", "Target UniProt ID", "PubMed ID"]
    )
    if len(versions) != 1:
        raise ValueError("GtoPdb exports belong to different releases")
    out = defaultdict(list)
    endo = {(r["Ligand ID"], r["Target ID"]) for r in pairs}

    def add(ligand, target_ids, ligand_ids, evidence, record):
        peptide = peptides.get(ligand)
        target = reference_partner(target_ids.split("|"), proteins)
        if not peptide or not target:
            return
        parts = (
            [peptides.get(p) for p in peptide["Subunit ids"].split("|")]
            if peptide["Subunit ids"]
            else [peptide]
        )
        if any(p is None for p in parts):
            return
        if peptide["Species"] != "Human" and not all(p["Species"] == "Human" for p in parts):
            return
        pids = ligand_ids or peptide["UniProt id"] or peptide.get("Subunit UniProt IDs", "")
        ids = sorted(set(pids.split("|")) - {""})
        if any(a not in proteins for a in ids):
            return
        seqs = ["".join(p["Single letter amino acid sequence"].split()).upper() for p in parts]
        has_seq = all(s and set(s) <= AMINO_ACIDS for s in seqs)
        if not ids and not has_seq:
            return
        partner = dict(
            partner_id="|".join(ids) if ids else f"GtoPdb:{ligand}",
            partner_name=peptide["Name"],
            partner_type="protein_complex" if len(parts) > 1 else "protein_or_peptide",
            partner_sequences="|".join(seqs) if has_seq else "",
            sequence_basis="GtoPdb_reported_chains" if has_seq else "sequence_unavailable",
        )
        target_subunits = evidence.pop("target_subunits", "")
        context = " | ".join(x for x in (peptide.get("Subunit names", ""), target_subunits) if x)
        common = dict(
            source="GtoPdb",
            record_id=record,
            status="parent_only",
            complex_context=context,
            modifications=" | ".join(
                filter(
                    None,
                    (peptide.get("Post-translational modification"), peptide.get("Chemical modification")),
                )
            ),
            citation=GTOP + "GRAC/LigandDisplayForward?ligandId=" + ligand,
            **evidence,
        )
        for a in target_ids.split("|"):
            out[a].append(
                dict(
                    **partner,
                    **common,
                    query_role="target",
                    reason="query_complex_subunit"
                    if target_subunits or "|" in target_ids
                    else "construct_unspecified",
                )
            )
        for a in ids:
            out[a].append(
                dict(
                    **target,
                    **common,
                    query_role="ligand",
                    reason="query_complex_subunit" if len(ids) > 1 else "construct_unspecified",
                )
            )

    for i, r in enumerate(pairs, 3):
        add(
            r["Ligand ID"],
            r["Target UniProt ID"],
            r["Ligand UniProt ID"],
            {"biological_role": "endogenous_pair"},
            f"endogenous:{i}",
        )
    for i, r in enumerate(measured, 3):
        if r["Ligand Type"] != "Peptide" or r["Target Species"] != "Human":
            continue
        quantitative = {k: v for k, v in r.items() if "Affinity" in k and v not in ("", "-")}
        has_value = any(
            re.search(r"\d", v) for k, v in quantitative.items() if "Units" not in k and "Relation" not in k
        )
        quantity = measurement(
            r.get("Original Affinity Units", ""),
            r.get("Original Affinity Median nm", ""),
            "nM",
            r.get("Original Affinity Relation", ""),
        )
        role = "endogenous_pair" if (r["Ligand ID"], r["Target ID"]) in endo else "reported_peptide"
        target_ids = r["Target UniProt ID"]
        if not target_ids and r.get("Target Ligand ID") in peptides:
            target_ids = r.get("Target Ligand UniProt ID", "")
        add(
            r["Ligand ID"],
            target_ids,
            "",
            dict(
                biological_role=role,
                method=r.get("Assay Description", ""),
                reference=r["PubMed ID"],
                target_subunits=r.get("Target Subunit IDs", ""),
                affinity=json.dumps(quantitative, sort_keys=True) if has_value else "",
                **quantity,
            ),
            f"interactions:{i}",
        )
    return out
