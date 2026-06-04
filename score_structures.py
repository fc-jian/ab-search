#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
import pandas as pd
from anarci import anarci
from Bio.Align import PairwiseAligner, substitution_matrices
from biotite.sequence import ProteinSequence


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovs",
]

MISSING_VALUES = {"", "nan", "none", "na", "n/a", "null", "-", "?"}
CHAIN_SPLIT_RE = re.compile(r"\s*\|\s*|[,;]\s*")
ERROR_COLUMNS = ["query_id", "pdb_id", "antigen_id", "error"]
CDR_DEFINITIONS = {
    "imgt": {
        "default": {
            "CDR1": (27, 38),
            "CDR2": (56, 65),
            "CDR3": (105, 117),
        },
    },
    "chothia": {
        "H": {
            "CDR1": (26, 32),
            "CDR2": (52, 56),
            "CDR3": (95, 102),
        },
        "L": {
            "CDR1": (24, 34),
            "CDR2": (50, 56),
            "CDR3": (89, 97),
        },
        "default": {
            "CDR1": (24, 34),
            "CDR2": (50, 56),
            "CDR3": (89, 97),
        },
    },
}


@dataclass(frozen=True)
class ResidueRecord:
    chain_id: str
    res_id: int
    ins_code: str
    res_name: str
    aa: str
    seq_index: int
    atoms: Any


@dataclass(frozen=True)
class EpitopeContact:
    residue: ResidueRecord
    ca_distance: float
    atom_distance: float
    chain_role: str


@dataclass(frozen=True)
class EpitopeResult:
    contacts: tuple[EpitopeContact, ...]
    multichain_chain_ids: tuple[str, ...]
    homodimer_excluded_chain_ids: tuple[str, ...]
    cdr_chain_ids: tuple[str, ...]
    cdr_warnings: tuple[str, ...]


@dataclass(frozen=True)
class AnarciDomain:
    chain_type: str
    sequence: str
    cdr1: str
    cdr2: str
    cdr3: str
    cdr_indices: tuple[int, ...]
    start: int | None
    end: int | None


@dataclass(frozen=True)
class StructureVariant:
    atoms: Any
    source: str
    warning: str


def clean_cell(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in MISSING_VALUES:
        return ""
    return text


def split_chains(value: object) -> list[str]:
    text = clean_cell(value)
    if not text:
        return []
    chains = []
    for part in CHAIN_SPLIT_RE.split(text):
        chain = part.strip()
        if chain and chain.lower() not in MISSING_VALUES:
            chains.append(chain)
    return chains


def unique_preserve_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value and value not in seen:
            result.append(value)
            seen.add(value)
    return result


def is_true(value: object) -> bool:
    return clean_cell(value).lower() in {"true", "t", "yes", "y", "1"}


def normalize_identity_threshold(value: float) -> float:
    return value * 100.0 if value <= 1.0 else value


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper().replace("*", "")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        records[name] = "".join(chunks).upper().replace("*", "")

    return records


def open_fasta_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_seqres_header(header: str) -> tuple[str, str]:
    record_id = header.split()[0]
    if "_" not in record_id:
        return "", ""
    pdb_id, chain_id = record_id.split("_", 1)
    return pdb_id.lower(), chain_id


def iter_fasta(path: Path) -> Iterable[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []

    with open_fasta_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)

    if header is not None:
        yield header, "".join(chunks)


def load_seqres_sequences(
    path: Path | None,
    keys: set[tuple[str, str]],
) -> dict[tuple[str, str], str]:
    if path is None or not path.exists() or not keys:
        return {}

    sequences: dict[tuple[str, str], str] = {}
    remaining = set(keys)
    for header, sequence in iter_fasta(path):
        key = parse_seqres_header(header)
        if key in remaining:
            sequences[key] = sequence.upper().replace("*", "")
            remaining.remove(key)
            if not remaining:
                break
    return sequences


def read_blast_hits(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    with path.open() as handle:
        first_line = handle.readline()

    column_count = len(first_line.rstrip("\n").split("\t"))
    columns = BLAST_COLUMNS[:column_count]
    if column_count > len(BLAST_COLUMNS):
        columns = BLAST_COLUMNS + [f"extra_{i}" for i in range(column_count - len(BLAST_COLUMNS))]

    hits = pd.read_csv(path, sep="\t", header=None, names=columns, dtype=str)
    for column in ("pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"):
        if column in hits:
            hits[column] = pd.to_numeric(hits[column], errors="coerce")

    if not hits.empty:
        parsed = hits["sseqid"].apply(parse_subject_id)
        hits["pdb_id"] = [item[0] for item in parsed]
        hits["antigen_id"] = [item[1] for item in parsed]

    return hits


def parse_subject_id(subject_id: object) -> tuple[str, str]:
    text = str(subject_id).strip()
    if "|" in text:
        parts = text.split("|")
        if len(parts) >= 3 and len(parts[1]) == 4:
            return parts[1].lower(), parts[2]

    record = text.split()[0]
    if "_" in record:
        pdb_id, chain_id = record.split("_", 1)
        return pdb_id[:4].lower(), chain_id

    return record[:4].lower(), ""


def load_pdbx_file(path: Path) -> Any:
    name = path.name.lower()
    if name.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            return pdbx.CIFFile.read(handle)
    elif name.endswith((".cif", ".mmcif", ".pdbx")):
        return pdbx.CIFFile.read(path)
    elif name.endswith(".bcif"):
        return pdbx.BinaryCIFFile.read(path)
    else:
        raise ValueError(f"unsupported structure format: {path}")


def get_asymmetric_unit(pdbx_file: Any, model: int = 1) -> Any:
    return pdbx.get_structure(
        pdbx_file,
        model=model,
        altloc="occupancy",
        use_author_fields=True,
        include_bonds=False,
    )


def get_first_assembly(pdbx_file: Any, model: int = 1) -> tuple[Any, str]:
    assemblies = pdbx.list_assemblies(pdbx_file)
    if not assemblies:
        raise ValueError("no biological assemblies listed in mmCIF")
    assembly_id = next(iter(assemblies))
    atoms = pdbx.get_assembly(
        pdbx_file,
        assembly_id=assembly_id,
        model=model,
        altloc="occupancy",
        use_author_fields=True,
        include_bonds=False,
    )
    return atoms, str(assembly_id)


def load_structure_variants(path: Path, model: int = 1) -> list[StructureVariant]:
    pdbx_file = load_pdbx_file(path)
    variants: list[StructureVariant] = []
    assembly_warning = ""

    try:
        assembly_atoms, assembly_id = get_first_assembly(pdbx_file, model=model)
        variants.append(
            StructureVariant(
                atoms=assembly_atoms,
                source=f"assembly:{assembly_id}",
                warning="",
            )
        )
    except Exception as exc:
        assembly_warning = f"first assembly unavailable: {exc}"

    asym_warning = assembly_warning
    if variants:
        asym_warning = "fallback to asymmetric unit after first assembly had no matching CDR ab-ag contacts"

    variants.append(
        StructureVariant(
            atoms=get_asymmetric_unit(pdbx_file, model=model),
            source="asymmetric_unit",
            warning=asym_warning,
        )
    )
    return variants


def find_structure_file(cache_dir: Path, pdb_id: str) -> Path | None:
    candidates = [
        cache_dir / f"{pdb_id.lower()}.cif.gz",
        cache_dir / f"{pdb_id.upper()}.cif.gz",
        cache_dir / f"{pdb_id.lower()}.cif",
        cache_dir / f"{pdb_id.upper()}.cif",
        cache_dir / f"{pdb_id.lower()}.mmcif",
        cache_dir / f"{pdb_id.upper()}.mmcif",
    ]
    for path in candidates:
        if path.exists() and path.stat().st_size > 0:
            return path
    return None


def residue_aa(res_name: str) -> str:
    try:
        return str(ProteinSequence.convert_letter_3to1(res_name.upper()))
    except Exception:
        return "X"


def chain_residues(atoms: Any, chain_id: str) -> list[ResidueRecord]:
    protein = atoms[struc.filter_amino_acids(atoms)]
    chain_atoms = protein[protein.chain_id == chain_id]
    if chain_atoms.array_length() == 0:
        return []

    starts = struc.get_residue_starts(chain_atoms, add_exclusive_stop=True)
    residues: list[ResidueRecord] = []
    for seq_index, (start, stop) in enumerate(zip(starts[:-1], starts[1:])):
        res_atoms = chain_atoms[start:stop]
        atom0 = res_atoms[0]
        ins_code = clean_cell(getattr(atom0, "ins_code", ""))
        res_name = str(atom0.res_name)
        residues.append(
            ResidueRecord(
                chain_id=str(atom0.chain_id),
                res_id=int(atom0.res_id),
                ins_code=ins_code,
                res_name=res_name,
                aa=residue_aa(res_name),
                seq_index=seq_index,
                atoms=res_atoms,
            )
        )
    return residues


def chain_sequence(atoms: Any, chain_id: str) -> str:
    return "".join(residue.aa for residue in chain_residues(atoms, chain_id))


def cdr_indices_from_domains(domains: list[AnarciDomain]) -> set[int]:
    indices: set[int] = set()
    for domain in domains:
        indices.update(domain.cdr_indices)
    return indices


def observed_cdr_indices(
    observed_sequence: str,
    reference_sequence: str,
    reference_cdr_indices: set[int],
    aligner: PairwiseAligner,
) -> set[int]:
    mapping = subject_to_query_mapping(reference_sequence, observed_sequence, aligner)
    return {
        observed_index
        for observed_index, reference_index in mapping.items()
        if reference_index in reference_cdr_indices
    }


def cdr_residue_indices_for_chain(
    atoms: Any,
    pdb_id: str,
    chain_id: str,
    scheme: str,
    seqres_sequences: dict[tuple[str, str], str],
    aligner: PairwiseAligner,
) -> tuple[set[int], str]:
    observed_sequence = chain_sequence(atoms, chain_id)
    if not observed_sequence:
        return set(), f"{chain_id}: chain sequence not found"

    reference_sequence = seqres_sequences.get((pdb_id, chain_id))
    if reference_sequence:
        domains, error = run_anarci_with_scfv_retry(
            f"{pdb_id}_{chain_id}_seqres",
            reference_sequence,
            scheme,
        )
        reference_cdr_indices = cdr_indices_from_domains(domains)
        if not reference_cdr_indices:
            return set(), f"{chain_id}: no CDR residues from SEQRES ANARCI ({error})"
        indices = observed_cdr_indices(
            observed_sequence,
            reference_sequence,
            reference_cdr_indices,
            aligner,
        )
        if not indices:
            return set(), f"{chain_id}: SEQRES CDR residues not present in structure"
        return indices, ""

    domains, error = run_anarci_with_scfv_retry(
        f"{pdb_id}_{chain_id}_structure",
        observed_sequence,
        scheme,
    )
    indices = cdr_indices_from_domains(domains)
    if not indices:
        return set(), f"{chain_id}: no CDR residues from structure ANARCI ({error})"
    return indices, ""


def cdr_atoms_for_antibody_chains(
    atoms: Any,
    pdb_id: str,
    antibody_chains: list[str],
    scheme: str,
    seqres_sequences: dict[tuple[str, str], str],
    aligner: PairwiseAligner,
) -> tuple[Any | None, list[str], list[str]]:
    residue_arrays: list[Any] = []
    cdr_chain_ids: list[str] = []
    warnings: list[str] = []

    for chain_id in antibody_chains:
        cdr_indices, warning = cdr_residue_indices_for_chain(
            atoms,
            pdb_id,
            chain_id,
            scheme,
            seqres_sequences,
            aligner,
        )
        if warning:
            warnings.append(warning)
        if not cdr_indices:
            continue

        selected = [
            residue.atoms
            for residue in chain_residues(atoms, chain_id)
            if residue.seq_index in cdr_indices
        ]
        if selected:
            residue_arrays.extend(selected)
            cdr_chain_ids.append(chain_id)

    if not residue_arrays:
        return None, [], unique_preserve_order(warnings)

    return struc.concatenate(residue_arrays), unique_preserve_order(cdr_chain_ids), unique_preserve_order(warnings)


def residue_label(residue: ResidueRecord) -> str:
    return f"{residue.chain_id}:{residue.res_id}{residue.ins_code}"


def min_distance(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    if coords_a.size == 0 or coords_b.size == 0:
        return math.nan
    distances = struc.distance(coords_a[:, None, :], coords_b[None, :, :])
    return float(np.min(distances))


def chain_ids_in_structure(atoms: Any) -> list[str]:
    protein = atoms[struc.filter_amino_acids(atoms)]
    return unique_preserve_order(str(chain_id) for chain_id in protein.chain_id)


def sequence_identity(seq_a: str, seq_b: str, aligner: PairwiseAligner) -> float:
    if not seq_a or not seq_b:
        return 0.0

    alignment = aligner.align(seq_a, seq_b)[0]
    blocks_a, blocks_b = alignment.aligned
    matches = 0
    for block_a, block_b in zip(blocks_a, blocks_b):
        start_a, end_a = map(int, block_a)
        start_b, end_b = map(int, block_b)
        span = min(end_a - start_a, end_b - start_b)
        for offset in range(span):
            if seq_a[start_a + offset] == seq_b[start_b + offset]:
                matches += 1

    return matches / max(len(seq_a), len(seq_b))


def contact_residues_for_chain(
    atoms: Any,
    chain_id: str,
    antibody_heavy: Any,
    antibody_ca: Any,
    ca_distance_threshold: float,
    atom_distance_threshold: float,
    chain_role: str,
) -> list[EpitopeContact]:
    contacts: list[EpitopeContact] = []

    for residue in chain_residues(atoms, chain_id):
        residue_ca = residue.atoms[residue.atoms.atom_name == "CA"]
        ca_min = min_distance(residue_ca.coord, antibody_ca.coord) if residue_ca.array_length() else math.nan

        residue_heavy = residue.atoms[residue.atoms.element != "H"]
        atom_min = (
            min_distance(residue_heavy.coord, antibody_heavy.coord)
            if residue_heavy.array_length()
            else math.nan
        )

        ca_contact = not math.isnan(ca_min) and ca_min <= ca_distance_threshold
        atom_contact = not math.isnan(atom_min) and atom_min <= atom_distance_threshold
        if ca_contact or atom_contact:
            contacts.append(
                EpitopeContact(
                    residue=residue,
                    ca_distance=ca_min,
                    atom_distance=atom_min,
                    chain_role=chain_role,
                )
            )

    return contacts


def identify_epitope_residues(
    atoms: Any,
    pdb_id: str,
    antigen_chain: str,
    antibody_chains: list[str],
    excluded_chains: list[str],
    scheme: str,
    seqres_sequences: dict[tuple[str, str], str],
    ca_distance_threshold: float,
    atom_distance_threshold: float,
    homodimer_identity_threshold: float,
    aligner: PairwiseAligner,
) -> EpitopeResult:
    if not chain_residues(atoms, antigen_chain):
        return EpitopeResult(
            contacts=(),
            multichain_chain_ids=(),
            homodimer_excluded_chain_ids=(),
            cdr_chain_ids=(),
            cdr_warnings=(f"{antigen_chain}: antigen chain not found",),
        )

    cdr_atoms, cdr_chain_ids, cdr_warnings = cdr_atoms_for_antibody_chains(
        atoms,
        pdb_id,
        antibody_chains,
        scheme,
        seqres_sequences,
        aligner,
    )
    if cdr_atoms is None or cdr_atoms.array_length() == 0:
        return EpitopeResult(
            contacts=(),
            multichain_chain_ids=(),
            homodimer_excluded_chain_ids=(),
            cdr_chain_ids=(),
            cdr_warnings=tuple(cdr_warnings or ["no antibody CDR atoms found"]),
        )

    cdr_heavy = cdr_atoms[cdr_atoms.element != "H"]
    cdr_ca = cdr_atoms[cdr_atoms.atom_name == "CA"]
    contacts = contact_residues_for_chain(
        atoms,
        antigen_chain,
        cdr_heavy,
        cdr_ca,
        ca_distance_threshold,
        atom_distance_threshold,
        chain_role="target",
    )

    target_sequence = chain_sequence(atoms, antigen_chain)
    excluded = set(excluded_chains) | set(antibody_chains) | {antigen_chain}
    multichain_chain_ids: list[str] = []
    homodimer_excluded_chain_ids: list[str] = []

    for chain_id in chain_ids_in_structure(atoms):
        if chain_id in excluded:
            continue

        other_contacts = contact_residues_for_chain(
            atoms,
            chain_id,
            cdr_heavy,
            cdr_ca,
            ca_distance_threshold,
            atom_distance_threshold,
            chain_role="multichain",
        )
        if not other_contacts:
            continue

        other_sequence = chain_sequence(atoms, chain_id)
        identity = sequence_identity(target_sequence, other_sequence, aligner)
        if identity >= homodimer_identity_threshold:
            homodimer_excluded_chain_ids.append(chain_id)
            continue

        contacts.extend(other_contacts)
        multichain_chain_ids.append(chain_id)

    return EpitopeResult(
        contacts=tuple(contacts),
        multichain_chain_ids=tuple(unique_preserve_order(multichain_chain_ids)),
        homodimer_excluded_chain_ids=tuple(unique_preserve_order(homodimer_excluded_chain_ids)),
        cdr_chain_ids=tuple(cdr_chain_ids),
        cdr_warnings=tuple(cdr_warnings),
    )


def make_pairwise_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def subject_to_query_mapping(
    query_sequence: str,
    subject_sequence: str,
    aligner: PairwiseAligner,
) -> dict[int, int]:
    if not query_sequence or not subject_sequence:
        return {}

    alignment = aligner.align(query_sequence, subject_sequence)[0]
    query_blocks, subject_blocks = alignment.aligned
    mapping: dict[int, int] = {}
    for query_block, subject_block in zip(query_blocks, subject_blocks):
        query_start, query_end = map(int, query_block)
        subject_start, subject_end = map(int, subject_block)
        span = min(query_end - query_start, subject_end - subject_start)
        for offset in range(span):
            mapping[subject_start + offset] = query_start + offset
    return mapping


def chain_family_from_type(chain_type: str) -> str:
    text = chain_type.upper()
    if text == "H" or "HEAVY" in text or text.startswith("VH"):
        return "H"
    if text in {"K", "L"} or "KAPPA" in text or "LAMBDA" in text or text.startswith(("VK", "VL")):
        return "L"
    return "default"


def cdr_definition(scheme: str, chain_type: str) -> dict[str, tuple[int, int]]:
    scheme_defs = CDR_DEFINITIONS.get(scheme.lower())
    if scheme_defs is None:
        scheme_defs = CDR_DEFINITIONS["chothia"]
    family = chain_family_from_type(chain_type)
    return scheme_defs.get(family) or scheme_defs["default"]


def extract_region(numbered: list[tuple[tuple[int, str], str]], start: int, end: int) -> str:
    residues = []
    for (position_number, _insertion_code), aa in numbered:
        if start <= position_number <= end and aa != "-":
            residues.append(aa)
    return "".join(residues)


def extract_cdr_sequences(
    numbered: list[tuple[tuple[int, str], str]],
    scheme: str,
    chain_type: str,
) -> dict[str, str]:
    regions = cdr_definition(scheme, chain_type)
    return {
        region: extract_region(numbered, start, end)
        for region, (start, end) in regions.items()
    }


def extract_cdr_indices(
    numbered: list[tuple[tuple[int, str], str]],
    domain_start: int | None,
    scheme: str,
    chain_type: str,
) -> tuple[int, ...]:
    regions = cdr_definition(scheme, chain_type)
    sequence_index = int(domain_start or 0)
    cdr_indices: list[int] = []

    for (position_number, _insertion_code), aa in numbered:
        if aa == "-":
            continue
        if any(start <= position_number <= end for start, end in regions.values()):
            cdr_indices.append(sequence_index)
        sequence_index += 1

    return tuple(cdr_indices)


def parse_chain_type(detail: object) -> str:
    if isinstance(detail, dict):
        for key in ("chain_type", "chain", "query_type"):
            value = clean_cell(detail.get(key))
            if value:
                return value
    return ""


def is_heavy_domain(domain: AnarciDomain) -> bool:
    return chain_family_from_type(domain.chain_type) == "H"


def is_light_domain(domain: AnarciDomain) -> bool:
    return chain_family_from_type(domain.chain_type) == "L"


def run_anarci_domains(sequence_id: str, sequence: str, scheme: str) -> tuple[list[AnarciDomain], str]:
    if len(sequence) < 50:
        return [], "sequence shorter than antibody variable domain"

    try:
        numbering, alignment_details, _hit_tables = anarci(
            [(sequence_id, sequence)],
            scheme=scheme,
            output=False,
        )
    except Exception as exc:
        return [], f"ANARCI failed: {exc}"

    if not numbering or numbering[0] is None:
        return [], "ANARCI found no antibody domains"

    seq_details = alignment_details[0] if alignment_details else None
    domains: list[AnarciDomain] = []

    for index, domain_result in enumerate(numbering[0]):
        domain_numbering, domain_start, domain_end = domain_result
        detail = None
        if isinstance(seq_details, list) and index < len(seq_details):
            detail = seq_details[index]
        elif isinstance(seq_details, dict):
            detail = seq_details

        chain_type = parse_chain_type(detail)
        cdr_sequences = extract_cdr_sequences(domain_numbering, scheme, chain_type)
        sequence_numbered = "".join(aa for _position, aa in domain_numbering if aa != "-")
        domains.append(
            AnarciDomain(
                chain_type=chain_type,
                sequence=sequence_numbered,
                cdr1=cdr_sequences["CDR1"],
                cdr2=cdr_sequences["CDR2"],
                cdr3=cdr_sequences["CDR3"],
                cdr_indices=extract_cdr_indices(
                    domain_numbering,
                    domain_start,
                    scheme,
                    chain_type,
                ),
                start=int(domain_start) if domain_start is not None else None,
                end=int(domain_end) if domain_end is not None else None,
            )
        )

    return domains, ""


def run_anarci_with_scfv_retry(
    sequence_id: str,
    sequence: str,
    scheme: str,
) -> tuple[list[AnarciDomain], str]:
    domains, error = run_anarci_domains(sequence_id, sequence, scheme)
    if len(domains) != 1 or domains[0].end is None:
        return domains, error

    retry_errors = []
    for cut in (domains[0].end + 1, domains[0].end, 110, 120):
        if cut <= 0 or cut >= len(sequence) - 50:
            continue
        tail_domains, tail_error = run_anarci_domains(f"{sequence_id}_tail_{cut}", sequence[cut:], scheme)
        if tail_domains:
            return domains + tail_domains, ""
        if tail_error:
            retry_errors.append(tail_error)

    return domains, error or "; ".join(retry_errors)


def choose_domain(domains: list[AnarciDomain], kind: str) -> AnarciDomain | None:
    if kind == "heavy":
        for domain in domains:
            if is_heavy_domain(domain):
                return domain
    elif kind == "light":
        for domain in domains:
            if is_light_domain(domain):
                return domain

    return domains[0] if domains else None


def antibody_chain_ids(row: pd.Series, antigen_chain: str) -> tuple[list[str], list[str], list[str]]:
    heavy_ids = split_chains(row.get("Hchain"))
    light_ids = split_chains(row.get("Lchain"))
    contact_ids = unique_preserve_order(
        chain for chain in heavy_ids + light_ids if chain and chain != antigen_chain
    )
    return heavy_ids, light_ids, contact_ids


def annotate_antibodies(
    atoms: Any,
    row: pd.Series,
    scheme: str,
    annotation_cache: dict[tuple[str, str], tuple[str, list[AnarciDomain], str]],
    seqres_sequences: dict[tuple[str, str], str],
) -> dict[str, str]:
    scfv = is_true(row.get("scfv"))
    heavy_ids = split_chains(row.get("Hchain"))
    light_ids = split_chains(row.get("Lchain"))
    all_ids = unique_preserve_order(heavy_ids + light_ids)

    chain_data: dict[str, tuple[str, list[AnarciDomain], str]] = {}
    pdb_id = clean_cell(row.get("pdb")).lower()
    for chain_id in all_ids:
        cache_key = (pdb_id, chain_id)
        if cache_key not in annotation_cache:
            sequence = seqres_sequences.get(cache_key) or chain_sequence(atoms, chain_id)
            domains, error = run_anarci_with_scfv_retry(f"{pdb_id}_{chain_id}", sequence, scheme)
            annotation_cache[cache_key] = (sequence, domains, error)
        chain_data[chain_id] = annotation_cache[cache_key]

    heavy_sequences = [chain_data[chain_id][0] for chain_id in heavy_ids if chain_id in chain_data]
    light_sequences = [chain_data[chain_id][0] for chain_id in light_ids if chain_id in chain_data]

    heavy_domains = [
        domain
        for chain_id in heavy_ids
        if chain_id in chain_data
        for domain in chain_data[chain_id][1]
    ]
    light_domains = [
        domain
        for chain_id in light_ids
        if chain_id in chain_data
        for domain in chain_data[chain_id][1]
    ]

    if scfv or set(heavy_ids) & set(light_ids):
        combined_domains = [
            domain
            for chain_id in all_ids
            if chain_id in chain_data
            for domain in chain_data[chain_id][1]
        ]
        heavy_domain = choose_domain(combined_domains, "heavy")
        light_domain = choose_domain(
            [domain for domain in combined_domains if domain is not heavy_domain],
            "light",
        )
    else:
        heavy_domain = choose_domain(heavy_domains, "heavy")
        light_domain = choose_domain(light_domains, "light")

    errors = [
        error
        for _sequence, domains, error in chain_data.values()
        if error and not domains
    ]

    return {
        "heavy_chain_id": ";".join(heavy_ids),
        "light_chain_id": "" if scfv else ";".join(light_ids),
        "seq_heavy_full": ";".join(heavy_sequences),
        "seq_light_full": "" if scfv else ";".join(light_sequences),
        "seq_heavy_anarci": heavy_domain.sequence if heavy_domain else "",
        "seq_light_anarci": light_domain.sequence if light_domain else "",
        "CDR-H1": heavy_domain.cdr1 if heavy_domain else "",
        "CDR-H2": heavy_domain.cdr2 if heavy_domain else "",
        "CDR-H3": heavy_domain.cdr3 if heavy_domain else "",
        "CDR-L1": light_domain.cdr1 if light_domain else "",
        "CDR-L2": light_domain.cdr2 if light_domain else "",
        "CDR-L3": light_domain.cdr3 if light_domain else "",
        "anarci_warnings": "; ".join(unique_preserve_order(errors)),
    }


def matching_sabdab_rows(summary: pd.DataFrame, pdb_id: str, antigen_chain: str) -> list[pd.Series]:
    subset = summary[summary["pdb_norm"] == pdb_id.lower()]
    matches: list[pd.Series] = []
    for _, row in subset.iterrows():
        antigen_type = clean_cell(row.get("antigen_type")).lower()
        if "protein" not in antigen_type:
            continue
        if antigen_chain in split_chains(row.get("antigen_chain")):
            matches.append(row)
    return matches


def all_antibody_chains_for_pdb(summary: pd.DataFrame, pdb_id: str) -> list[str]:
    subset = summary[summary["pdb_norm"] == pdb_id.lower()]
    chains: list[str] = []
    for _, row in subset.iterrows():
        chains.extend(split_chains(row.get("Hchain")))
        chains.extend(split_chains(row.get("Lchain")))
    return unique_preserve_order(chains)


def score_epitope(
    query_sequence: str,
    antigen_residues: list[ResidueRecord],
    epitope: EpitopeResult,
    aligner: PairwiseAligner,
) -> dict[str, object]:
    subject_sequence = "".join(residue.aa for residue in antigen_residues)
    mapping = subject_to_query_mapping(query_sequence, subject_sequence, aligner)

    identical = 0
    mutated = 0
    covered = 0
    mutation_labels: list[str] = []
    target_contacts = [contact for contact in epitope.contacts if contact.chain_role == "target"]
    multichain_contacts = [contact for contact in epitope.contacts if contact.chain_role == "multichain"]

    for contact in epitope.contacts:
        residue = contact.residue
        if contact.chain_role == "multichain":
            mutated += 1
            mutation_labels.append(f"{residue_label(residue)}:{residue.aa}>missing_chain")
            continue

        query_index = mapping.get(residue.seq_index)
        if query_index is None:
            mutated += 1
            mutation_labels.append(f"{residue_label(residue)}:{residue.aa}>-")
            continue

        covered += 1
        query_aa = query_sequence[query_index]
        if query_aa == residue.aa:
            identical += 1
        else:
            mutated += 1
            mutation_labels.append(f"{residue_label(residue)}:{residue.aa}>{query_aa}")

    epitope_count = len(epitope.contacts)
    return {
        "number_epitope_residues": epitope_count,
        "number_target_epitope_residues": len(target_contacts),
        "number_multichain_epitope_residues": len(multichain_contacts),
        "number_identical_epitope_residues": identical,
        "number_mutated_epitope_residues": mutated,
        "epitope_identity": identical / epitope_count if epitope_count else math.nan,
        "epitope_covered_residues": covered,
        "epitope_residue_ids": ";".join(residue_label(contact.residue) for contact in epitope.contacts),
        "target_epitope_residue_ids": ";".join(residue_label(contact.residue) for contact in target_contacts),
        "multichain_epitope_residue_ids": ";".join(
            residue_label(contact.residue) for contact in multichain_contacts
        ),
        "multichain_epitope_chain_ids": "|".join(epitope.multichain_chain_ids),
        "homodimer_excluded_chain_ids": "|".join(epitope.homodimer_excluded_chain_ids),
        "cdr_contact_chain_ids": "|".join(epitope.cdr_chain_ids),
        "cdr_contact_warnings": "; ".join(epitope.cdr_warnings),
        "mutated_epitope_residues": ";".join(mutation_labels),
        "subject_structure_sequence": subject_sequence,
    }


def base_output_columns() -> list[str]:
    return [
        "query_id",
        "pdb_id",
        "antigen_id",
        "heavy_chain_id",
        "light_chain_id",
        "number_epitope_residues",
        "number_target_epitope_residues",
        "number_multichain_epitope_residues",
        "number_identical_epitope_residues",
        "number_mutated_epitope_residues",
        "epitope_identity",
        "seq_heavy_full",
        "seq_light_full",
        "seq_heavy_anarci",
        "seq_light_anarci",
        "CDR-H1",
        "CDR-H2",
        "CDR-H3",
        "CDR-L1",
        "CDR-L2",
        "CDR-L3",
        "pident",
        "qcovs",
        "evalue",
        "bitscore",
        "antigen_name",
        "epitope_covered_residues",
        "epitope_residue_ids",
        "target_epitope_residue_ids",
        "multichain_epitope_residue_ids",
        "multichain_epitope_chain_ids",
        "homodimer_excluded_chain_ids",
        "cdr_contact_chain_ids",
        "cdr_contact_warnings",
        "mutated_epitope_residues",
        "subject_structure_sequence",
        "structure_source",
        "structure_warnings",
        "anarci_warnings",
    ]


def merged_antibody_columns() -> list[str]:
    return [
        "antibody_group_id",
        "support_pdb_ids",
        "support_pdb_antigen_ids",
        "support_structure_count",
        "epitope_identity_range",
        "query_id_list",
        "heavy_chain_id_list",
        "light_chain_id_list",
        "seq_heavy_full",
        "seq_light_full",
        "seq_heavy_anarci",
        "seq_light_anarci",
        "CDR-H1",
        "CDR-H2",
        "CDR-H3",
        "CDR-L1",
        "CDR-L2",
        "CDR-L3",
        "antigen_name_list",
    ]


def join_unique(values: Iterable[object], separator: str = "|") -> str:
    cleaned = [clean_cell(value) for value in values]
    return separator.join(unique_preserve_order(value for value in cleaned if value))


def support_pdb_antigen_ids(group: pd.DataFrame) -> str:
    values = []
    for _, row in group.iterrows():
        pdb_id = clean_cell(row.get("pdb_id"))
        antigen_id = clean_cell(row.get("antigen_id"))
        if pdb_id and antigen_id:
            values.append(f"{pdb_id}:{antigen_id}")
    return "|".join(unique_preserve_order(values))


def format_float(value: float) -> str:
    return f"{value:.6g}"


def epitope_identity_range(group: pd.DataFrame) -> str:
    identities = pd.to_numeric(group["epitope_identity"], errors="coerce").dropna()
    if identities.empty:
        return ""

    min_identity = float(identities.min())
    max_identity = float(identities.max())
    if math.isclose(min_identity, max_identity):
        return format_float(min_identity)
    return f"{format_float(min_identity)}-{format_float(max_identity)}"


def make_merged_antibody_table(output: pd.DataFrame) -> pd.DataFrame:
    if output.empty:
        return pd.DataFrame(columns=merged_antibody_columns())

    work = output.copy()
    for column in ("seq_heavy_full", "seq_light_full"):
        work[column] = work[column].fillna("").astype(str)

    rows: list[dict[str, object]] = []
    for index, ((_heavy_full, _light_full), group) in enumerate(
        work.groupby(["seq_heavy_full", "seq_light_full"], sort=False, dropna=False),
        start=1,
    ):
        first = group.iloc[0]
        rows.append(
            {
                "antibody_group_id": f"AB{index:05d}",
                "support_pdb_ids": join_unique(group["pdb_id"]),
                "support_pdb_antigen_ids": support_pdb_antigen_ids(group),
                "support_structure_count": len(group),
                "epitope_identity_range": epitope_identity_range(group),
                "query_id_list": join_unique(group["query_id"]),
                "heavy_chain_id_list": join_unique(group["heavy_chain_id"]),
                "light_chain_id_list": join_unique(group["light_chain_id"]),
                "seq_heavy_full": clean_cell(first.get("seq_heavy_full")),
                "seq_light_full": clean_cell(first.get("seq_light_full")),
                "seq_heavy_anarci": clean_cell(first.get("seq_heavy_anarci")),
                "seq_light_anarci": clean_cell(first.get("seq_light_anarci")),
                "CDR-H1": clean_cell(first.get("CDR-H1")),
                "CDR-H2": clean_cell(first.get("CDR-H2")),
                "CDR-H3": clean_cell(first.get("CDR-H3")),
                "CDR-L1": clean_cell(first.get("CDR-L1")),
                "CDR-L2": clean_cell(first.get("CDR-L2")),
                "CDR-L3": clean_cell(first.get("CDR-L3")),
                "antigen_name_list": join_unique(group["antigen_name"]),
            }
        )

    return pd.DataFrame(rows, columns=merged_antibody_columns())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score SAbDab antibody-antigen structures against query sequence epitopes."
    )
    parser.add_argument("--query-fasta", required=True, type=Path)
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--sabdab", required=True, type=Path)
    parser.add_argument("--pdb-seqres", type=Path, default=None)
    parser.add_argument("--pdb-cache", "--pdb_cache", dest="pdb_cache", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--merged-antibody-output", type=Path, default=None)
    parser.add_argument("--errors-output", type=Path, default=None)
    parser.add_argument("--min-identity", type=float, default=0.7)
    parser.add_argument("--max-hits", type=int, default=200)
    parser.add_argument("--ca-distance-threshold", type=float, default=8.0)
    parser.add_argument("--atom-distance-threshold", type=float, default=4.5)
    parser.add_argument(
        "--homodimer-identity-threshold",
        type=float,
        default=0.95,
        help="Exclude other contacting chains as homodimers when sequence identity to target antigen is at least this value.",
    )
    parser.add_argument("--scheme", default="chothia")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    queries = read_fasta(args.query_fasta)
    hits = read_blast_hits(args.hits)
    summary = pd.read_csv(args.sabdab, sep="\t", dtype=str)
    summary["pdb_norm"] = summary["pdb"].str.lower()

    output_rows: list[dict[str, object]] = []
    error_rows: list[dict[str, object]] = []

    if hits.empty:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=base_output_columns()).to_csv(args.output, index=False)
        merged_output = args.merged_antibody_output or args.output.with_name("merged_antibodies.csv")
        pd.DataFrame(columns=merged_antibody_columns()).to_csv(merged_output, index=False)
        print("No BLAST hits found", file=sys.stderr)
        return

    threshold = normalize_identity_threshold(args.min_identity)
    hits = hits[hits["pident"] >= threshold].copy()
    hits = hits.sort_values(["evalue", "bitscore"], ascending=[True, False])
    hits = hits.drop_duplicates(["qseqid", "pdb_id", "antigen_id"])
    if args.max_hits:
        hits = hits.head(args.max_hits)

    seqres_keys: set[tuple[str, str]] = set()
    for _, hit in hits.iterrows():
        pdb_id = clean_cell(hit.get("pdb_id")).lower()
        antigen_chain = clean_cell(hit.get("antigen_id"))
        for sabdab_row in matching_sabdab_rows(summary, pdb_id, antigen_chain):
            for chain_id in split_chains(sabdab_row.get("Hchain")) + split_chains(sabdab_row.get("Lchain")):
                seqres_keys.add((pdb_id, chain_id))
    seqres_sequences = load_seqres_sequences(args.pdb_seqres, seqres_keys)

    aligner = make_pairwise_aligner()
    structure_cache: dict[str, list[StructureVariant]] = {}
    residue_cache: dict[tuple[str, str, str], list[ResidueRecord]] = {}
    epitope_cache: dict[tuple[str, str, str, tuple[str, ...], tuple[str, ...], str], EpitopeResult] = {}
    annotation_cache: dict[tuple[str, str], tuple[str, list[AnarciDomain], str]] = {}

    for _, hit in hits.iterrows():
        query_id = clean_cell(hit.get("qseqid"))
        pdb_id = clean_cell(hit.get("pdb_id")).lower()
        antigen_chain = clean_cell(hit.get("antigen_id"))

        if query_id not in queries:
            error_rows.append({"query_id": query_id, "pdb_id": pdb_id, "error": "query not found in FASTA"})
            continue
        if not pdb_id or not antigen_chain:
            error_rows.append({"query_id": query_id, "pdb_id": pdb_id, "error": "could not parse subject id"})
            continue

        structure_path = find_structure_file(args.pdb_cache, pdb_id)
        if structure_path is None:
            error_rows.append({"query_id": query_id, "pdb_id": pdb_id, "error": "structure file not found"})
            continue

        try:
            if pdb_id not in structure_cache:
                structure_cache[pdb_id] = load_structure_variants(structure_path)
            structure_variants = structure_cache[pdb_id]
        except Exception as exc:
            error_rows.append({"query_id": query_id, "pdb_id": pdb_id, "error": f"structure load failed: {exc}"})
            continue

        row_matches = matching_sabdab_rows(summary, pdb_id, antigen_chain)
        if not row_matches:
            error_rows.append(
                {
                    "query_id": query_id,
                    "pdb_id": pdb_id,
                    "antigen_id": antigen_chain,
                    "error": "no matching SAbDab protein antigen row",
                }
            )
            continue

        for sabdab_row in row_matches:
            heavy_ids, light_ids, contact_ids = antibody_chain_ids(sabdab_row, antigen_chain)
            contact_ids = unique_preserve_order(contact_ids)
            excluded_chain_ids = all_antibody_chains_for_pdb(summary, pdb_id)
            if not contact_ids:
                error_rows.append(
                    {
                        "query_id": query_id,
                        "pdb_id": pdb_id,
                        "antigen_id": antigen_chain,
                        "error": "no antibody chains for contact scoring",
                    }
                )
                continue

            selected_variant: StructureVariant | None = None
            selected_antigen_residues: list[ResidueRecord] = []
            selected_epitope: EpitopeResult | None = None
            variant_reasons: list[str] = []

            for variant in structure_variants:
                antigen_key = (pdb_id, variant.source, antigen_chain)
                if antigen_key not in residue_cache:
                    residue_cache[antigen_key] = chain_residues(variant.atoms, antigen_chain)
                antigen_residues = residue_cache[antigen_key]
                if not antigen_residues:
                    variant_reasons.append(f"{variant.source}: antigen chain {antigen_chain} not found")
                    continue

                epitope_key = (
                    pdb_id,
                    variant.source,
                    antigen_chain,
                    tuple(contact_ids),
                    tuple(excluded_chain_ids),
                    args.scheme,
                )
                if epitope_key not in epitope_cache:
                    epitope_cache[epitope_key] = identify_epitope_residues(
                        variant.atoms,
                        pdb_id,
                        antigen_chain,
                        contact_ids,
                        excluded_chain_ids,
                        args.scheme,
                        seqres_sequences,
                        args.ca_distance_threshold,
                        args.atom_distance_threshold,
                        args.homodimer_identity_threshold,
                        aligner,
                    )
                epitope = epitope_cache[epitope_key]
                if not epitope.contacts:
                    reason = "; ".join(epitope.cdr_warnings) or "no CDR-mediated contacts"
                    variant_reasons.append(f"{variant.source}: {reason}")
                    continue

                selected_variant = variant
                selected_antigen_residues = antigen_residues
                selected_epitope = epitope
                break

            if selected_variant is None or selected_epitope is None:
                error_rows.append(
                    {
                        "query_id": query_id,
                        "pdb_id": pdb_id,
                        "antigen_id": antigen_chain,
                        "error": "no matching CDR-mediated ab-ag complex in first assembly or asymmetric unit",
                        "details": "; ".join(variant_reasons),
                    }
                )
                continue

            antibody_annotation = annotate_antibodies(
                selected_variant.atoms,
                sabdab_row,
                args.scheme,
                annotation_cache,
                seqres_sequences,
            )
            epitope_scores = score_epitope(
                queries[query_id],
                selected_antigen_residues,
                selected_epitope,
                aligner,
            )
            structure_warnings = selected_variant.warning
            if selected_variant.source == "asymmetric_unit" and structure_warnings:
                print(f"{pdb_id}:{antigen_chain}: {structure_warnings}", file=sys.stderr)

            output_rows.append(
                {
                    "query_id": query_id,
                    "pdb_id": pdb_id,
                    "antigen_id": antigen_chain,
                    **antibody_annotation,
                    **epitope_scores,
                    "pident": hit.get("pident"),
                    "qcovs": hit.get("qcovs"),
                    "evalue": hit.get("evalue"),
                    "bitscore": hit.get("bitscore"),
                    "antigen_name": sabdab_row.get("antigen_name"),
                    "structure_source": selected_variant.source,
                    "structure_warnings": structure_warnings,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output = pd.DataFrame(output_rows)
    if output.empty:
        output = pd.DataFrame(columns=base_output_columns())
    else:
        columns = base_output_columns() + [column for column in output.columns if column not in base_output_columns()]
        output = output.reindex(columns=columns)
    output.to_csv(args.output, index=False)

    merged_output = args.merged_antibody_output or args.output.with_name("merged_antibodies.csv")
    make_merged_antibody_table(output).to_csv(merged_output, index=False)

    errors_output = args.errors_output or args.output.with_suffix(".errors.tsv")
    errors = pd.DataFrame(error_rows)
    if errors.empty:
        errors = pd.DataFrame(columns=ERROR_COLUMNS)
    else:
        columns = ERROR_COLUMNS + [column for column in errors.columns if column not in ERROR_COLUMNS]
        errors = errors.reindex(columns=columns)
    errors.to_csv(errors_output, sep="\t", index=False)
    print(f"wrote {len(output_rows)} scored rows to {args.output}", file=sys.stderr)
    if error_rows:
        print(f"wrote {len(error_rows)} error rows to {errors_output}", file=sys.stderr)


if __name__ == "__main__":
    main()
