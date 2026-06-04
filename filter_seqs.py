#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


MISSING_VALUES = {"", "nan", "none", "na", "n/a", "null", "-", "?"}
CHAIN_SPLIT_RE = re.compile(r"\s*\|\s*|[,;]\s*")


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


def normalize_pdb_id(value: object) -> str:
    return clean_cell(value).lower()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def iter_fasta(path: Path) -> Iterable[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []

    with open_text(path) as handle:
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


def parse_seqres_header(header: str) -> tuple[str, str, str, str]:
    record_id = header.split()[0]
    if "_" not in record_id:
        return "", "", record_id, ""

    pdb_id, chain_id = record_id.split("_", 1)
    mol_match = re.search(r"(?:^|\s)mol:([^\s]+)", header)
    mol_type = mol_match.group(1).lower() if mol_match else ""
    return pdb_id.lower(), chain_id, record_id, mol_type


def sabdab_targets(summary_path: Path, scope: str) -> tuple[set[str], set[tuple[str, str]]]:
    summary = pd.read_csv(summary_path, sep="\t", dtype=str)
    pdb_ids = {normalize_pdb_id(pdb) for pdb in summary["pdb"]}
    pdb_ids.discard("")

    chain_keys: set[tuple[str, str]] = set()

    if scope == "pdb":
        return pdb_ids, chain_keys

    for _, row in summary.iterrows():
        pdb_id = normalize_pdb_id(row.get("pdb"))
        if not pdb_id:
            continue

        if scope in {"antigen", "antibody-antigen"}:
            antigen_type = clean_cell(row.get("antigen_type")).lower()
            if "protein" in antigen_type:
                for chain in split_chains(row.get("antigen_chain")):
                    chain_keys.add((pdb_id, chain))

        if scope == "antibody-antigen":
            for column in ("Hchain", "Lchain"):
                for chain in split_chains(row.get(column)):
                    chain_keys.add((pdb_id, chain))

    return pdb_ids, chain_keys


def write_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def filter_fasta(
    fasta_path: Path,
    summary_path: Path,
    output_path: Path,
    scope: str,
) -> dict[str, int]:
    pdb_ids, chain_keys = sabdab_targets(summary_path, scope)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    stats = {
        "records_read": 0,
        "protein_records_read": 0,
        "records_written": 0,
        "target_pdb_ids": len(pdb_ids),
        "target_chain_keys": len(chain_keys),
    }

    with output_path.open("w") as out_handle:
        for header, sequence in iter_fasta(fasta_path):
            stats["records_read"] += 1
            pdb_id, chain_id, _, mol_type = parse_seqres_header(header)
            if mol_type != "protein":
                continue
            stats["protein_records_read"] += 1

            if scope == "pdb":
                keep = pdb_id in pdb_ids
            else:
                keep = (pdb_id, chain_id) in chain_keys

            if keep:
                write_record(out_handle, header, sequence)
                stats["records_written"] += 1

    return stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter PDB seqres FASTA to SAbDab protein chains."
    )
    parser.add_argument("--fasta", required=True, type=Path, help="PDB seqres FASTA, optionally .gz")
    parser.add_argument("--sabdab", required=True, type=Path, help="SAbDab summary TSV")
    parser.add_argument("--output", required=True, type=Path, help="Output FASTA")
    parser.add_argument(
        "--scope",
        choices=("antigen", "antibody-antigen", "pdb"),
        default="antigen",
        help=(
            "antigen: protein antigen chains only; antibody-antigen: protein antigen "
            "plus listed H/L chains; pdb: all protein chains from SAbDab PDB entries"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stats = filter_fasta(args.fasta, args.sabdab, args.output, args.scope)
    for key, value in stats.items():
        print(f"{key}: {value}", file=sys.stderr)


if __name__ == "__main__":
    main()
