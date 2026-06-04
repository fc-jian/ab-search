#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd


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


def normalize_identity_threshold(value: float) -> float:
    return value * 100.0 if value <= 1.0 else value


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
    if "pident" in hits:
        hits["pident"] = pd.to_numeric(hits["pident"], errors="coerce")
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


def unique_pdb_ids_from_hits(
    hits: pd.DataFrame,
    min_identity: float,
    max_structures: int | None,
) -> list[str]:
    if hits.empty:
        return []

    threshold = normalize_identity_threshold(min_identity)
    filtered = hits[hits["pident"] >= threshold].copy()
    if filtered.empty:
        return []

    filtered[["pdb_id", "chain_id"]] = filtered["sseqid"].apply(
        lambda value: pd.Series(parse_subject_id(value))
    )
    filtered["bitscore"] = pd.to_numeric(filtered.get("bitscore"), errors="coerce")
    filtered["evalue"] = pd.to_numeric(filtered.get("evalue"), errors="coerce")
    filtered = filtered.sort_values(["evalue", "bitscore"], ascending=[True, False])

    pdb_ids: list[str] = []
    seen: set[str] = set()
    for pdb_id in filtered["pdb_id"]:
        if not pdb_id or pdb_id in seen:
            continue
        pdb_ids.append(pdb_id)
        seen.add(pdb_id)
        if max_structures is not None and len(pdb_ids) >= max_structures:
            break
    return pdb_ids


def download_pdb_cif(
    pdb_id: str,
    cache_dir: Path,
    url_template: str,
    timeout: int,
) -> tuple[Path, str]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.lower()
    output_path = cache_dir / f"{pdb_id}.cif.gz"
    if output_path.exists() and output_path.stat().st_size > 0:
        return output_path, "cached"

    url = url_template.format(pdb_id=pdb_id, PDB_ID=pdb_id.upper())
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    request = urllib.request.Request(url, headers={"User-Agent": "ab-search/0.1"})

    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            tmp_path.write_bytes(response.read())
        tmp_path.replace(output_path)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        if tmp_path.exists():
            tmp_path.unlink()
        raise RuntimeError(f"failed to download {pdb_id} from {url}: {exc}") from exc

    return output_path, "downloaded"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download PDB mmCIF files for BLAST hits.")
    parser.add_argument("--hits", required=True, type=Path, help="BLAST tabular output")
    parser.add_argument("--pdb_cache", "--pdb-cache", dest="pdb_cache", required=True, type=Path)
    parser.add_argument("--min-identity", type=float, default=0.7)
    parser.add_argument("--max-structures", type=int, default=None)
    parser.add_argument(
        "--url-template",
        default="https://files.rcsb.org/download/{PDB_ID}.cif.gz",
        help="Download URL template with {pdb_id} and {PDB_ID} placeholders",
    )
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--manifest", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hits = read_blast_hits(args.hits)
    pdb_ids = unique_pdb_ids_from_hits(hits, args.min_identity, args.max_structures)

    rows = []
    failures = 0
    for pdb_id in pdb_ids:
        try:
            path, status = download_pdb_cif(
                pdb_id,
                args.pdb_cache,
                args.url_template,
                args.timeout,
            )
            rows.append({"pdb_id": pdb_id, "path": str(path), "status": status})
            print(f"{pdb_id}: {status}", file=sys.stderr)
        except RuntimeError as exc:
            failures += 1
            rows.append({"pdb_id": pdb_id, "path": "", "status": "failed", "error": str(exc)})
            print(str(exc), file=sys.stderr)

    manifest = args.manifest or (args.pdb_cache / "download_manifest.tsv")
    pd.DataFrame(rows).to_csv(manifest, sep="\t", index=False)

    if failures:
        raise SystemExit(f"{failures} PDB downloads failed")


if __name__ == "__main__":
    main()
