#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from pathlib import Path


BLAST_COLUMNS = [
    "query_id",
    "subject_id",
    "pident",
    "length",
    "query_length",
    "subject_length",
    "query_coverage",
    "evalue",
    "bitscore",
    "subject_description",
]


def normalize_identity_threshold(value: float) -> float:
    return value * 100.0 if value <= 1.0 else value


def read_sequences(path: Path) -> list[str]:
    sequences: list[str] = []
    with path.open() as handle:
        for line in handle:
            sequence = line.strip().upper().replace("*", "")
            if sequence:
                sequences.append(sequence)
    return sequences


def write_query_fasta(sequences: list[str], path: Path) -> list[str]:
    query_ids = []
    with path.open("w") as handle:
        for index, sequence in enumerate(sequences, start=1):
            query_id = f"seq_{index}"
            query_ids.append(query_id)
            handle.write(f">{query_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    return query_ids


def run_blastp(
    blastp_exe: str,
    query_fasta: Path,
    blastdb: str,
    output_path: Path,
    evalue: str,
    max_target_seqs: int,
    threads: int,
) -> None:
    outfmt = (
        "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle"
    )
    subprocess.run(
        [
            blastp_exe,
            "-query",
            str(query_fasta),
            "-db",
            blastdb,
            "-out",
            str(output_path),
            "-outfmt",
            outfmt,
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            str(max_target_seqs),
            "-num_threads",
            str(threads),
        ],
        check=True,
    )


def parse_blast_hits(path: Path) -> dict[str, list[dict[str, object]]]:
    hits_by_query: dict[str, list[dict[str, object]]] = {}
    if not path.exists() or path.stat().st_size == 0:
        return hits_by_query

    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if len(fields) < len(BLAST_COLUMNS):
                continue
            row = dict(zip(BLAST_COLUMNS, fields[: len(BLAST_COLUMNS)]))
            row["pident"] = float(row["pident"])
            row["length"] = int(float(row["length"]))
            row["query_length"] = int(float(row["query_length"]))
            row["subject_length"] = int(float(row["subject_length"]))
            row["query_coverage"] = float(row["query_coverage"])
            row["evalue"] = float(row["evalue"])
            row["bitscore"] = float(row["bitscore"])
            row["subject_description"] = strip_subject_id_from_title(
                str(row["subject_id"]),
                str(row["subject_description"]),
            )
            hits_by_query.setdefault(str(row["query_id"]), []).append(row)
    return hits_by_query


def strip_subject_id_from_title(subject_id: str, title: str) -> str:
    title = title.strip()
    if title == subject_id:
        return ""
    prefix = f"{subject_id} "
    if title.startswith(prefix):
        return title[len(prefix) :].strip()
    return title


def best_hit(
    hits: list[dict[str, object]],
    min_identity: float,
    min_query_coverage: float,
) -> dict[str, object] | None:
    threshold = normalize_identity_threshold(min_identity)
    coverage_threshold = min_query_coverage * 100.0 if min_query_coverage <= 1.0 else min_query_coverage
    filtered = [
        hit
        for hit in hits
        if float(hit["pident"]) >= threshold
        and float(hit["query_coverage"]) >= coverage_threshold
    ]
    if not filtered:
        return None

    return sorted(
        filtered,
        key=lambda hit: (
            float(hit["evalue"]),
            -float(hit["bitscore"]),
            -float(hit["pident"]),
            -float(hit["query_coverage"]),
        ),
    )[0]


def write_output(
    output_path: Path,
    sequences: list[str],
    query_ids: list[str],
    hits_by_query: dict[str, list[dict[str, object]]],
    min_identity: float,
    min_query_coverage: float,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "query_index",
        "query_id",
        "query_sequence",
        "matched",
        "subject_id",
        "subject_description",
        "pident",
        "query_coverage",
        "evalue",
        "bitscore",
        "alignment_length",
        "query_length",
        "subject_length",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for index, (query_id, sequence) in enumerate(zip(query_ids, sequences), start=1):
            hit = best_hit(
                hits_by_query.get(query_id, []),
                min_identity=min_identity,
                min_query_coverage=min_query_coverage,
            )
            if hit is None:
                writer.writerow(
                    {
                        "query_index": index,
                        "query_id": query_id,
                        "query_sequence": sequence,
                        "matched": False,
                    }
                )
                continue

            writer.writerow(
                {
                    "query_index": index,
                    "query_id": query_id,
                    "query_sequence": sequence,
                    "matched": True,
                    "subject_id": hit["subject_id"],
                    "subject_description": hit["subject_description"],
                    "pident": hit["pident"],
                    "query_coverage": hit["query_coverage"],
                    "evalue": hit["evalue"],
                    "bitscore": hit["bitscore"],
                    "alignment_length": hit["length"],
                    "query_length": hit["query_length"],
                    "subject_length": hit["subject_length"],
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find best matching BLAST database sequence descriptions for plain "
            "one-sequence-per-line protein input."
        )
    )
    parser.add_argument("--input", required=True, type=Path, help="TXT with one protein sequence per line")
    parser.add_argument("--blastdb", required=True, help="BLAST database prefix")
    parser.add_argument("--output", required=True, type=Path, help="Output CSV")
    parser.add_argument("--blastp", default=None, help="blastp executable path; defaults to PATH lookup")
    parser.add_argument("--min-identity", type=float, default=0.7)
    parser.add_argument("--min-query-coverage", type=float, default=0.0)
    parser.add_argument("--evalue", default="1e-5")
    parser.add_argument("--max-target-seqs", type=int, default=25)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--keep-temp",
        type=Path,
        default=None,
        help="Optional directory to keep generated query FASTA and BLAST TSV",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sequences = read_sequences(args.input)
    if not sequences:
        raise SystemExit(f"No sequences found in {args.input}")
    blastp_exe = args.blastp or shutil.which("blastp")
    if blastp_exe is None:
        raise SystemExit("blastp not found on PATH; pass --blastp /path/to/blastp")

    if args.keep_temp:
        args.keep_temp.mkdir(parents=True, exist_ok=True)
        query_fasta = args.keep_temp / "find_chain_name_queries.fasta"
        blast_output = args.keep_temp / "find_chain_name_blast.tsv"
        query_ids = write_query_fasta(sequences, query_fasta)
        run_blastp(
            blastp_exe,
            query_fasta,
            args.blastdb,
            blast_output,
            args.evalue,
            args.max_target_seqs,
            args.threads,
        )
        hits_by_query = parse_blast_hits(blast_output)
    else:
        with tempfile.TemporaryDirectory(prefix="find_chain_name_") as tmpdir:
            tmpdir_path = Path(tmpdir)
            query_fasta = tmpdir_path / "queries.fasta"
            blast_output = tmpdir_path / "blast.tsv"
            query_ids = write_query_fasta(sequences, query_fasta)
            run_blastp(
                blastp_exe,
                query_fasta,
                args.blastdb,
                blast_output,
                args.evalue,
                args.max_target_seqs,
                args.threads,
            )
            hits_by_query = parse_blast_hits(blast_output)

    write_output(
        args.output,
        sequences,
        query_ids,
        hits_by_query,
        min_identity=args.min_identity,
        min_query_coverage=args.min_query_coverage,
    )


if __name__ == "__main__":
    main()
