#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from filter_seqs import filter_fasta
from score_structures import make_merged_antibody_table, merged_antibody_columns


BLAST_OUTFMT = (
    "6 qseqid sseqid pident length qlen slen qstart qend sstart send "
    "evalue bitscore qcovs"
)
SCRIPT_DIR = Path(__file__).resolve().parent
SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass(frozen=True)
class QueryRecord:
    name: str
    description: str
    sequence: str


def run_command(command: list[str], log_path: Path | None = None) -> None:
    printable = " ".join(command)
    print(printable, file=sys.stderr)
    if log_path is None:
        subprocess.run(command, check=True)
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_handle:
        log_handle.write(f"$ {printable}\n")
        log_handle.flush()
        subprocess.run(
            command,
            check=True,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )


def blast_db_exists(prefix: Path) -> bool:
    return any(prefix.with_suffix(suffix).exists() for suffix in (".pin", ".psq", ".phr"))


def query_stem(path: Path) -> str:
    name = path.name
    for suffix in (".fasta", ".fa", ".faa", ".fas"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def read_query_records(path: Path) -> list[QueryRecord]:
    records: list[QueryRecord] = []
    description: str | None = None
    chunks: list[str] = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if description is not None:
                    records.append(make_query_record(description, chunks))
                description = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    if description is not None:
        records.append(make_query_record(description, chunks))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    names = [record.name for record in records]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate FASTA record IDs: {', '.join(duplicates)}")
    return records


def make_query_record(description: str, chunks: list[str]) -> QueryRecord:
    if not description:
        raise ValueError("FASTA record has an empty header")
    name = description.split()[0]
    sequence = "".join(chunks).upper().replace("*", "")
    if not sequence:
        raise ValueError(f"FASTA record {name!r} has an empty sequence")
    return QueryRecord(name=name, description=description, sequence=sequence)


def safe_record_dir_name(name: str) -> str:
    safe_name = SAFE_NAME_RE.sub("_", name).strip("._")
    if not safe_name:
        raise ValueError(f"FASTA record ID {name!r} cannot be converted to a directory name")
    return safe_name


def record_directory_names(records: list[QueryRecord]) -> dict[str, str]:
    directory_names = {record.name: safe_record_dir_name(record.name) for record in records}
    reverse: dict[str, list[str]] = {}
    for record_name, directory_name in directory_names.items():
        reverse.setdefault(directory_name, []).append(record_name)

    collisions = {
        directory_name: record_names
        for directory_name, record_names in reverse.items()
        if len(record_names) > 1
    }
    if collisions:
        details = "; ".join(
            f"{directory_name}: {', '.join(record_names)}"
            for directory_name, record_names in collisions.items()
        )
        raise ValueError(f"FASTA record IDs produce colliding directory names: {details}")
    return directory_names


def format_query_record(record: QueryRecord) -> str:
    lines = [f">{record.description}"]
    lines.extend(
        record.sequence[start : start + 80]
        for start in range(0, len(record.sequence), 80)
    )
    return "\n".join(lines) + "\n"


def write_query_record(record: QueryRecord, path: Path) -> bool:
    content = format_query_record(record)
    changed = not path.exists() or path.read_text() != content
    path.parent.mkdir(parents=True, exist_ok=True)
    if changed:
        path.write_text(content)
    return changed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the SAbDab antigen epitope search pipeline.")
    parser.add_argument("--query", default=Path("test_input.fasta"), type=Path)
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference"))
    parser.add_argument("--sabdab-summary", type=Path, default=None)
    parser.add_argument("--pdb-seqres", type=Path, default=None)
    parser.add_argument("--sabdab-fasta", type=Path, default=None)
    parser.add_argument("--blastdb-prefix", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-target-seqs", type=int, default=10000)
    parser.add_argument("--evalue", default="1e-5")
    parser.add_argument("--min-identity", type=float, default=0.7)
    parser.add_argument("--max-structures", type=int, default=None)
    parser.add_argument("--max-score-hits", type=int, default=200)
    parser.add_argument("--ca-distance-threshold", type=float, default=8.0)
    parser.add_argument("--atom-distance-threshold", type=float, default=4.5)
    parser.add_argument("--homodimer-identity-threshold", type=float, default=0.95)
    parser.add_argument("--scheme", default="chothia")
    parser.add_argument(
        "--scope",
        choices=("antigen", "antibody-antigen", "pdb"),
        default="antigen",
        help="Reference FASTA scope passed to filter_seqs.py",
    )
    parser.add_argument("--force-reference", action="store_true")
    parser.add_argument("--force-blast", action="store_true")
    return parser.parse_args()


def prepare_reference(
    args: argparse.Namespace,
    sabdab_summary: Path,
    pdb_seqres: Path,
    sabdab_fasta: Path,
    blastdb_prefix: Path,
    logs_dir: Path,
) -> None:
    if args.force_reference or not sabdab_fasta.exists() or sabdab_fasta.stat().st_size == 0:
        stats = filter_fasta(pdb_seqres, sabdab_summary, sabdab_fasta, args.scope)
        for key, value in stats.items():
            print(f"{key}: {value}", file=sys.stderr)

    if args.force_reference or not blast_db_exists(blastdb_prefix):
        if shutil.which("makeblastdb") is None:
            raise SystemExit("makeblastdb not found on PATH")
        run_command(
            [
                "makeblastdb",
                "-in",
                str(sabdab_fasta),
                "-dbtype",
                "prot",
                "-out",
                str(blastdb_prefix),
            ],
            logs_dir / "makeblastdb.log",
        )


def run_single_query(
    args: argparse.Namespace,
    query_path: Path,
    outdir: Path,
    sabdab_summary: Path,
    pdb_seqres: Path,
    blastdb_prefix: Path,
    pdb_cache: Path,
    force_blast: bool = False,
) -> Path:
    logs_dir = outdir / "logs"
    outdir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    hits_path = outdir / "blastp_hits.tsv"
    if force_blast or args.force_blast or not hits_path.exists() or hits_path.stat().st_size == 0:
        if shutil.which("blastp") is None:
            raise SystemExit("blastp not found on PATH")
        run_command(
            [
                "blastp",
                "-query",
                str(query_path),
                "-db",
                str(blastdb_prefix),
                "-out",
                str(hits_path),
                "-outfmt",
                BLAST_OUTFMT,
                "-max_target_seqs",
                str(args.max_target_seqs),
                "-evalue",
                str(args.evalue),
                "-num_threads",
                str(args.threads),
            ],
            logs_dir / "blastp.log",
        )

    run_command(
        [
            sys.executable,
            str(SCRIPT_DIR / "download_pdb.py"),
            "--hits",
            str(hits_path),
            "--pdb-cache",
            str(pdb_cache),
            "--min-identity",
            str(args.min_identity),
            "--manifest",
            str(outdir / "download_manifest.tsv"),
        ]
        + (
            ["--max-structures", str(args.max_structures)]
            if args.max_structures is not None
            else []
        ),
        logs_dir / "download_pdb.log",
    )

    run_command(
        [
            sys.executable,
            str(SCRIPT_DIR / "score_structures.py"),
            "--query-fasta",
            str(query_path),
            "--hits",
            str(hits_path),
            "--sabdab",
            str(sabdab_summary),
            "--pdb-seqres",
            str(pdb_seqres),
            "--pdb-cache",
            str(pdb_cache),
            "--output",
            str(outdir / "structure_scores.csv"),
            "--merged-antibody-output",
            str(outdir / "merged_antibodies.csv"),
            "--errors-output",
            str(outdir / "structure_score_errors.tsv"),
            "--min-identity",
            str(args.min_identity),
            "--max-hits",
            str(args.max_score_hits),
            "--ca-distance-threshold",
            str(args.ca_distance_threshold),
            "--atom-distance-threshold",
            str(args.atom_distance_threshold),
            "--homodimer-identity-threshold",
            str(args.homodimer_identity_threshold),
            "--scheme",
            str(args.scheme),
        ],
        logs_dir / "score_structures.log",
    )

    print(f"Query pipeline complete: {outdir}", file=sys.stderr)
    return outdir / "structure_scores.csv"


def write_multi_query_manifest(rows: list[dict[str, str]], output_path: Path) -> None:
    pd.DataFrame(
        rows,
        columns=[
            "query_id",
            "record_description",
            "output_subdir",
            "query_fasta",
            "structure_scores",
            "merged_antibodies",
        ],
    ).to_csv(output_path, sep="\t", index=False)


def write_all_query_merged_summary(score_paths: list[Path], output_path: Path) -> None:
    tables = []
    for score_path in score_paths:
        if score_path.exists():
            tables.append(pd.read_csv(score_path))

    all_scores = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    if all_scores.empty:
        summary = pd.DataFrame(columns=merged_antibody_columns())
    else:
        summary = make_merged_antibody_table(all_scores)
    summary.to_csv(output_path, index=False)


def main() -> None:
    args = parse_args()
    reference_dir = args.reference_dir
    sabdab_summary = args.sabdab_summary or (reference_dir / "sabdab_summary_all.tsv")
    pdb_seqres = args.pdb_seqres or (reference_dir / "pdb_seqres.fasta.gz")
    sabdab_fasta = args.sabdab_fasta or (reference_dir / "sabdab_seq.fasta")
    blastdb_prefix = args.blastdb_prefix or (reference_dir / "blastdb")
    pdb_cache = reference_dir / "pdb_files"
    outdir = args.outdir or (Path("results") / query_stem(args.query))

    if not sabdab_summary.exists():
        raise SystemExit(f"Missing SAbDab summary: {sabdab_summary}")
    if not pdb_seqres.exists():
        raise SystemExit(f"Missing PDB seqres FASTA: {pdb_seqres}")
    if not args.query.exists():
        raise SystemExit(f"Missing query FASTA: {args.query}")

    try:
        records = read_query_records(args.query)
        directory_names = record_directory_names(records)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    outdir.mkdir(parents=True, exist_ok=True)
    prepare_reference(
        args,
        sabdab_summary,
        pdb_seqres,
        sabdab_fasta,
        blastdb_prefix,
        outdir / "logs",
    )

    if len(records) == 1:
        run_single_query(
            args,
            args.query,
            outdir,
            sabdab_summary,
            pdb_seqres,
            blastdb_prefix,
            pdb_cache,
        )
        print(f"Pipeline complete: {outdir}", file=sys.stderr)
        return

    score_paths: list[Path] = []
    manifest_rows: list[dict[str, str]] = []
    for index, record in enumerate(records, start=1):
        record_outdir = outdir / directory_names[record.name]
        record_query_path = record_outdir / "query.fasta"
        query_changed = write_query_record(record, record_query_path)
        print(
            f"Running query {index}/{len(records)}: {record.name} -> {record_outdir}",
            file=sys.stderr,
        )
        score_path = run_single_query(
            args,
            record_query_path,
            record_outdir,
            sabdab_summary,
            pdb_seqres,
            blastdb_prefix,
            pdb_cache,
            force_blast=query_changed,
        )
        score_paths.append(score_path)
        manifest_rows.append(
            {
                "query_id": record.name,
                "record_description": record.description,
                "output_subdir": str(record_outdir),
                "query_fasta": str(record_query_path),
                "structure_scores": str(score_path),
                "merged_antibodies": str(record_outdir / "merged_antibodies.csv"),
            }
        )

    write_multi_query_manifest(manifest_rows, outdir / "query_manifest.tsv")
    write_all_query_merged_summary(score_paths, outdir / "merged_antibodies.csv")
    print(
        f"Multi-query pipeline complete: {len(records)} records; summary: "
        f"{outdir / 'merged_antibodies.csv'}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
