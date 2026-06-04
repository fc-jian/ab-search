#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from filter_seqs import filter_fasta


BLAST_OUTFMT = (
    "6 qseqid sseqid pident length qlen slen qstart qend sstart send "
    "evalue bitscore qcovs"
)


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


def main() -> None:
    args = parse_args()
    reference_dir = args.reference_dir
    sabdab_summary = args.sabdab_summary or (reference_dir / "sabdab_summary_all.tsv")
    pdb_seqres = args.pdb_seqres or (reference_dir / "pdb_seqres.fasta.gz")
    sabdab_fasta = args.sabdab_fasta or (reference_dir / "sabdab_seq.fasta")
    blastdb_prefix = args.blastdb_prefix or (reference_dir / "blastdb")
    pdb_cache = reference_dir / "pdb_files"
    outdir = args.outdir or (Path("results") / query_stem(args.query))
    logs_dir = outdir / "logs"

    outdir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    if not sabdab_summary.exists():
        raise SystemExit(f"Missing SAbDab summary: {sabdab_summary}")
    if not pdb_seqres.exists():
        raise SystemExit(f"Missing PDB seqres FASTA: {pdb_seqres}")
    if not args.query.exists():
        raise SystemExit(f"Missing query FASTA: {args.query}")

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

    hits_path = outdir / "blastp_hits.tsv"
    if args.force_blast or not hits_path.exists() or hits_path.stat().st_size == 0:
        if shutil.which("blastp") is None:
            raise SystemExit("blastp not found on PATH")
        run_command(
            [
                "blastp",
                "-query",
                str(args.query),
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
            "download_pdb.py",
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
            "score_structures.py",
            "--query-fasta",
            str(args.query),
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

    print(f"Pipeline complete: {outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
