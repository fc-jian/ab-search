#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter
from contextlib import ExitStack, closing
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from filter_seqs import filter_fasta
from binding_search import Calibration, CANDIDATE_FIELDS, open_catalog, search_record, write_table
from binding_sites import PROTEIN_ALPHABET
from prepare_hil6 import RELEASE as AVIDA_RELEASE
from prepare_patents import digest


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
    return all(Path(f"{prefix}{suffix}").is_file() for suffix in (".pin", ".psq", ".phr")) or Path(f"{prefix}.pal").is_file()


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

    with path.open(encoding="utf-8-sig") as handle:
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
                if description is None:
                    raise ValueError("Sequence before first FASTA header")
                chunks.append(line)

    if description is not None:
        records.append(make_query_record(description, chunks))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    names = Counter(record.name for record in records)
    duplicates = sorted(name for name, count in names.items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate FASTA record IDs: {', '.join(duplicates)}")
    return records


def make_query_record(description: str, chunks: list[str]) -> QueryRecord:
    if not description:
        raise ValueError("FASTA record has an empty header")
    name = description.split()[0]
    sequence = "".join("".join(chunks).split()).upper().removesuffix("*")
    if not sequence:
        raise ValueError(f"FASTA record {name!r} has an empty sequence")
    invalid = set(sequence) - PROTEIN_ALPHABET
    if invalid:
        raise ValueError(f"Invalid protein characters in {name!r}: {''.join(sorted(invalid))!r}")
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
        reverse.setdefault(directory_name.casefold(), []).append(record_name)

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
    changed = not path.exists() or path.read_text(encoding="utf-8") != content
    path.parent.mkdir(parents=True, exist_ok=True)
    if changed:
        path.write_text(content, encoding="utf-8")
    return changed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Find published antibody sequences for a protein FASTA.")
    parser.add_argument("--query", default=Path("test_input.fasta"), type=Path)
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference"))
    parser.add_argument("--sources", nargs="+", choices=("sabdab", "patents", "avida"),
                        help="Default: use all prepared sources")
    parser.add_argument("--patent-db", type=Path, help="SQLite index built with prepare_patents.py")
    parser.add_argument("--avida-db", type=Path, help="Human IL-6 assay index built with prepare_hil6.py")
    parser.add_argument("--max-site-substitutions", type=int, default=0)
    parser.add_argument("--min-site-identity", type=float, default=0.90)
    parser.add_argument("--calibration-assays", type=Path, help="Independent target- and assay-matched binding labels")
    parser.add_argument("--free-antibody-nm", type=float, help="Free antibody concentration for 1:1 equilibrium occupancy")
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
    args = parser.parse_args()
    if args.max_site_substitutions < 0 or not 0 <= args.min_site_identity <= 1 or not 0 <= args.min_identity <= 100:
        parser.error("Invalid site substitution limit or identity threshold")
    if args.free_antibody_nm is not None and (not math.isfinite(args.free_antibody_nm) or args.free_antibody_nm <= 0):
        parser.error("--free-antibody-nm must be finite and positive")
    return args


def prepare_reference(
    args: argparse.Namespace,
    sabdab_summary: Path,
    pdb_seqres: Path,
    sabdab_fasta: Path,
    blastdb_prefix: Path,
    logs_dir: Path,
) -> None:
    source_stamp = Path(f"{sabdab_fasta}.sources.json")
    sources = json.dumps({"scope": args.scope, "sabdab_sha256": digest(sabdab_summary),
                          "seqres_sha256": digest(pdb_seqres)}, sort_keys=True)
    changed = source_stamp.exists() and source_stamp.read_text() != sources
    rebuild = args.force_reference or changed or not sabdab_fasta.exists() or sabdab_fasta.stat().st_size == 0
    if not args.sabdab_fasta and not source_stamp.exists():
        rebuild = True
    if rebuild:
        sabdab_fasta.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=sabdab_fasta.parent) as staging:
            new_fasta = Path(staging) / "antigens.fasta"
            stats = filter_fasta(pdb_seqres, sabdab_summary, new_fasta, args.scope)
            new_fasta.replace(sabdab_fasta)
        for key, value in stats.items():
            print(f"{key}: {value}", file=sys.stderr)
    index_stamp = Path(f"{blastdb_prefix}.fasta.sha256")
    fasta_digest = digest(sabdab_fasta)
    index_changed = index_stamp.exists() and index_stamp.read_text() != fasta_digest
    if rebuild or index_changed or not blast_db_exists(blastdb_prefix):
        if shutil.which("makeblastdb") is None:
            raise SystemExit("makeblastdb not found on PATH")
        blastdb_prefix.parent.mkdir(parents=True, exist_ok=True)
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
        index_stamp.write_text(fasta_digest)
    if rebuild:
        source_stamp.write_text(sources)


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
    stamp = logs_dir / "blast_inputs.json"
    db_files = sorted(blastdb_prefix.parent.glob(blastdb_prefix.name + ".*"))
    fingerprint = json.dumps({
        "query": hashlib.sha256(query_path.read_bytes()).hexdigest(),
        "db": str(blastdb_prefix.resolve()), "files": [(p.name, p.stat().st_size, p.stat().st_mtime_ns) for p in db_files],
        "format": BLAST_OUTFMT, "evalue": args.evalue, "max_targets": args.max_target_seqs,
    }, sort_keys=True)
    if (force_blast or args.force_blast or not hits_path.exists()
            or not stamp.exists() or stamp.read_text() != fingerprint):
        if shutil.which("blastp") is None:
            raise SystemExit("blastp not found on PATH")
        with tempfile.TemporaryDirectory(dir=outdir) as staging:
            new_hits = Path(staging) / "hits.tsv"
            run_command([
                "blastp",
                "-query",
                str(query_path),
                "-db",
                str(blastdb_prefix),
                "-out",
                str(new_hits),
                "-outfmt",
                BLAST_OUTFMT,
                "-max_target_seqs",
                str(args.max_target_seqs),
                "-evalue",
                str(args.evalue),
                "-num_threads",
                str(args.threads),
            ], logs_dir / "blastp.log")
            new_hits.replace(hits_path)
            stamp.write_text(fingerprint)

    with tempfile.TemporaryDirectory(dir=outdir) as staging:
        manifest = Path(staging) / "download_manifest.tsv"
        download_command = [
            sys.executable,
            str(SCRIPT_DIR / "download_pdb.py"),
            "--hits",
            str(hits_path),
            "--pdb-cache",
            str(pdb_cache),
            "--min-identity",
            str(args.min_identity),
            "--manifest",
            str(manifest),
        ]
        if args.max_structures is not None:
            download_command.extend(["--max-structures", str(args.max_structures)])
        try:
            run_command(download_command, logs_dir / "download_pdb.log")
        except subprocess.CalledProcessError:
            if not manifest.is_file():
                raise
            with manifest.open(newline="", encoding="utf-8") as handle:
                failures = sum(r["status"] == "failed" for r in csv.DictReader(handle, delimiter="\t"))
            if not failures:
                raise
            print(f"{failures} PDB downloads failed. Scoring available structures. See {manifest.name}.", file=sys.stderr)
        manifest.replace(outdir / "download_manifest.tsv")

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
            "binding_candidates",
            "binding_evidence",
            "plabdab_references",
        ],
    ).to_csv(output_path, sep="\t", index=False)


def write_all_query_merged_summary(score_paths: list[Path], output_path: Path) -> None:
    from score_structures import make_merged_antibody_table, merged_antibody_columns
    tables = []
    for score_path in score_paths:
        if score_path.exists():
            tables.append(pd.read_csv(score_path, dtype=str, keep_default_na=False))

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
    patent_db = args.patent_db or reference_dir / "patents.sqlite"
    panel_db = args.avida_db or reference_dir / "avida_hil6.sqlite"
    if args.sources is None:
        use_sabdab = sabdab_summary.is_file() and pdb_seqres.is_file()
        use_patents = args.patent_db is not None or patent_db.is_file()
        use_panel = args.avida_db is not None or panel_db.is_file()
    else:
        use_sabdab, use_patents, use_panel = (s in args.sources for s in ("sabdab", "patents", "avida"))
    if not (use_patents or use_sabdab or use_panel):
        raise SystemExit("No prepared sources. Run python prepare_patents.py --download or prepare SAbDab references.")
    if use_patents and not patent_db.is_file():
        raise SystemExit(f"Missing patent index: {patent_db}. Build it with prepare_patents.py.")
    if use_panel and not panel_db.is_file():
        raise SystemExit(f"Missing human IL-6 assay index: {panel_db}. Build it with prepare_hil6.py.")
    if use_sabdab and not sabdab_summary.exists():
        raise SystemExit(f"Missing SAbDab summary: {sabdab_summary}")
    if use_sabdab and not pdb_seqres.exists():
        raise SystemExit(f"Missing PDB seqres FASTA: {pdb_seqres}")
    if not args.query.exists():
        raise SystemExit(f"Missing query FASTA: {args.query}")

    try:
        records = read_query_records(args.query)
        directory_names = record_directory_names(records)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    outdir.mkdir(parents=True, exist_ok=True)
    sources = [s for s, enabled in (("SAbDab", use_sabdab), ("patents", use_patents), ("AVIDa-hIL6", use_panel)) if enabled]
    print(f"Searching {' + '.join(sources)}", file=sys.stderr)
    if use_sabdab:
        prepare_reference(args, sabdab_summary, pdb_seqres, sabdab_fasta, blastdb_prefix, outdir / "logs")
    calibration = Calibration(args.calibration_assays)
    metadata = {
        "sources": sources, "query_sha256": digest(args.query),
        "max_site_substitutions": args.max_site_substitutions, "min_site_identity": args.min_site_identity,
        "min_alignment_identity": args.min_identity / 100 if args.min_identity > 1 else args.min_identity,
        "free_antibody_nM": args.free_antibody_nm,
        "calibration_assays_sha256": digest(args.calibration_assays) if args.calibration_assays else None,
    }
    if use_sabdab:
        metadata["sabdab"] = {"summary_sha256": digest(sabdab_summary), "seqres_sha256": digest(pdb_seqres),
                              "antigen_fasta_sha256": digest(sabdab_fasta), "scope": args.scope, "scheme": args.scheme,
                              "ca_distance": args.ca_distance_threshold, "atom_distance": args.atom_distance_threshold,
                              "homodimer_identity": args.homodimer_identity_threshold, "blast_evalue": args.evalue,
                              "max_targets": args.max_target_seqs, "max_structures": args.max_structures,
                              "max_score_hits": args.max_score_hits}
    multi = len(records) > 1
    score_paths: list[Path] = []
    manifest_rows: list[dict[str, str]] = []
    all_candidates = []
    with ExitStack() as resources:
        db = resources.enter_context(closing(open_catalog(patent_db))) if use_patents else None
        if db is not None:
            metadata["patents"] = {k: json.loads(v) for k, v in db.execute("SELECT key,value FROM metadata")}
        panel = None
        if use_panel:
            panel = resources.enter_context(closing(sqlite3.connect(panel_db.resolve().as_uri() + "?mode=ro", uri=True)))
            if panel.execute("SELECT revision FROM release").fetchone() != (AVIDA_RELEASE,):
                raise ValueError("Unexpected AVIDa panel release. Rebuild with prepare_hil6.py.")
            metadata["avida_revision"] = AVIDA_RELEASE
        for record in records:
            record_outdir = outdir / directory_names[record.name] if multi else outdir
            record_outdir.mkdir(parents=True, exist_ok=True)
            record_query_path = record_outdir / "query.fasta"
            if record_query_path.resolve() == args.query.resolve():
                record_query_path = record_outdir / "normalized_query.fasta"
            changed = write_query_record(record, record_query_path)
            score_path = None
            if use_sabdab:
                score_path = run_single_query(args, record_query_path, record_outdir, sabdab_summary,
                                              pdb_seqres, blastdb_prefix, pdb_cache, force_blast=changed)
                score_paths.append(score_path)
            candidates = search_record(record, db, panel, score_path, record_outdir, args, calibration)
            print(f"{record.name}: {len(candidates)} candidate antibody sequences", file=sys.stderr)
            if multi:
                all_candidates.extend(candidates)
            structure_errors = 0
            if score_path:
                with (record_outdir / "structure_score_errors.tsv").open(encoding="utf-8", newline="") as errors:
                    structure_errors = sum(1 for _ in csv.DictReader(errors, delimiter="\t"))
                if structure_errors:
                    print(f"{record.name}: {structure_errors} structure records require review", file=sys.stderr)
            manifest_rows.append({
                "query_id": record.name,
                "record_description": record.description,
                "output_subdir": str(record_outdir),
                "query_fasta": str(record_query_path),
                "structure_scores": str(score_path) if score_path else "",
                "structure_errors": structure_errors,
                "merged_antibodies": str(record_outdir / "merged_antibodies.csv") if score_path else "",
                "binding_candidates": str(record_outdir / "binding_candidates.csv"),
                "binding_evidence": str(record_outdir / "binding_evidence.csv"),
                "plabdab_references": str(record_outdir / "plabdab_references.csv"),
            })
    if multi:
        write_multi_query_manifest(manifest_rows, outdir / "query_manifest.tsv")
        write_table(outdir / "binding_candidates.csv", CANDIDATE_FIELDS, all_candidates)
        if use_sabdab:
            write_all_query_merged_summary(score_paths, outdir / "merged_antibodies.csv")
    metadata["records"] = manifest_rows
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Pipeline complete: {outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
