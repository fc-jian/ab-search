#!/usr/bin/env python3
"""Build a reusable, versioned index of published patent antibody evidence."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import sys
import tempfile
import urllib.request
from collections import Counter
from contextlib import closing
from pathlib import Path

from Bio.Data.IUPACData import protein_letters_3to1

from binding_sites import AMINO_ACIDS, HOSTS, PROTEIN_ALPHABET


POISK_URL = "https://zenodo.org/records/21775623/files/poisk_dataset.csv?download=1"
POISK_MD5 = "f0c1f098abf1756e4b9891aa4dc028ee"
PLABDAB_URL = "https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/paired_sequences.csv.gz"
CURATIONS = Path(__file__).resolve().parent / "data" / "verified_patent_sites.tsv"
SCHEMA_VERSION = 2
POISK_FIELDS = """unique_record_id patent_id antibody_identifier antibody_type vh.sequence_from_patent
vl.sequence_from_patent vh.seq_id vl.seq_id vh.is_full_heavy_chain vl.is_full_light_chain
target.organism target.target_type target.protein_name target.gene_names
target.uniprot_id target.sequence target.sequence_from_patent target.sequence_from_uniprot
target.sequence_accession target.seq_id target.domain_start target.domain_end target.domain_coordinates
target.binding_components target.binding_affinity.affinity_value target.binding_affinity.method
target.binding_conditions epitope.type epitope.linear_epitope_range.start epitope.linear_epitope_range.end
epitope.structural_contact_residues epitope.evidence_types audit.verdict audit.flags""".split()
PLABDAB_FIELDS = "ID heavy_sequence light_sequence pairing reference_title targets_mentioned".split()


def digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            hasher.update(block)
    return hasher.hexdigest()


def download(url: str, destination: Path, md5: str = "") -> None:
    if destination.is_file() and destination.stat().st_size:
        if not md5 or digest(destination, "md5") == md5:
            return
        raise ValueError(f"Cached file has a different checksum: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as tmp:
        staging = Path(tmp.name)
    try:
        with urllib.request.urlopen(url, timeout=120) as response, staging.open("wb") as output:
            shutil.copyfileobj(response, output)
        if not staging.stat().st_size or (md5 and digest(staging, "md5") != md5):
            raise ValueError(f"Incomplete download or checksum mismatch: {url}")
        if destination.suffix == ".gz":
            # Read through EOF to verify the compressed stream's length and CRC.
            with gzip.open(staging, "rb") as stream:
                while stream.read(1024 * 1024):
                    pass
        staging.replace(destination)
    finally:
        staging.unlink(missing_ok=True)


def csv_rows(path: Path, required: list[str]):
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, strict=True)
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(sorted(missing))}")
        for row in reader:
            if None in row or any(row[k] is None for k in required):
                raise ValueError(f"{path}: malformed CSV at line {reader.line_num}")
            yield {key: row[key].strip() for key in required}


def sequence(value: str) -> str:
    value = "".join(value.upper().split())
    if value in {"NA", "NAN", "NONE", "NULL", "UNKNOWN", "NOTSPECIFIED", "NOTAPPLICABLE",
                 "NOTFOUND", "UNAVAILABLE"} or not set(value) <= PROTEIN_ALPHABET:
        return ""
    return value


def pair_id(heavy: str, light: str) -> str:
    return hashlib.sha256(f"{heavy}\0{light}".encode("ascii")).hexdigest()


def checked_contacts(row: dict, reference: str) -> list[int]:
    if "Structure Determination" not in row["epitope.evidence_types"]:
        return []
    try:
        contacts = json.loads(row["epitope.structural_contact_residues"])
        if not isinstance(contacts, list):
            return []
        positions = set()
        for item in contacts:
            position = int(str(item["position"]))
            aa = item["residue"].strip().upper()
            if len(aa) == 3:
                aa = protein_letters_3to1.get(aa.title(), "")
            if aa not in AMINO_ACIDS or not (1 <= position <= len(reference)) or reference[position - 1] != aa:
                return []
            positions.add(position)
        # One matching letter alone offers little protection against a numbering offset.
        return sorted(positions) if len(positions) >= 2 else []
    except (AttributeError, KeyError, TypeError, ValueError):
        return []


def patent_record(row: dict, verified: dict | None) -> tuple[str, dict]:
    reference = sequence(row["target.sequence"])
    patent = sequence(row["target.sequence_from_patent"])
    uniprot = sequence(row["target.sequence_from_uniprot"])
    origin = "release_only"
    if reference and reference == patent:
        origin = "patent_listing"
    elif patent:
        origin = "listing_discrepancy"
    elif reference and reference == uniprot:
        origin = "uniprot"
    positions, site_type, citation, free_n = [], "", "", False
    if verified:
        if (row["patent_id"] != verified["patent_id"]
                or hashlib.sha256(reference.encode()).hexdigest() != verified["release_target_sha256"]):
            raise ValueError(f"Curated anchor no longer matches POISK {row['unique_record_id']}")
        peptide = verified["tested_peptide"]
        if peptide:
            if peptide not in reference or reference.find(peptide) != reference.rfind(peptide):
                raise ValueError(f"Tested peptide is ambiguous in {row['unique_record_id']}")
            reference = peptide
            origin = "verified_assay_peptide"
        region = verified["mapped_region"]
        if not region or region not in reference or reference.find(region) != reference.rfind(region):
            raise ValueError(f"Mapped region is ambiguous in {row['unique_record_id']}")
        start = reference.index(region) + 1
        positions = list(range(start, start + len(region)))
        site_type, citation = "verified_region", verified["source"]
        free_n = verified["free_n_terminus"] == "yes"
    elif origin == "patent_listing":
        positions = checked_contacts(row, reference)
        if positions:
            site_type = "structural_contacts"
    if not positions and (row["epitope.linear_epitope_range.start"] or row["epitope.linear_epitope_range.end"]
                          or row["epitope.structural_contact_residues"]):
        site_type = "unverified_coordinates"
    return reference, {
        "record_id": row["unique_record_id"], "document_id": row["patent_id"],
        "antibody": row["antibody_identifier"], "format": row["antibody_type"],
        "heavy": sequence(row["vh.sequence_from_patent"]), "light": sequence(row["vl.sequence_from_patent"]),
        "heavy_seq_id": row["vh.seq_id"], "light_seq_id": row["vl.seq_id"],
        "full_heavy_reported": row["vh.is_full_heavy_chain"], "full_light_reported": row["vl.is_full_light_chain"],
        "chain_basis": "patent_listing",
        "target": row["target.protein_name"], "organism": row["target.organism"],
        "gene_names": row["target.gene_names"], "uniprot": row["target.uniprot_id"],
        "target_type": row["target.target_type"], "origin": origin,
        "seq_id": row["target.seq_id"], "accession": row["target.sequence_accession"],
        "has_domain": any(row[f"target.{k}"] for k in ("domain_start", "domain_end", "domain_coordinates")),
        "components": row["target.binding_components"], "conditions": row["target.binding_conditions"],
        "affinity": row["target.binding_affinity.affinity_value"], "method": row["target.binding_affinity.method"],
        "audit": row["audit.verdict"], "flags": row["audit.flags"], "positions": positions,
        "site_type": site_type, "citation": citation, "free_n": free_n,
    }


def build_catalog(poisk: Path, plabdab: Path | None, output: Path, curations: Path | None = None) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    if curations is None and CURATIONS.is_file():
        curations = CURATIONS
    verified = {}
    if curations is not None:
        with curations.open(encoding="utf-8", newline="") as handle:
            verified = {row["record_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    fd, name = tempfile.mkstemp(dir=output.parent, suffix=".sqlite.tmp")
    os.close(fd)
    staging = Path(name)
    counts = Counter()
    targets = {}
    projectable, localized = set(), set()
    try:
        with closing(sqlite3.connect(staging)) as db, db:
            db.executescript("""
                CREATE TABLE targets(id INTEGER PRIMARY KEY, sequence TEXT NOT NULL UNIQUE,
                                     host INTEGER DEFAULT 0, localized INTEGER DEFAULT 0);
                CREATE TABLE records(id TEXT PRIMARY KEY, target_id INTEGER NOT NULL, data TEXT NOT NULL);
                CREATE TABLE publications(id TEXT PRIMARY KEY, pair_key TEXT, heavy TEXT, light TEXT,
                                          pairing TEXT, title TEXT, mentions TEXT);
                CREATE TABLE mentions(name TEXT, publication_id TEXT);
                CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
            """)
            for row in csv_rows(poisk, POISK_FIELDS):
                counts["poisk_rows"] += 1
                reference, data = patent_record(row, verified.get(row["unique_record_id"]))
                if not reference or not (data["heavy"] or data["light"]):
                    counts["no_target_or_chain"] += 1
                    continue
                if not data["record_id"] or not data["document_id"]:
                    raise ValueError("POISK row has no source identifier")
                if reference not in targets:
                    cursor = db.execute("INSERT INTO targets(sequence) VALUES(?)", (reference,))
                    targets[reference] = cursor.lastrowid
                target_id = targets[reference]
                if data["organism"].casefold() in HOSTS and data["target_type"] == "Protein":
                    projectable.add(target_id)
                    if data["positions"]:
                        localized.add(target_id)
                db.execute("INSERT INTO records VALUES(?,?,?)",
                           (data["record_id"], target_id, json.dumps(data, ensure_ascii=False)))
                counts["poisk_indexed"] += 1
            db.executemany("UPDATE targets SET host=1 WHERE id=?", ((i,) for i in projectable))
            db.executemany("UPDATE targets SET localized=1 WHERE id=?", ((i,) for i in localized))
            db.execute("CREATE INDEX records_by_target ON records(target_id)")
            if plabdab:
                for row in csv_rows(plabdab, PLABDAB_FIELDS):
                    counts["plabdab_rows"] += 1
                    heavy, light = sequence(row["heavy_sequence"]), sequence(row["light_sequence"])
                    if not heavy or not light:
                        continue
                    rid = row["ID"]
                    if not rid:
                        raise ValueError("PLAbDab row has no ID")
                    db.execute("INSERT INTO publications VALUES(?,?,?,?,?,?,?)",
                               (rid, pair_id(heavy, light), heavy, light, row["pairing"],
                                row["reference_title"], row["targets_mentioned"]))
                    names = {n.strip().casefold() for n in row["targets_mentioned"].split(";") if n.strip()}
                    db.executemany("INSERT INTO mentions VALUES(?,?)", ((n, rid) for n in names))
                    counts["plabdab_indexed"] += 1
            db.execute("CREATE INDEX publications_by_pair ON publications(pair_key)")
            db.execute("CREATE INDEX publications_by_mention ON mentions(name)")
            if not counts["poisk_indexed"] or (plabdab and not counts["plabdab_indexed"]):
                raise ValueError("No searchable records in a supplied source")
            counts["distinct_reference_sequences"] = len(targets)
            metadata = {"schema": SCHEMA_VERSION, "counts": dict(counts), "poisk_sha256": digest(poisk),
                        "plabdab_sha256": digest(plabdab) if plabdab else "",
                        "curations_sha256": digest(curations) if curations else "", "poisk_release": "1.1-compatible",
                        "poisk_url": POISK_URL, "plabdab_url": PLABDAB_URL if plabdab else ""}
            if digest(poisk, "md5") == POISK_MD5:
                metadata["poisk_release"] = "1.1"
            db.executemany("INSERT INTO metadata VALUES(?,?)",
                           ((k, json.dumps(v)) for k, v in metadata.items()))
        staging.replace(output)
        return metadata
    finally:
        staging.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference"))
    parser.add_argument("--download", action="store_true", help="Fetch the official POISK and PLAbDab tables")
    parser.add_argument("--poisk", type=Path, help="POISK 1.1 schema CSV or CSV.GZ")
    parser.add_argument("--plabdab", type=Path, help="Oxford paired-sequences CSV or CSV.GZ")
    parser.add_argument("--curations", type=Path, help="Optional TSV of independently verified patent regions")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    poisk = args.poisk or args.reference_dir / "poisk_v1.1.csv"
    plabdab = args.plabdab or args.reference_dir / "plabdab_paired.csv.gz"
    output = args.output or args.reference_dir / "patents.sqlite"
    if args.download:
        download(POISK_URL, poisk, POISK_MD5)
        download(PLABDAB_URL, plabdab)
    if not args.plabdab and not plabdab.exists():
        plabdab = None
    metadata = build_catalog(poisk, plabdab, output, curations=args.curations)
    print(json.dumps(metadata, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
