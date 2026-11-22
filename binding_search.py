"""Combine structural and patent evidence for one fixed query at a time."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sqlite3
import tempfile
from collections import defaultdict
from pathlib import Path

from binding_sites import AMINO_ACIDS, HOSTS, assess_site, sequence_match
from prepare_hil6 import assay_evidence
from prepare_patents import SCHEMA_VERSION, pair_id


GRADES = {"exact_structure": 4, "exact_assay": 4, "exact_report": 3, "site_conserved": 2, "site_similar": 1}
EVIDENCE_FIELDS = [
    "query_id", "source", "record_id", "document_id", "antibody", "heavy", "light", "target",
    "organism", "uniprot", "accession", "seq_id", "origin", "match", "alignment_identity",
    "site_type", "site_size", "site_covered", "site_coverage", "site_substitutions", "site_identity",
    "tier", "reason", "audit", "flags", "method", "affinity", "conditions", "citation",
    "free_antibody_nM", "equilibrium_bound_fraction",
    "binding_probability", "calibration_n", "calibration_assay", "probability_status",
    "chain_basis", "heavy_seq_id", "light_seq_id", "full_heavy_reported", "full_light_reported",
]
CANDIDATE_FIELDS = [
    "query_id", "query_sha256", "antibody_key", "heavy", "light", "antibody_names",
    "evidence_tier", "evidence_grade", "sources", "source_documents", "evidence_records",
    "best_evidence_record", "chain_basis", "site_type", "site_coverage", "site_substitutions", "site_identity",
    "binding_probability", "calibration_n", "calibration_assay", "probability_status",
]
PUBLICATION_FIELDS = ["query_id", "relationship", "ID", "heavy_sequence", "light_sequence",
                      "pairing", "reference_title", "targets_mentioned"]
NEGATIVE = re.compile(r"\b(no (?:detectable |measurable )?(?:binding|interaction)|non[- ]?bind\w*|"
                      r"not detect\w*|(?:does|did|could) not bind|not bind\w*)\b", re.I)
KD = re.compile(r"K\s*[_]?D\s*[:=]\s*(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*(pM|nM|[uµμ]M|mM|M)",
                re.IGNORECASE | re.ASCII)


def open_catalog(path: Path) -> sqlite3.Connection:
    db = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
    try:
        schema = db.execute("SELECT value FROM metadata WHERE key='schema'").fetchone()
        if not schema or json.loads(schema[0]) != SCHEMA_VERSION:
            raise ValueError("Unsupported patent catalog schema. Rebuild with prepare_patents.py.")
    except Exception:
        db.close()
        raise
    return db


def classify_patent(data, query, reference, match, max_substitutions, min_site_identity):
    common = {k: data[k] for k in EVIDENCE_FIELDS if k in data}
    common.update(source="POISK", match=match.kind, alignment_identity=match.identity, tier="review")
    if not common["citation"]:
        common["citation"] = f"https://patents.google.com/patent/US{data['document_id']}"
    if data["audit"] != "CONFIRMED" or data["flags"]:
        return {**common, "reason": "patent_audit_requires_review"}
    paired = (data["format"] == "Conventional" and data["heavy"] and data["light"])
    vhh = data["format"] == "VHH" and data["heavy"] and not data["light"]
    if (not (paired or vhh) or not set(data["heavy"] + data["light"]) <= AMINO_ACIDS
            or any(len(chain) < 65 for chain in (data["heavy"], data["light"]) if chain)):
        return {**common, "reason": "incomplete_or_ambiguous_antibody"}
    if data["target_type"] != "Protein" or data["components"]:
        return {**common, "reason": "target_context_required"}
    curated = data["site_type"] == "verified_region"
    if not curated and (data["origin"] != "patent_listing" or data["has_domain"]):
        return {**common, "reason": "assayed_construct_unresolved"}
    if data["site_type"] == "unverified_coordinates":
        return {**common, "reason": "coordinate_basis_unverified"}
    if match.ambiguous:
        return {**common, "reason": "ambiguous_alignment"}
    if data["free_n"] and match.mapping.get(0) != 0:
        return {**common, "reason": "requires_free_n_terminus"}
    if NEGATIVE.search(data["affinity"]):
        return {**common, "reason": "reported_negative_or_conditional"}
    if match.kind == "exact" and not set(reference) <= AMINO_ACIDS:
        return {**common, "reason": "ambiguous_reference_sequence"}
    if data["positions"]:
        checked = assess_site(query, reference, match, data["positions"],
                              max_substitutions, min_site_identity, linear=curated)
        common.update(checked)
        if checked.get("reason"):
            return common
    elif match.kind != "exact":
        return {**common, "reason": "site_not_localized"}
    if match.kind == "exact":
        common["tier"] = "exact_report"
    return common


def bound_fraction(affinity: str, method: str, free_nm: float) -> float | None:
    match = KD.fullmatch(affinity)
    if not match or not re.search(r"\b(SPR|BLI|Biacore|KinExA|Octet)\b", method, re.I):
        return None
    units = {"pm": .001, "nm": 1, "um": 1000, "µm": 1000, "μm": 1000, "mm": 1e6, "m": 1e9}
    kd_nm = float(match[1]) * units[match[2].lower()]
    return free_nm / (free_nm + kd_nm) if kd_nm > 0 and math.isfinite(kd_nm) else None


def patent_evidence(db, query, min_identity, max_substitutions, min_site_identity, free_nm):
    rows, aliases = [], set()
    for target_id, reference, host, localized in db.execute("SELECT * FROM targets ORDER BY id"):
        if not host and query != reference:
            continue
        match = sequence_match(query, reference, allow_local=bool(localized), min_identity=min_identity)
        if match is None:
            continue
        for (serialized,) in db.execute("SELECT data FROM records WHERE target_id=? ORDER BY id", (target_id,)):
            data = json.loads(serialized)
            if match.kind != "exact" and (data["organism"].casefold() not in HOSTS or data["target_type"] != "Protein"):
                continue
            row = classify_patent(data, query, reference, match, max_substitutions, min_site_identity)
            if match.kind == "exact" and data["audit"] == "CONFIRMED" and not data["flags"]:
                aliases.update(s.strip().casefold() for s in re.split(r"[,;|]", data["gene_names"]) if s.strip())
                if data["target"]:
                    aliases.add(data["target"].casefold())
            if free_nm is not None and row["tier"] == "exact_report" and data["origin"] == "patent_listing":
                fraction = bound_fraction(data["affinity"], data["method"], free_nm)
                if fraction is not None:
                    row.update(free_antibody_nM=free_nm, equilibrium_bound_fraction=fraction)
            rows.append(row)
    return rows, aliases


def structure_evidence(path: Path, query: str, min_identity, max_substitutions, min_site_identity):
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for r in csv.DictReader(handle):
            ref = r.get("subject_seqres_sequence", "")
            host = r.get("antigen_species", "").casefold() in HOSTS
            positions = [int(p) for p in r.get("host_contact_positions", "").split("|") if p]
            # Keep the chain pair together when choosing variable domains or full chains.
            heavy, light = r["seq_heavy_full"], r["seq_light_full"]
            chain_basis = "deposited_chains"
            if r["seq_heavy_anarci"] and (r["seq_light_anarci"] or not light):
                heavy, light = r["seq_heavy_anarci"], r["seq_light_anarci"]
                chain_basis = "numbered_domains"
            rid = ":".join(r[k] for k in ("pdb_id", "antigen_id", "heavy_chain_id", "light_chain_id"))
            row = dict(source="SAbDab", record_id=rid, document_id=r["pdb_id"], heavy=heavy, light=light,
                       target=r["antigen_name"], organism=r.get("antigen_species", ""),
                       origin="PDB_SEQRES", site_type="CDR_contacts", tier="review",
                       chain_basis=chain_basis, heavy_seq_id=r["heavy_chain_id"], light_seq_id=r["light_chain_id"],
                       citation=f"https://www.rcsb.org/structure/{r['pdb_id']}")
            rows.append(row)
            vhh = re.search(r"\b(VHH|nanobody|single[- ]domain)\b", r.get("heavy_chain_description", ""), re.I)
            if not heavy or not (light or vhh) or not set(heavy + light) <= AMINO_ACIDS:
                row["reason"] = "incomplete_or_ambiguous_antibody"
                continue
            if r.get("cdr_contact_warnings") or r.get("anarci_warnings"):
                row["reason"] = "incomplete_CDR_annotation"
                continue
            if r.get("homodimer_excluded_chain_ids") or int(r["number_multichain_epitope_residues"]) > 0:
                row["reason"] = "multichain_site"
                continue
            if not ref or (query != ref and not host):
                row["reason"] = "reference_or_projection_unavailable"
                continue
            match = sequence_match(query, ref, allow_local=bool(positions), min_identity=min_identity)
            if not match:
                row["reason"] = "no_supported_alignment"
                continue
            row.update(match=match.kind, alignment_identity=match.identity)
            if match.kind == "exact" and set(ref) <= AMINO_ACIDS and int(r["number_target_epitope_residues"]) > 0:
                row["tier"] = "exact_structure"
            else:
                row.update(assess_site(query, ref, match, positions, max_substitutions, min_site_identity, linear=False))
    return rows


class Calibration:
    """Estimate tier precision from independent labelled assays of the same query."""

    def __init__(self, path: Path | None):
        self.assay = ""
        self.groups = defaultdict(dict)
        if path is None:
            return
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = {"query_sha256", "antibody_key", "evidence_tier", "assay", "study_id", "bound"}
            if fields - set(reader.fieldnames or []):
                raise ValueError(f"Calibration file must contain {', '.join(sorted(fields))}")
            for line, r in enumerate(reader, 2):
                if (not all(r[k] and r[k].strip() for k in fields) or r["bound"] not in ("0", "1")
                        or r["evidence_tier"] not in GRADES
                        or not all(re.fullmatch("[0-9a-f]{64}", r[k]) for k in ("query_sha256", "antibody_key"))):
                    raise ValueError(f"Invalid calibration row {line}")
                if self.assay and self.assay != r["assay"]:
                    raise ValueError("Use one assay endpoint and protocol per calibration file")
                self.assay = r["assay"]
                group = self.groups[r["query_sha256"], r["evidence_tier"]]
                outcome, studies = group.setdefault(r["antibody_key"], (int(r["bound"]), set()))
                if outcome != int(r["bound"]):
                    raise ValueError(f"Conflicting calibration labels at row {line}")
                studies.add(r["study_id"])

    def estimate(self, query_sha, tier, key, documents):
        entries = [value for ab, value in self.groups.get((query_sha, tier), {}).items()
                   if ab != key and not value[1].intersection(documents)]
        n = len(entries)
        if not self.assay:
            return {"probability_status": "no_calibration"}
        result = {"calibration_n": n, "calibration_assay": self.assay}
        positives = sum(v[0] for v in entries)
        if n < 10 or positives in (0, n) or len(set().union(*(v[1] for v in entries))) < 2:
            return {**result, "probability_status": "insufficient_matched_assays"}
        return {**result, "probability_status": "assay_calibrated", "binding_probability": (positives + 1) / (n + 2)}


def rank_candidates(rows, query_id, query, calibration):
    groups = defaultdict(list)
    negative_pairs = set()
    for row in rows:
        key = (row["heavy"], row["light"])
        if row["tier"] in GRADES:
            groups[key].append(row)
        elif row.get("reason") == "observed_negative":
            negative_pairs.add(key)
        elif (row.get("reason") == "reported_negative_or_conditional" and row["match"] == "exact"
              and row["audit"] == "CONFIRMED" and not row["flags"]):
            negative_pairs.add(key)
    candidates = []
    query_sha = hashlib.sha256(query.encode()).hexdigest()
    for (heavy, light), supporting in groups.items():
        if (heavy, light) in negative_pairs:
            for row in supporting:
                row.update(tier="review", reason="conflicting_binding_reports")
            continue
        best = max(supporting, key=lambda r: (GRADES[r["tier"]], r["tier"] == "exact_assay", r.get("site_identity", 0)))
        key = pair_id(heavy, light)
        documents = {r["document_id"] for r in supporting}
        evidence_ids = {f"{r['source']}:{r['record_id']}" for r in supporting}
        prediction = calibration.estimate(query_sha, best["tier"], key, documents)
        if best.get("probability_status") == "panel_calibrated" or (
            best.get("probability_status") and "binding_probability" not in prediction
        ):
            prediction = {k: best[k] for k in ("binding_probability", "calibration_assay", "calibration_n", "probability_status")}
        candidates.append({
            "query_id": query_id, "query_sha256": query_sha, "antibody_key": key,
            "heavy": heavy, "light": light, "antibody_names": "|".join(sorted({r["antibody"] for r in supporting if r.get("antibody")})),
            "evidence_tier": best["tier"], "evidence_grade": GRADES[best["tier"]],
            "sources": "|".join(sorted({r["source"] for r in supporting})),
            "best_evidence_record": f"{best['source']}:{best['record_id']}",
            **{k: best[k] for k in ("chain_basis", "site_type", "site_coverage", "site_substitutions", "site_identity") if k in best},
            "source_documents": "|".join(sorted(documents)), "evidence_records": "|".join(sorted(evidence_ids)),
            **prediction,
        })
    def rank(row):
        probability = row.get("binding_probability")
        return (-row["evidence_grade"], -probability if probability is not None else 1,
                -row.get("site_identity", 0), row["antibody_key"])
    return sorted(candidates, key=rank)


def publications(db, candidates, aliases, query_id):
    found = {}
    columns = "p.id,p.heavy,p.light,p.pairing,p.title,p.mentions"
    for name in sorted(aliases):
        for row in db.execute(f"SELECT {columns} FROM mentions m JOIN publications p ON p.id=m.publication_id WHERE m.name=?", (name,)):
            found[row[0]] = ("target_mention_only", row)
    for candidate in candidates:
        for row in db.execute(f"SELECT {columns} FROM publications p WHERE p.pair_key=?", (candidate["antibody_key"],)):
            found[row[0]] = ("exact_chain_pair", row)
    return [dict(zip(PUBLICATION_FIELDS, (query_id, relation, *row)))
            for relation, row in sorted(found.values(), key=lambda item: (item[0], item[1][0]))]


def write_table(path, fields, rows):
    with tempfile.TemporaryDirectory(dir=path.parent) as tempdir:
        staging = Path(tempdir) / path.name
        with staging.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        staging.replace(path)


def search_record(record, db, panel, scores, outdir, args, calibration):
    query, name = record.sequence, record.name
    min_identity = args.min_identity / 100 if args.min_identity > 1 else args.min_identity
    rows, aliases = [], set()
    if db is not None:
        rows, aliases = patent_evidence(db, query, min_identity, args.max_site_substitutions,
                                        args.min_site_identity, args.free_antibody_nm)
    if scores:
        structures = structure_evidence(scores, query, min_identity, args.max_site_substitutions, args.min_site_identity)
        aliases.update(r["target"].casefold() for r in structures if r["tier"] == "exact_structure" and r["target"])
        rows.extend(structures)
    if panel is not None:
        rows.extend(assay_evidence(panel, query, rows))
    for row in rows:
        row["query_id"] = name
    candidates = rank_candidates(rows, name, query, calibration)
    refs = publications(db, candidates, aliases, name) if db is not None else []
    write_table(outdir / "binding_evidence.csv", EVIDENCE_FIELDS, rows)
    write_table(outdir / "binding_candidates.csv", CANDIDATE_FIELDS, candidates)
    write_table(outdir / "plabdab_references.csv", PUBLICATION_FIELDS, refs)
    return candidates
