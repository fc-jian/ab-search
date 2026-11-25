"""Optional protein and peptide partner search for a human protein FASTA."""

from __future__ import annotations

import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

from binding_search import write_table
from partner_sources import Cache, UNIPROT, cellphonedb, gtopdb, intact, read_proteome, resolve_query, synbip


INDEXES = {"cellphonedb": cellphonedb, "gtopdb": gtopdb, "synbip": synbip}
SOURCES = ("pdbe", "cellphonedb", "gtopdb", "intact", "synbip")
GRADES = {"exact_structure": 3, "site_conserved": 2, "site_similar": 1}
BASE = [
    "query_id",
    "query_sha256",
    "query_uniprot",
    "gene",
    "query_range",
    "partner_id",
    "partner_name",
    "partner_type",
    "partner_sequences",
    "sequence_basis",
    "status",
    "source",
    "record_id",
]
EVIDENCE_FIELDS = BASE + [
    "reason",
    "query_role",
    "biological_role",
    "assembly_id",
    "target_chains",
    "partner_chains",
    "partner_chain_sequences",
    "complex_context",
    "assembly_components",
    "query_contact_positions",
    "site_size",
    "site_covered",
    "site_coverage",
    "site_substitutions",
    "site_identity",
    "method",
    "affinity",
    "kd_nM",
    "ki_nM",
    "ic50_nM",
    "measurement_relation",
    "assay_parameters",
    "query_feature",
    "partner_feature",
    "modifications",
    "target_modifications",
    "partner_stoichiometry",
    "reference",
    "citation",
    "structure_title",
]
CANDIDATE_FIELDS = BASE + [
    "evidence_records",
    "site_size",
    "site_covered",
    "site_substitutions",
    "partner_chains",
    "partner_chain_sequences",
    "partner_stoichiometry",
    "modifications",
    "assembly_components",
    "citation",
    "form_key",
    "assembly_id",
    "target_chains",
    "query_contact_positions",
    "site_coverage",
    "site_identity",
    "method",
    "reference",
    "structure_title",
    "complex_context",
    "target_modifications",
    "sources",
    "support_assembly_count",
    "support_evidence_json",
    "related_evidence_count",
    "related_evidence_json",
    "reported_measurements_json",
]
BINDER_FIELDS = [
    "query_id", "gene", "query_uniprot", "query_range", "partner_id", "partner_name", "partner_type",
    "partner_sequences", "sequence_basis", "partner_stoichiometry", "modifications",
    "status", "site_size", "site_covered", "site_coverage", "site_substitutions", "site_identity",
    "query_contact_positions", "source", "record_id", "assembly_id", "method", "reference", "citation",
    "structure_title", "target_chains", "partner_chains", "partner_chain_sequences", "complex_context",
    "assembly_components", "target_modifications", "form_key", "supported_form_count",
    "support_assembly_count", "evidence_records", "sources", "related_evidence_count",
    "reported_measurements_json", "supported_forms_json", "related_evidence_json", "query_sha256",
]
SUMMARY_FIELDS = [
    "query_id",
    "query_uniprot",
    "gene",
    "query_range",
    "partner_id",
    "partner_names",
    "best_status",
    "query_supported",
    "supported_forms",
    "sources",
    "evidence_count",
    "reasons",
    "best_citation",
]


def add_options(parser):
    parser.add_argument(
        "--protein-partners", action="store_true", help="Also search published human protein/peptide partners"
    )
    parser.add_argument(
        "--partner-sources",
        nargs="+",
        choices=SOURCES,
        help="Partner sources to search (default: all five when enabled)",
    )
    parser.add_argument(
        "--partner-cache", type=Path, help="Public source snapshots (default: reference/partners)"
    )
    parser.add_argument(
        "--partner-human-fasta",
        type=Path,
        help="Local reviewed human UniProt FASTA including isoforms, plain or gzip",
    )
    parser.add_argument(
        "--partner-offline", action="store_true", help="Use previously cached public partner data"
    )
    parser.add_argument(
        "--partner-pdb",
        nargs="+",
        type=str.lower,
        metavar="ID",
        help="Restrict PDBe to the intersection of these PDB IDs and the discovered target structures",
    )
    parser.add_argument(
        "--partner-max-structures",
        type=int,
        metavar="N",
        help="Inspect at most N distinct mapped PDB entries per query (default: all)",
    )


def structure_rank(row):
    return (-GRADES[row["status"]], int(row["site_substitutions"]), -int(row["site_size"]),
            row["record_id"], row["assembly_id"])


def json_records(rows, fields):
    # Preserve distinct observations of a single accession or assembly, including zero values.
    return json.dumps([{k: r[k] for k in fields if r.get(k) not in (None, "")} for r in rows],
                      ensure_ascii=False, separators=(",", ":"))


OBSERVATION_FIELDS = [
    "source", "record_id", "assembly_id", "status", "reason", "reference", "citation", "method",
    "target_chains", "partner_chains", "site_size", "site_covered", "site_substitutions", "site_coverage",
    "site_identity", "query_contact_positions", "complex_context", "assembly_components",
    "structure_title", "query_role", "biological_role", "sequence_basis", "modifications",
    "target_modifications", "query_feature", "partner_feature", "affinity", "kd_nM", "ki_nM", "ic50_nM",
    "measurement_relation", "assay_parameters",
]
MEASUREMENT_FIELDS = [
    "source", "record_id", "status", "reason", "query_role", "biological_role", "method", "affinity",
    "kd_nM", "ki_nM", "ic50_nM", "measurement_relation", "assay_parameters", "reference", "citation",
]


def ranked(rows):
    groups, related, sources = defaultdict(list), defaultdict(list), defaultdict(set)
    for r in rows:
        if not r.get("partner_id"):
            continue
        pair = (r.get("query_id", ""), r["partner_id"])
        sources[pair].add(r["source"])
        if r["status"] in GRADES and r.get("partner_sequences"):
            form = (r["partner_id"], r["partner_sequences"], r.get("partner_stoichiometry") or "",
                    r.get("modifications") or "")
            groups[pair, form].append(r)
        else:
            related[pair].append(r)
    result = []
    pair_context = {}
    for (pair, form), rs in groups.items():
        if pair not in pair_context:
            context = related[pair]
            measured = [r for r in context if any(r.get(k) not in (None, "")
                        for k in ("affinity", "kd_nM", "ki_nM", "ic50_nM"))]
            pair_context[pair] = dict(sources="|".join(sorted(sources[pair])), related_evidence_count=len(context),
                                      related_evidence_json=json_records(context, OBSERVATION_FIELDS),
                                      reported_measurements_json=json_records(measured, MEASUREMENT_FIELDS))
        best = min(rs, key=structure_rank)
        records = sorted({f"{r['source']}:{r['record_id']}:{r['assembly_id']}" for r in rs})
        key = hashlib.sha256(json.dumps(form, separators=(",", ":")).encode("utf-8")).hexdigest()
        result.append(dict(best, **pair_context[pair], evidence_records="|".join(records),
                           form_key=key, support_assembly_count=len(records),
                           support_evidence_json=json_records(rs, OBSERVATION_FIELDS)))
    return sorted(result, key=lambda r: (r.get("query_id", ""), -GRADES[r["status"]], r["partner_id"], r["partner_sequences"]))


def write_binder_table(outdir, candidates):
    """Write one fragment-supported binder per query and partner ID with its deposited forms."""
    pairs = defaultdict(list)
    for row in candidates:
        pairs[row["query_id"], row["partner_id"]].append(row)
    binders = []
    form_fields = ["form_key", "partner_sequences", "partner_stoichiometry", "modifications", "status",
                   "source", "record_id", "assembly_id", "site_size", "site_covered", "site_substitutions",
                   "evidence_records", "reference", "citation"]
    for forms in pairs.values():
        forms.sort(key=lambda r: (structure_rank(r), r["form_key"]))
        records = sorted({rec for r in forms for rec in r["evidence_records"].split("|")})
        binders.append(dict(forms[0], supported_form_count=len(forms), support_assembly_count=len(records),
                            evidence_records="|".join(records), supported_forms_json=json_records(forms, form_fields)))
    write_table(outdir / "protein_binders.csv", BINDER_FIELDS, binders)


def summarize(rows, candidates):
    groups = defaultdict(list)
    for row in rows:
        if row.get("partner_id"):
            groups[row["partner_id"]].append(row)
    priority = {**GRADES, "structural_context": 0, "parent_only": -1, "reported_negative": -2}
    result = []
    for partner_id, records in groups.items():
        best = max(records, key=lambda r: priority[r["status"]])
        forms = sum(r["partner_id"] == partner_id for r in candidates)
        item = {key: best.get(key, "") for key in ("query_id", "query_uniprot", "gene", "query_range")}
        result.append(
            dict(
                **item,
                partner_id=partner_id,
                best_status=best["status"],
                query_supported=bool(forms),
                supported_forms=forms,
                evidence_count=len(records),
                best_citation=best.get("citation", ""),
                partner_names=" | ".join(sorted({r["partner_name"] for r in records})),
                reasons="|".join(sorted({r["reason"] for r in records if r.get("reason")})),
                sources="|".join(sorted({r["source"] for r in records})),
            )
        )
    return sorted(result, key=lambda r: (-priority[r["best_status"]], r["partner_id"]))


class PartnerSearch:
    def __init__(self, args):
        self.args = args
        self.sources = list(dict.fromkeys(args.partner_sources or SOURCES))
        self.cache = Cache(args.partner_cache or args.reference_dir / "partners", args.partner_offline)
        if args.partner_human_fasta:
            data = args.partner_human_fasta.read_bytes()
            self.proteins = read_proteome(data)
            self.cache.record("human.fasta", str(args.partner_human_fasta), data)
        else:
            self.proteins = self.cache.read("human-reviewed-isoforms.fasta.gz", UNIPROT, read_proteome)
        self.indexes = {}
        self.queries = []

    def structures(self, parent, query, match, report):
        from partner_structures import load_entry, mapped_structures, structure_rows

        ids = mapped_structures(self.cache, parent.accession, match.mapping)
        if parent.canonical != parent.accession:
            ids = list(dict.fromkeys(ids + mapped_structures(self.cache, parent.canonical)))
        report["discovered"] = len(ids)
        if self.args.partner_pdb:
            allowed = set(self.args.partner_pdb)
            ids = [i for i in ids if i in allowed]
        if self.args.partner_max_structures:
            ids = ids[: self.args.partner_max_structures]
        report["selected"] = ids
        rows, errors = [], []
        for pdb_id in ids:
            try:
                block = load_entry(self.cache, pdb_id)
                rows.extend(
                    structure_rows(
                        block,
                        pdb_id,
                        self.proteins,
                        parent,
                        query,
                        match,
                        self.args.atom_distance_threshold,
                        self.args.max_site_substitutions,
                        self.args.min_site_identity,
                    )
                )
            except Exception as exc:
                errors.append(pdb_id)
                rows.append(
                    dict(
                        source="PDBe",
                        record_id=pdb_id,
                        status="source_error",
                        reason=f"{type(exc).__name__}: {exc}",
                    )
                )
        if errors:
            report.update(status="partial" if len(errors) < len(ids) else "failed", failed_entries=errors)
        return rows

    def run(self, record, outdir):
        query = record.sequence
        parent, match, accessions, reason = resolve_query(query, self.proteins)
        report = dict(
            query_id=record.name, matching_accessions=accessions, resolution=reason or match.kind, sources={}
        )
        info = dict(query_id=record.name, query_sha256=hashlib.sha256(query.encode()).hexdigest())
        rows = []
        if parent:
            info.update(
                query_uniprot=parent.accession,
                gene=parent.gene,
                query_range=f"{min(match.mapping) + 1}-{max(match.mapping) + 1}",
            )
            for name in self.sources:
                state = {"status": "ok"}
                try:
                    if name == "pdbe":
                        found = self.structures(parent, query, match, state)
                    elif name == "intact":
                        found = intact(self.cache, parent.canonical, self.proteins)
                    else:
                        if name not in self.indexes:
                            self.indexes[name] = INDEXES[name](self.cache, self.proteins)
                        found = self.indexes[name].get(parent.canonical, [])
                    rows.extend(found)
                    state["records"] = len(found)
                except Exception as exc:
                    state.update(status="failed", error=f"{type(exc).__name__}: {exc}")
                    rows.append(dict(source=name, status="source_error", reason=state["error"]))
                if state["status"] != "ok":
                    print(
                        f"{record.name}: {name} {state['status']}: {state.get('error', state.get('failed_entries'))}",
                        file=sys.stderr,
                    )
                report["sources"][name] = state
        else:
            rows.append(
                dict(source="UniProt", status="query_review", reason=reason, reference="|".join(accessions))
            )
        rows = [dict(r, **info) for r in rows]
        candidates = ranked(rows)
        summary = summarize(rows, candidates)
        write_table(outdir / "partner_evidence.csv", EVIDENCE_FIELDS, rows)
        write_table(outdir / "partner_candidates.csv", CANDIDATE_FIELDS, candidates)
        write_binder_table(outdir, candidates)
        write_table(outdir / "partner_summary.csv", SUMMARY_FIELDS, summary)
        states = [s["status"] for s in report["sources"].values()]
        outcome = "unresolved" if not parent else "complete"
        if any(s != "ok" for s in states):
            outcome = "failed" if all(s == "failed" for s in states) else "partial"
        report.update(
            status=outcome,
            candidate_count=len(candidates),
            evidence_count=len(rows),
            reported_partners=len(summary),
        )
        self.queries.append(report)
        return candidates, rows, summary

    def metadata(self):
        return dict(
            sources=self.sources,
            offline=self.cache.offline,
            organism="Homo sapiens (9606)",
            pdb_filter=self.args.partner_pdb,
            max_structures=self.args.partner_max_structures,
            contact_angstrom=self.args.atom_distance_threshold,
            snapshots=self.cache.used,
            queries=self.queries,
        )
