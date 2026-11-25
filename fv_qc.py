"""Evidence-preserving IMGT variable-domain checks on reported antibody sequences."""
from __future__ import annotations

import csv
import os
import re
import shutil
import sys
import time
import urllib.parse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from binding_sites import AMINO_ACIDS
from partner_sources import Cache
from prepare_patents import pair_id


REFERENCES = Path(__file__).resolve().parent / "data" / "verified_fv_sequences.tsv"
CHAIN_COLUMNS = ("status", "input_issues", "input_imgt_span", "input_span_1based",
                 "input_cdr3_imgt", "sequence", "cdr3_imgt", "source", "reference_checks")
FV_FIELDS = ["fv_status", "fv_format", "fv_scheme", "fv_warnings", "fv_pair_sha256"] + [
    f"fv_{chain}_{col}" for chain in ("heavy", "light") for col in CHAIN_COLUMNS
]
SINGLE_DOMAIN = re.compile(r"\b(?:VHH|nanobody|single[- ]domain)\b", re.I)
TAG_OR_LINKER = re.compile(r"H{6,}|(?:GGGGS){2,}")


def structure_context(row):
    """Recover raw SEQRES only when the score carries metadata from that author chain."""
    heavy = row.get("seq_heavy_full", "") if row.get("heavy_chain_description") else ""
    light = row.get("seq_light_full", "") if row.get("light_chain_description") else ""
    if not light and not row.get("light_chain_id") and row.get("seq_light_anarci"):
        form = "scFv"
    elif not row.get("light_chain_id") and SINGLE_DOMAIN.search(row.get("heavy_chain_description", "")):
        form = "single_domain"
    else:
        form = "paired" if row.get("light_chain_id") else ""
    return {"antibody_format": form, "_deposited_heavy": heavy,
            "_deposited_light": heavy if form == "scFv" else light}


@dataclass(frozen=True)
class Domain:
    kind: str
    sequence: str
    span: tuple[int, int]  # zero-based, inclusive, on the input chain
    imgt: tuple[int, int]
    cdr3: str
    problems: tuple[str, ...]


def read_domain(raw, result, detail, j_genes, occupied_v):
    numbered, start, end = result
    residues = [(pos, aa) for (pos, _ins), aa in numbered if aa != "-"]
    # Keep the supplied order of all IMGT insertion positions, especially in CDR3.
    seq = "".join(aa for _, aa in residues)
    if not seq or start is None or end is None or raw[start:end + 1] != seq:
        raise ValueError("numbering_coordinate_mismatch")
    kind = detail.get("chain_type", "")
    expected_end = 128 if kind == "H" else 127
    if kind in {"K", "L"}:
        # A following constant-domain residue can be aligned to HMM position 128.
        residues = [(p, aa) for p, aa in residues if p <= 127]
        seq = "".join(aa for _, aa in residues)
        end = start + len(seq) - 1
    if not residues or raw[start:end + 1] != seq:
        raise ValueError("numbering_coordinate_mismatch")
    first, last = residues[0][0], residues[-1][0]
    points = dict(residues)
    def region(a, b):
        return "".join(aa for p, aa in residues if a <= p <= b)
    problems = []
    if first != 1:
        problems.append(f"n_boundary_starts_at_{first}")
    if last < expected_end:
        # A few alleles have a genuinely shorter J region. Require its complete,
        # exact FR4 in the HMM-assigned species before accepting that endpoint.
        fr4 = region(118, 128)
        species_alleles = j_genes.get(kind, {})
        alleles = species_alleles.get(detail.get("species"), {}).values()
        known_short = any(len(j.rstrip("-")) == last and j[117:].replace("-", "") == fr4 for j in alleles)
        could_be_truncated = any(
            j[117:].replace("-", "").startswith(fr4) and len(j.rstrip("-")) > last
            for variants in species_alleles.values() for j in variants.values()
        )
        if not known_short or could_be_truncated:
            problems.append(f"c_boundary_ends_at_{last}")
    if not set(seq) <= AMINO_ACIDS:
        problems.append("ambiguous_domain_residue")
    if any(m.start() <= start < m.end() or m.start() <= end < m.end() for m in TAG_OR_LINKER.finditer(raw)):
        problems.append("tag_or_linker_overlaps_boundary")
    # Ordinary IMGT alignment gaps are common. Flag a substantial internal FR
    # deletion only at positions occupied by every available V allele of this type.
    expected = occupied_v.get((kind, detail.get("species")), set())
    for a, b in ((1, 26), (39, 55), (66, 104)):
        missing = "".join("1" if first < p < last and p in expected and p not in points else "0"
                          for p in range(a, b + 1))
        for gap in re.finditer(r"1{3,}", missing):
            problems.append(f"framework_gap_{a + gap.start()}_{a + gap.end() - 1}")
    for a, b, label in ((27, 38, "cdr1"), (56, 65, "cdr2"), (105, 117, "cdr3")):
        if not region(a, b):
            problems.append(f"{label}_missing")
    if 104 not in points or 118 not in points:
        problems.append("cdr3_boundary_unresolved")
    for pos, allowed in ((23, "C"), (41, "W"), (104, "C"), (118, "FW")):
        if pos in points and points[pos] not in allowed:
            problems.append(f"atypical_anchor_{pos}_{points[pos]}")
    if any(p not in points for p in range(118, last + 1)):
        problems.append("fr4_internal_gap")
    return Domain(kind, seq, (int(start), int(end)), (first, last), region(105, 117), tuple(problems))


def supported_extension(raw, full, *, internal=False):
    """Allow exact terminal containment or one documented internal gap with unique flanks."""
    if len(raw) < 20 or len(full) <= len(raw):
        return False
    if raw in full:
        return full.find(raw) == full.rfind(raw)
    if internal:
        missing = len(full) - len(raw)
        for i in range(8, len(raw) - 7):
            if (full.startswith(raw[:i]) and full[i + missing:] == raw[i:]
                    and full.rfind(raw[:i]) == 0 and full.find(raw[i:]) == i + missing):
                return True
    return False


class FvChecker:
    def __init__(self, references=None, threads=1, patent_cache=None, offline=False):
        self.threads = max(1, threads)
        self.cache = {}
        self.patents = Cache(patent_cache, offline) if patent_cache is not None else None
        self.listings = {}
        self.last_request = 0.0
        self.references = defaultdict(list)
        if references is None and REFERENCES.is_file():
            references = REFERENCES
        self.reference_path = Path(references) if references is not None else None
        if self.reference_path is not None:
            with self.reference_path.open(encoding="utf-8-sig", newline="") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    if (row["chain"] not in {"heavy", "light"} or not row["citation"] or not row["reference"]
                            or row["relationship"] not in {"same_record", "same_clone"}
                            or not re.fullmatch(r"[0-9a-f]{64}", row["raw_pair_sha256"])
                            or not row["sequence"] or not set(row["sequence"]) <= AMINO_ACIDS):
                        raise ValueError(f"Invalid Fv reference for {row['record_id']}")
                    self.references[row["record_id"]].append(row)
        self.engine_error = ""
        try:
            from anarci import anarci
            from anarci.germlines import all_germlines
            self.number, self.j_genes = anarci, all_germlines["J"]
            self.occupied_v = {
                (kind, species): {i for i, aas in enumerate(zip(*alleles.values()), 1) if "-" not in aas}
                for kind, species_map in all_germlines["V"].items() if kind in {"H", "K", "L"}
                for species, alleles in species_map.items()
            }
        except (ImportError, OSError):
            self.engine_error = "anarci_unavailable"
        bindir = Path(sys.executable).parent
        binary = shutil.which("hmmscan") or shutil.which(
            "hmmscan", path=os.pathsep.join(str(p) for p in (bindir, bindir / "Library" / "bin")))
        self.hmmerpath = str(Path(binary).parent) if binary else ""
        if not binary and not self.engine_error:
            self.engine_error = "hmmscan_unavailable"

    def patent_listing(self, patent, seq_id):
        """Fetch one protein SEQ ID, verifying both patent number and sequence identifier."""
        if not all(re.fullmatch(r"[1-9]\d*", value) for value in (patent, seq_id)):
            raise ValueError("Expected numeric US patent and protein SEQ ID")
        key = (patent, seq_id)
        if key in self.listings:
            return self.listings[key]
        name = f"US{patent}_{seq_id}.fasta"
        source = f"US{patent} SEQ ID {seq_id}"
        params = {"db": "protein", "id": f"pat|US|{patent}|{seq_id}", "rettype": "fasta", "retmode": "text"}
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
        def parse(data):
            lines = data.decode("ascii").strip().splitlines()
            head = re.fullmatch(r">([\w]+\.\d+) Sequence (\d+) from patent US (\d+)\s*", lines[0] if lines else "")
            seq = "".join(lines[1:])
            if not head or head.group(2, 3) != (seq_id, patent) or not seq or not set(seq) <= AMINO_ACIDS | {"X"}:
                raise ValueError("patent FASTA identifier or sequence mismatch")
            return seq, f"{source} https://www.ncbi.nlm.nih.gov/protein/{head[1]}"
        if not self.patents.offline and not (self.patents.root / name).is_file():
            time.sleep(max(0.0, .36 - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()
        try:
            result = self.patents.read(name, url, parse)
        except (OSError, ValueError) as exc:
            # A failed lookup must not turn an incomplete sequence into a pass.
            result = ("", f"{source} lookup_unavailable ({type(exc).__name__})")
        self.listings[key] = result
        return result

    def prime(self, sequences):
        pending = []
        for seq in sorted(set(sequences) - self.cache.keys()):
            if not seq:
                self.cache[seq] = ([], "missing_chain")
            elif not set(seq) <= AMINO_ACIDS | {"X"}:
                self.cache[seq] = ([], "invalid_chain_characters")
            elif self.engine_error:
                self.cache[seq] = ([], self.engine_error)
            else:
                pending.append(seq)
        for start in range(0, len(pending), 256):
            self._number_batch(pending[start:start + 256])

    def _number_batch(self, seqs):
        try:
            results, details, _ = self.number([(str(i), s) for i, s in enumerate(seqs)],
                                               scheme="imgt", output=False, ncpu=self.threads,
                                               hmmerpath=self.hmmerpath)
            if len(results) != len(seqs) or len(details) != len(seqs):
                raise ValueError("incomplete ANARCI batch")
        except Exception:
            if len(seqs) > 1:
                mid = len(seqs) // 2
                self._number_batch(seqs[:mid])
                self._number_batch(seqs[mid:])
            else:
                self.cache[seqs[0]] = ([], "anarci_failed")
            return
        for seq, assignments, info in zip(seqs, results, details):
            if not assignments:
                self.cache[seq] = ([], "no_numbered_domain")
                continue
            try:
                domains = [read_domain(seq, item, meta, self.j_genes, self.occupied_v)
                           for item, meta in zip(assignments, info, strict=True)]
                self.cache[seq] = (domains, "")
            except (IndexError, KeyError, ValueError, TypeError, AttributeError):
                self.cache[seq] = ([], "numbering_coordinate_mismatch")

    def select(self, seq, chain, scfv=False):
        domains, error = self.cache[seq]
        if error:
            return None, error
        if len(domains) > 1:
            if not scfv or sorted(d.kind for d in domains) not in (["H", "K"], ["H", "L"]):
                return None, "multiple_domains"
        expected = {"H"} if chain == "heavy" else {"K", "L"}
        hits = [d for d in domains if d.kind in expected]
        if len(hits) != 1:
            return None, "unexpected_chain_type"
        return hits[0], ""

    def chain_result(self, raw, chain, scfv, donors):
        domain, error = self.select(raw, chain, scfv)
        out = {"status": "warning"}
        issues = [error] if error else list(domain.problems)
        if domain:
            lo, hi = domain.span
            out.update(input_imgt_span=f"{domain.imgt[0]}-{domain.imgt[1]}",
                       input_span_1based=f"{lo + 1}-{hi + 1}", input_cdr3_imgt=domain.cdr3)
            if lo:
                issues.append(f"n_outside_numbering_{lo}")
            if hi < len(raw) - 1:
                issues.append(f"c_outside_numbering_{len(raw) - hi - 1}")
        out["input_issues"] = "|".join(issues)
        if domain and not domain.problems:
            # A deletion inside CDR3 can leave both anchors and a plausible
            # shorter loop. Resolve it only against the actual source chain.
            originals = defaultdict(list)
            for seq, rel, source, is_fusion in donors:
                if rel == "same_record" and seq:
                    deposited, _ = self.select(seq, chain, is_fusion)
                    if deposited and not deposited.problems and deposited.kind == domain.kind:
                        originals[deposited.sequence, deposited.cdr3].append(source)
            if any(seq != domain.sequence for seq, _ in originals):
                issues.append("same_record_domain_disagrees")
                out.update(input_issues="|".join(issues), reference_checks=" | ".join(
                    sorted({src for refs in originals.values() for src in refs})))
                if len(originals) == 1:
                    (seq, cdr3), refs = next(iter(originals.items()))
                    if supported_extension(domain.sequence, seq, internal=True):
                        out.update(status="recovered_same_record", sequence=seq, cdr3_imgt=cdr3,
                                   source=" | ".join(sorted(set(refs))))
                return out
            out.update(sequence=domain.sequence, cdr3_imgt=domain.cdr3,
                       status="pass" if domain.sequence == raw else "trimmed", source="reported_chain")
            return out
        if not raw or error in {"multiple_domains", "unexpected_chain_type", "invalid_chain_characters"}:
            return out
        options = defaultdict(list)
        checks = []
        for seq, rel, source, is_fusion in donors:
            if not seq:
                checks.append(source)
                continue
            donor, _ = self.select(seq, chain, is_fusion)
            if not donor or donor.problems or (domain and donor.kind != domain.kind):
                checks.append(source + " [no verified complete domain]")
                continue
            # A truncated patent listing may still carry part of its leader.
            # Accept the aligned intersection only if the *entire* raw input occurs
            # unchanged at a unique location in the very same deposited sequence.
            placed = seq.find(raw)
            overlap = min(placed + len(raw), donor.span[1] + 1) - max(placed, donor.span[0])
            same_construct = (rel == "same_record" and placed >= 0 and placed == seq.rfind(raw)
                              and overlap >= 20 and overlap < len(donor.sequence))
            # An explicit terminal His tag can mask a truncated FR4 by filling
            # empty HMM positions. Compare the intact remainder to the same record.
            untagged = re.sub(r"^H{6,}|H{6,}$", "", raw) if rel == "same_record" else raw
            if same_construct or supported_extension(untagged, donor.sequence, internal=rel == "same_record"):
                options[donor.sequence, donor.cdr3].append((rel, source))
                checks.append(source + " [verified overlap]")
            else:
                checks.append(source + " [reported residues do not match]")
        out["reference_checks"] = " | ".join(sorted(set(checks)))
        if len(options) == 1:
            (sequence, cdr3), refs = next(iter(options.items()))
            status = "recovered_same_record" if any(rel == "same_record" for rel, _ in refs) else "same_clone_reference"
            out.update(status=status, sequence=sequence, cdr3_imgt=cdr3,
                       source=" | ".join(sorted({source for _, source in refs})))
        elif len(options) > 1:
            out["input_issues"] += "|conflicting_repair_sources"
        return out

    def annotate(self, candidates, evidence):
        by_id = {r["source"] + ":" + r["record_id"]: r for r in evidence}
        tasks = []
        seqs = set()
        for row in candidates:
            heavy, light = row["heavy"], row["light"]
            raw_key = pair_id(heavy, light)
            if row["antibody_key"] != raw_key:
                raise ValueError(f"Fv input sequence key mismatch for {row['query_id']}")
            donors = {"heavy": [], "light": []}
            patent_ids = defaultdict(set)
            formats, warnings = set(), []
            for rid in row["evidence_records"].split("|"):
                e = by_id.get(rid, {})
                if e and (e["heavy"], e["light"]) != (heavy, light):
                    warnings.append("evidence_pair_mismatch:" + rid)
                    continue
                if e.get("antibody_format"):
                    formats.add(e["antibody_format"])
                for chain in donors:
                    if (e.get("source") == "POISK" and re.fullmatch(r"[1-9]\d*", e.get("document_id", ""))
                            and re.fullmatch(r"[1-9]\d*", e.get(chain + "_seq_id", ""))):
                        patent_ids[chain].add((e["document_id"], e[chain + "_seq_id"]))
                    seq = e.get("_deposited_" + chain, "")
                    if seq:
                        donors[chain].append((seq, "same_record", rid + " SEQRES", e.get("antibody_format") == "scFv"))
                for ref in self.references.get(rid, ()):
                    if ref["raw_pair_sha256"] != raw_key:
                        warnings.append("source_reference_key_mismatch:" + rid)
                        continue
                    donors[ref["chain"]].append((ref["sequence"], ref["relationship"],
                                                 ref["reference"] + " " + ref["citation"], False))
            single = bool(formats & {"VHH", "single_domain"})
            fusion = "scFv" in formats
            if single and (light or fusion or formats & {"paired", "Conventional"}):
                warnings.append("conflicting_antibody_format")
            form = "scFv" if fusion else "single_domain" if single else "paired" if light else "unconfirmed_pair"
            if fusion and not light and heavy:
                light = heavy
            tasks.append((row, heavy, light, form, donors, patent_ids, warnings))
            seqs.update((heavy, light))
            seqs.update(seq for choices in donors.values() for seq, *_ in choices)
        self.prime(seqs)
        extra = set()
        if self.patents and not self.engine_error:
            for _, heavy, light, form, donors, patent_ids, _ in tasks:
                for chain, raw in (("heavy", heavy), ("light", light)):
                    domain, error = self.select(raw, chain, form == "scFv")
                    if error == "no_numbered_domain" or (domain and (domain.problems or len(domain.cdr3) <= 4)):
                        for pat, seq_id in sorted(patent_ids[chain]):
                            seq, source = self.patent_listing(pat, seq_id)
                            donors[chain].append((seq, "same_record", source, form == "scFv"))
                            extra.add(seq)
            self.prime(extra)
        for row, heavy, light, form, donors, _, warnings in tasks:
            statuses, checked = [], []
            for chain, raw in (("heavy", heavy), ("light", light)):
                if chain == "light" and form == "single_domain" and not raw:
                    result = {"status": "not_applicable"}
                else:
                    result = self.chain_result(raw, chain, form == "scFv", donors[chain])
                statuses.append(result["status"])
                checked.append(result.get("sequence", ""))
                if result["status"] == "warning":
                    warnings.extend(f"{chain}:{why}" for why in result["input_issues"].split("|") if why)
                row.update({f"fv_{chain}_{col}": result.get(col, "") for col in CHAIN_COLUMNS})
            status = "pass"
            if warnings:
                status = "warning"
            elif any(s in {"recovered_same_record", "same_clone_reference"} for s in statuses):
                status = "recovered"
            elif "trimmed" in statuses:
                status = "trimmed"
            row.update(fv_format=form, fv_scheme="imgt", fv_status=status,
                       fv_warnings="|".join(dict.fromkeys(warnings)),
                       fv_pair_sha256=pair_id(*checked) if status != "warning" else "")
