import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from binding_search import Calibration, bound_fraction, classify_patent, rank_candidates, structure_evidence
from binding_sites import Match, assess_site, sequence_match
from download_pdb import download_pdb_cif
from prepare_hil6 import panel_predictions
from prepare_patents import PLABDAB_FIELDS, POISK_FIELDS, build_catalog, pair_id, patent_record
from run_pipeline import read_query_records, run_single_query


ROOT = Path(__file__).resolve().parents[1]
REF = "ACDEFGHIKLMNPQRSTVWY" + "VYSQMEPFDGIWCKNTLRHA"
HEAVY = "ACDFHIKLMN" * 8
LIGHT = "STVYEGQPRA" * 8


def source_row():
    row = dict.fromkeys(POISK_FIELDS, "")
    row.update({
        "unique_record_id": "000001", "patent_id": "7000000", "antibody_identifier": "control",
        "antibody_type": "Conventional", "vh.sequence_from_patent": HEAVY, "vl.sequence_from_patent": LIGHT,
        "target.sequence": REF, "target.sequence_from_patent": REF, "target.organism": "Human",
        "target.target_type": "Protein", "target.protein_name": "Control protein", "target.gene_names": "CTRL",
        "audit.verdict": "CONFIRMED", "epitope.evidence_types": "Structure Determination",
        "epitope.structural_contact_residues": json.dumps([{"position": p, "residue": REF[p - 1]} for p in (12, 15, 20)]),
    })
    return row


def assess(row, query, limit=0):
    reference, data = patent_record(row, None)
    match = sequence_match(query, reference, allow_local=True, min_identity=0.7)
    return classify_patent(data, query, reference, match, limit, 0.6)


class SiteEvidence(unittest.TestCase):
    def test_fragment_contains_all_contacts(self):
        row = assess(source_row(), REF[7:26])
        self.assertEqual((row["tier"], row["site_covered"], row["site_substitutions"]), ("site_conserved", 3, 0))
        missed = assess(source_row(), REF[16:])
        self.assertEqual((missed["reason"], missed["site_covered"]), ("site_not_fully_covered", 1))

    def test_differences_inside_and_outside(self):
        offsite = REF[:25] + "W" + REF[26:]
        self.assertEqual(assess(source_row(), offsite)["tier"], "site_conserved")
        onsite = REF[:14] + "W" + REF[15:]
        self.assertEqual(assess(source_row(), onsite)["reason"], "site_difference_limit")
        self.assertEqual(assess(source_row(), onsite, 1)["tier"], "site_similar")
        unknown = REF[:14] + "X" + REF[15:]
        self.assertEqual(assess(source_row(), unknown, 1)["reason"], "unknown_site_residue")

    def test_gap_insertion_and_repeated_placement(self):
        segment = REF[10:24]
        ambiguous = sequence_match(segment, REF + REF, allow_local=False, min_identity=.7)
        self.assertTrue(ambiguous.ambiguous)
        self.assertEqual(assess_site(segment, REF + REF, ambiguous, [12], 0, 1, linear=False)["reason"], "ambiguous_alignment")
        inserted = REF[:14] + "WW" + REF[14:]
        match = sequence_match(inserted, REF, allow_local=True, min_identity=.7)
        self.assertEqual(assess_site(inserted, REF, match, list(range(12, 20)), 0, 1, linear=True)["reason"],
                         "insertion_in_mapped_region")
        missing = Match("local", {10: 0, 12: 1}, 1)
        self.assertEqual(assess_site(REF, REF, missing, [11, 12, 13], 2, .5, linear=False)["reason"],
                         "site_not_fully_covered")

    def test_unverified_numbering_is_not_silently_truncated(self):
        for contact in ({"position": 9000, "residue": "T"}, {"position": 1, "residue": "Q"}):
            r = source_row()
            coords = json.loads(r["epitope.structural_contact_residues"]) + [contact]
            r["epitope.structural_contact_residues"] = json.dumps(coords)
            self.assertEqual(assess(r, REF[5:30])["reason"], "coordinate_basis_unverified")
        r = source_row()
        r["epitope.structural_contact_residues"] = ""
        self.assertEqual(assess(r, REF[5:30])["reason"], "site_not_localized")
        r["epitope.linear_epitope_range.start"] = "12"
        r["epitope.linear_epitope_range.end"] = "20"
        self.assertEqual(assess(r, REF[5:30])["reason"], "coordinate_basis_unverified")

    def test_context_and_origin(self):
        for changes, reason in [
            ({"target.sequence_from_patent": "", "target.sequence_from_uniprot": REF}, "assayed_construct_unresolved"),
            ({"target.target_type": "Protein Complex"}, "target_context_required"),
            ({"audit.verdict": "REVIEW"}, "patent_audit_requires_review"),
            ({"target.binding_components": "two subunits"}, "target_context_required"),
            ({"vl.sequence_from_patent": ""}, "incomplete_or_ambiguous_antibody"),
        ]:
            self.assertEqual(assess(source_row() | changes, REF)["reason"], reason)
        vhh = source_row() | {"antibody_type": "VHH", "vl.sequence_from_patent": ""}
        self.assertEqual(assess(vhh, REF)["tier"], "exact_report")

    def test_conflicting_measurement_suppresses_candidate(self):
        yes = assess(source_row(), REF)
        no = assess(source_row() | {"target.binding_affinity.affinity_value": "No detectable binding"}, REF)
        self.assertEqual(rank_candidates([yes, no], "q", REF, Calibration(None)), [])
        self.assertEqual(yes["reason"], "conflicting_binding_reports")

    def test_physical_quantity_requires_a_point_kd(self):
        self.assertAlmostEqual(bound_fraction("Kd = 500 pM", "BLI", 1), 2 / 3)
        for text, method in [("IC50 = 0.5 nM", "SPR"), ("<1 nM", "SPR"), ("1 nM", "ELISA"),
                             ("KD = 5 nM", "FACS"), ("KD = 0 nM", "SPR")]:
            self.assertIsNone(bound_fraction(text, method, 1))

    def test_verified_peptide_requires_the_correct_terminus(self):
        seq = REF[2:12]
        row = source_row() | {"target.sequence": "MKW" + seq + "RK", "target.sequence_from_patent": "MKW" + seq + "RK"}
        verified = dict(patent_id=row["patent_id"], release_target_sha256=hashlib.sha256(row["target.sequence"].encode()).hexdigest(),
                        tested_peptide=seq, mapped_region=seq[:5], free_n_terminus="yes", source="test citation")
        ref, data = patent_record(row, verified)
        self.assertEqual(ref, seq)
        for q, expected in [(seq, "exact_report"), (seq + "WW", "site_conserved"), ("MKW" + seq, "requires_free_n_terminus")]:
            match = sequence_match(q, ref, allow_local=True, min_identity=.7)
            result = classify_patent(data, q, ref, match, 0, .90)
            self.assertEqual(result.get("reason") or result["tier"], expected)
        match = sequence_match(seq[:-1] + "A", seq, allow_local=True, min_identity=.7)
        self.assertEqual((match.kind, match.mapping[9]), ("full_length", 9))

    def test_structure_rechecks_coverage_against_full_reference(self):
        base = dict(pdb_id="test", antigen_id="A", heavy_chain_id="H", light_chain_id="L", antigen_name="control",
                    antigen_species="Homo sapiens", seq_heavy_full=HEAVY, seq_light_full=LIGHT,
                    seq_heavy_anarci=HEAVY, seq_light_anarci=LIGHT, subject_seqres_sequence=REF,
                    number_target_epitope_residues="3", number_multichain_epitope_residues="0",
                    host_contact_positions="12|15|20")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "structure_scores.csv"
            with path.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=base)
                w.writeheader()
                w.writerow(base)
                w.writerow(base | {"number_multichain_epitope_residues": "1"})
            supported, multichain = structure_evidence(path, REF[7:26], .7, 0, .9)
            self.assertEqual((supported["tier"], multichain["reason"]), ("site_conserved", "multichain_site"))
            missed, _ = structure_evidence(path, REF[15:30], .7, 0, .9)
            self.assertEqual(missed["reason"], "site_not_fully_covered")

    def test_corrupt_structure_cache_recovers_before_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "cache"
            cache.mkdir()
            source = root / "source.gz"
            source.write_bytes(gzip.compress(b"# comment\ndata_test\n#\n"))
            (cache / "test.cif.gz").write_bytes(source.read_bytes()[:-5])
            path, status = download_pdb_cif("TEST", cache, source.as_uri(), 5)
            self.assertEqual((path.read_bytes(), status), (source.read_bytes(), "downloaded"))
            self.assertEqual(download_pdb_cif("TEST", cache, source.as_uri(), 5)[1], "cached")


class CatalogAndCLI(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.poisk = self.work / "p.csv.gz"
        self.oxford = self.work / "o.csv"
        self.db = self.work / "patents.sqlite"
        self.write_poisk([source_row(), source_row() | {"unique_record_id": "000002", "target.sequence": ""}])
        with self.oxford.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=PLABDAB_FIELDS)
            writer.writeheader()
            writer.writerow(dict(zip(PLABDAB_FIELDS, ("pair1", HEAVY, LIGHT, "Patent text", "Example", "CTRL"))))
            writer.writerow(dict(zip(PLABDAB_FIELDS, ("pair2", HEAVY, HEAVY, "Unique word", "Example", "CTRL"))))
        build_catalog(self.poisk, self.oxford, self.db)

    def write_poisk(self, rows):
        with gzip.open(self.poisk, "wt", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=POISK_FIELDS)
            writer.writeheader()
            writer.writerows(rows)

    def run_search(self, content):
        fasta = self.work / "input.fa"
        fasta.write_text(content)
        root = self.work / "out"
        subprocess.run([sys.executable, str(ROOT / "run_pipeline.py"), "--sources", "patents", "--query", str(fasta),
                        "--patent-db", str(self.db), "--outdir", str(root)], check=True, capture_output=True)
        return root

    def csv_read(self, path):
        with path.open(encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))

    def test_cli_multi_query_provenance_and_unranked_mentions(self):
        out = self.run_search(f">full\n{REF}\n>partial\n{REF[7:26]}\n>elsewhere\n{REF[25:]}\n")
        rows = self.csv_read(out / "binding_candidates.csv")
        self.assertEqual([(r["query_id"], r["evidence_tier"]) for r in rows], [("full", "exact_report"), ("partial", "site_conserved")])
        self.assertEqual(rows[0]["evidence_records"], "POISK:000001")
        self.assertEqual(rows[0]["binding_probability"], "")
        self.assertEqual(self.csv_read(out / "elsewhere" / "binding_candidates.csv"), [])
        self.assertEqual({r["ID"]: r["relationship"] for r in self.csv_read(out / "full" / "plabdab_references.csv")},
                         {"pair1": "exact_chain_pair", "pair2": "target_mention_only"})
        self.assertEqual(json.loads((out / "run_metadata.json").read_text())["patents"]["counts"]["no_target_or_chain"], 1)

    def test_single_query_rerun_replaces_results(self):
        self.run_search(f">q\n{REF}\n")
        out = self.run_search(">q\nQQQQQQQQQQQQQQQQ\n")
        self.assertEqual(self.csv_read(out / "binding_candidates.csv"), [])
        self.assertEqual(self.csv_read(out / "binding_evidence.csv"), [])

    def test_failed_rebuild_preserves_catalog(self):
        old = self.db.read_bytes()
        self.write_poisk([source_row(), source_row()])
        with self.assertRaises(sqlite3.IntegrityError):
            build_catalog(self.poisk, None, self.db)
        self.assertEqual(self.db.read_bytes(), old)

    def test_curation_file_is_optional_and_explicit_paths_are_checked(self):
        missing = self.work / "missing.tsv"
        with patch("prepare_patents.CURATIONS", missing):
            metadata = build_catalog(self.poisk, None, self.db)
        self.assertEqual(metadata["curations_sha256"], "")
        with self.assertRaises(FileNotFoundError):
            build_catalog(self.poisk, None, self.db, curations=missing)

    def test_input_does_not_join_internal_stops(self):
        fasta = self.work / "q.fa"
        fasta.write_text(">q\nACD*EF\n")
        with self.assertRaisesRegex(ValueError, "Invalid protein"):
            read_query_records(fasta)
        fasta.write_text(">q\n ac dE f*\n")
        self.assertEqual(read_query_records(fasta)[0].sequence, "ACDEF")

    def test_one_failed_structure_download_allows_scoring_to_continue(self):
        fasta = self.work / "q.fa"
        fasta.write_text(f">q\n{REF}\n")
        out = self.work / "s"
        args = SimpleNamespace(force_blast=False, max_target_seqs=10, evalue="1e-5", threads=1,
                               min_identity=.7, max_structures=None, max_score_hits=10,
                               ca_distance_threshold=8., atom_distance_threshold=4.5,
                               homodimer_identity_threshold=.95, scheme="chothia")
        def command(argv, log):
            if argv[0] == "blastp":
                Path(argv[argv.index("-out") + 1]).write_text("")
            elif "download_pdb.py" in argv[1]:
                Path(argv[argv.index("--manifest") + 1]).write_text("pdb_id\tstatus\npdb1\tfailed\n")
                raise subprocess.CalledProcessError(1, argv)
            else:
                Path(argv[argv.index("--output") + 1]).write_text("query_id,pdb_id\n")
        with patch("run_pipeline.run_command", side_effect=command), patch("run_pipeline.shutil.which", return_value="blastp"):
            scores = run_single_query(args, fasta, out, self.work / "s.tsv", self.work / "s.gz",
                                      self.work / "blastdb", self.work / "structures")
        self.assertTrue(scores.exists())
        self.assertIn("failed", (out / "download_manifest.tsv").read_text())


class CalibratedProbabilities(unittest.TestCase):
    def test_external_calibration_deduplicates_and_holds_out_candidate_and_source(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "cal.csv"
            q = "f" * 64
            key = pair_id(HEAVY, LIGHT)
            lines = [f"{q},{i:064x},exact_report,SPR,s{i % 2},{int(i < 7)}" for i in range(10)]
            path.write_text("query_sha256,antibody_key,evidence_tier,assay,study_id,bound\n" +
                            "\n".join(lines + [lines[0], f"{q},{key},exact_report,SPR,other,1",
                                               f"{q},{20:064x},exact_report,SPR,7000000,1"]))
            calibration = Calibration(path)
            p = calibration.estimate(q, "exact_report", key, {"7000000"})
            self.assertEqual(p["calibration_n"], 10)
            self.assertAlmostEqual(p["binding_probability"], 8 / 12)
            self.assertNotIn("binding_probability", calibration.estimate("0" * 64, "exact_report", key, set()))

    def test_panel_holdout_excludes_own_calibration_labels(self):
        profile = np.tile([1, 1, 0, 0, -1], (40, 1))
        fit = np.tile([0, 1], (40, 1))
        p, n, _ = panel_predictions(profile, fit)
        fit[0, :] = 0
        p2, n2, _ = panel_predictions(profile, fit)
        self.assertEqual(n[0], 78)
        self.assertEqual((p[0], n[0]), (p2[0], n2[0]))
        profile[0, :] = -1
        self.assertTrue(np.isnan(panel_predictions(profile, fit)[0][0]))

@unittest.skipUnless(importlib.util.find_spec("biotite") and importlib.util.find_spec("anarci"),
                     "requires structural analysis dependencies")
class StructuralCoordinates(unittest.TestCase):
    def test_compressed_coordinates_preserve_author_insertion_codes(self):
        from biotite.structure import AtomArray
        from biotite.structure.io import pdbx
        from score_structures import chain_residues, load_structure_variants, min_distance

        atoms = AtomArray(2)
        atoms.chain_id[:] = "A"
        atoms.res_id[:] = 42
        atoms.ins_code[:] = ["", "A"]
        atoms.res_name[:] = ["ALA", "GLY"]
        atoms.atom_name[:] = "CA"
        atoms.element[:] = "C"
        atoms.coord = np.array([[0., 0., 0.], [3., 0., 0.]])
        cif = pdbx.CIFFile()
        pdbx.set_structure(cif, atoms)
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "structure.cif.gz"
            with gzip.open(path, "wt") as handle:
                cif.write(handle)
            variant, = load_structure_variants(path)
        self.assertEqual(variant.source, "asymmetric_unit")
        self.assertTrue(variant.warning)
        residues = chain_residues(variant.atoms, "A")
        self.assertEqual([(r.res_id, r.ins_code) for r in residues], [(42, ""), (42, "A")])
        self.assertEqual(min_distance(residues[0].atoms.coord, residues[1].atoms.coord), 3.)


if __name__ == "__main__":
    unittest.main()
