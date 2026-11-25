import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fv_qc import FvChecker, read_domain, supported_extension
from prepare_patents import pair_id


VH = "ACDEFGHIKLMNPQRSTVWY" + "YVTSRQPNMLKIHGFEDCWA"
VL = VH[::-1]


def candidate(h, l, rid="SAbDab:test", fmt="paired", **extra):
    source, record = rid.split(":", 1)
    row = dict(query_id="T", heavy=h, light=l, antibody_key=pair_id(h, l), evidence_records=rid)
    evidence = dict(source=source, record_id=record, heavy=h, light=l, antibody_format=fmt, **extra)
    return row, evidence


class FvRecords(unittest.TestCase):
    def test_reference_table_is_optional_and_explicit_paths_are_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "references.tsv"
            with patch("fv_qc.REFERENCES", path):
                self.assertIsNone(FvChecker().reference_path)
            with self.assertRaises(FileNotFoundError):
                FvChecker(path)
            path.write_text("record_id\traw_pair_sha256\tchain\tsequence\trelationship\treference\tcitation\n"
                            f"source:1\t{pair_id(VH, VL)}\theavy\t{VH}\tsame_record\ttest\tcitation\n")
            checker = FvChecker(path)
            self.assertEqual(checker.reference_path, path)
            self.assertEqual(checker.references["source:1"][0]["sequence"], VH)

    def test_fetch_requires_exact_patent_and_seq_id_and_uses_disk_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            checker = FvChecker(patent_cache=tmp, offline=True)
            path = Path(tmp) / "US12345_6.fasta"
            for header in ("ABC123.1 Sequence 5 from patent US 12345", "ABC123.1 Sequence 6 from patent US 24567"):
                path.write_text(">" + header + "\n" + VH + "\n")
                checker.listings.clear()
                self.assertEqual(checker.patent_listing("12345", "6")[0], "")
            path.write_text(">ABC123.1 Sequence 6 from patent US 12345\n" + VH + "\n")
            checker.listings.clear()
            seq, source = checker.patent_listing("12345", "6")
            self.assertEqual(seq, VH)
            self.assertIn("/protein/ABC123.1", source)

    def test_engine_absent_leaves_warning_and_original_pair(self):
        checker = FvChecker()
        checker.engine_error = "anarci_unavailable"
        row, evidence = candidate(VH, VL)
        checker.annotate([row], [evidence])
        self.assertEqual((row["fv_status"], row["fv_pair_sha256"], row["heavy"]), ("warning", "", VH))
        self.assertIn("anarci_unavailable", row["fv_warnings"])

    def test_reference_overlap_rejects_substitutions_and_unrelated_chains(self):
        self.assertTrue(supported_extension(VH[4:], VH))
        self.assertFalse(supported_extension("W" + VH[4:], VH))
        deleted = VH[:14] + VH[20:]
        self.assertTrue(supported_extension(deleted, VH, internal=True))
        self.assertFalse(supported_extension(deleted, VH))
        self.assertFalse(supported_extension("WGQG", VH))
        self.assertFalse(supported_extension("A" * 20, "A" * 21))

    def test_short_j_allele_must_be_distinguishable_from_truncation(self):
        fr4 = "WGQGTLVTVS"
        seq = "A" * 117 + fr4
        numbering = ([((i, " "), aa) for i, aa in enumerate(seq, 1)], 0, len(seq) - 1)
        genes = {"H": {"alpaca": {"short": "-" * 117 + fr4 + "-"},
                       "human": {"long": "-" * 117 + fr4 + "S"}}}
        domain = read_domain(seq, numbering, {"chain_type": "H", "species": "alpaca"}, genes, {})
        self.assertIn("c_boundary_ends_at_127", domain.problems)

if __name__ == "__main__":
    unittest.main()
