import csv
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from binding_sites import Match
from partner_sources import Cache, Protein, cellphonedb, gtopdb, intact, read_proteome, resolve_query, synbip
from protein_partners import ranked
from run_pipeline import read_query_records


ROOT = Path(__file__).resolve().parents[1]
SEQ = "MLRCSHAFDGTEPVKIQNWFYPLDEKMAAT"
OTHER = "TWCDERALGMNHIVKPEFQWYTASLC"


def table(path, fields, rows, banner=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as out:
        if banner:
            out.write(f'"# GtoPdb Version: {banner}"\n')
        writer = csv.DictWriter(out, fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


class PartnerSources(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.cache = Cache(self.root, offline=True)
        self.fasta = self.root / "reference.fa"
        self.fasta.write_text(
            f">sp|P00001|ONE_HUMAN Receptor OS=Homo sapiens OX=9606 GN=ONE\n{SEQ}\n"
            f">sp|P00002|TWO_HUMAN Hormone OS=Homo sapiens OX=9606 GN=TWO\n{OTHER}\n"
            f">sp|Q00001|EXAMPLE_MOUSE Other protein OS=Mus musculus OX=10090\n{SEQ}\n"
        )
        self.proteins = read_proteome(self.fasta.read_bytes())
        cp = self.root / "cellphonedb/v5.0.0"
        table(
            cp / "complex_input.csv",
            ["complex_name", "uniprot_1", "uniprot_2"],
            [{"complex_name": "R", "uniprot_1": "P00001", "uniprot_2": "P00002"}],
        )
        cols = ["partner_a", "partner_b", "is_ppi", "source", "directionality", "annotation_strategy"]
        rows = [("P00002", "R", "True"), ("X", "P00002", "False"), ("P00002", "X", "True")]
        table(
            cp / "interaction_input.csv",
            cols,
            [
                dict(
                    partner_a=a,
                    partner_b=b,
                    is_ppi=ppi,
                    source="PMID:1",
                    directionality="Ligand-Receptor",
                    annotation_strategy="curated",
                )
                for a, b, ppi in rows
            ],
        )

    def test_identity_fragments_variants_and_ambiguity(self):
        self.assertNotIn("Q00001", self.proteins)
        parent, m, ids, error = resolve_query(SEQ[2:22], self.proteins)
        self.assertEqual((parent.accession, m.mapping[2], ids, error), ("P00001", 0, ["P00001"], ""))
        query = "W" + SEQ[2:27]
        parent, m, _, _ = resolve_query(query, self.proteins)
        self.assertEqual((parent.accession, m.kind, m.mapping[1]), ("P00001", "variant", 0))
        repeated = {"A": Protein("A", "A", "protein", SEQ + SEQ)}
        self.assertEqual(resolve_query(SEQ, repeated)[-1], "ambiguous_query_position")
        same = self.proteins | {"P3": Protein("P3", "THREE", "protein", SEQ)}
        self.assertEqual(resolve_query(SEQ, same)[-1], "ambiguous_human_parent")

    def test_cellphone_both_directions_and_complex(self):
        result = cellphonedb(self.cache, self.proteins)
        (receptor,) = result["P00001"]
        self.assertEqual(
            (receptor["query_role"], receptor["reason"], receptor["partner_id"]),
            ("Receptor", "query_complex_subunit", "P00002"),
        )
        ligand = next(r for r in result["P00002"] if r["query_role"] == "Ligand")
        self.assertEqual((ligand["partner_type"], ligand["partner_id"]), ("protein_complex", "P00001|P00002"))
        self.assertTrue(all(r["status"] == "parent_only" for rows in result.values() for r in rows))

    def test_gtopdb_uses_peptides_and_retains_mature_sequence(self):
        d = self.root / "gtopdb"
        peptides = [
            "Ligand id",
            "Name",
            "Species",
            "Single letter amino acid sequence",
            "Subunit ids",
            "UniProt id",
        ]
        table(
            d / "peptides.csv",
            peptides,
            [dict(zip(peptides, ["12", "peptide", "Human", OTHER[4:12], "", "P00002"]))],
            "1",
        )
        endo = ["Ligand ID", "Ligand UniProt ID", "Target UniProt ID", "Target ID"]
        table(
            d / "endogenous_ligands.csv",
            endo,
            [dict(zip(endo, ["12", "P00002", "P00001", "1"])), dict(zip(endo, ["99", "", "P00001", "1"]))],
            "1",
        )
        fields = [
            "Ligand ID",
            "Ligand Type",
            "Target Species",
            "Target UniProt ID",
            "Target ID",
            "PubMed ID",
            "Affinity Units",
        ]
        table(
            d / "interactions.csv",
            fields,
            [
                dict(zip(fields, ["12", "Peptide", "Human", "P00001", "1", "200", "pKd"])),
                dict(zip(fields, ["99", "Synthetic organic", "Human", "P00001", "1", "201", "pKi"])),
            ],
            "1",
        )
        records = gtopdb(self.cache, self.proteins)
        self.assertEqual({r["partner_sequences"] for r in records["P00001"]}, {OTHER[4:12]})
        self.assertEqual({r["partner_id"] for r in records["P00002"]}, {"P00001"})
        self.assertEqual(ranked(records["P00001"]), [])
        table(d / "interactions.csv", fields, [], "2")
        with self.assertRaisesRegex(ValueError, "different releases"):
            gtopdb(self.cache, self.proteins)

    def test_intact_preserves_negative_and_expanded_evidence(self):
        d = self.root / "intact"
        d.mkdir()
        base = ["-"] * 42
        base[:2] = ["uniprotkb:P00002", "uniprotkb:P00001"]
        base[6], base[8] = "SPR", "pubmed:123"
        base[9:11] = ["taxid:9606(human)"] * 2
        base[11], base[13] = 'psi-mi:"MI:0407"(direct interaction)', "intact:EBI-123"
        base[36:38] = ["ligand feature", "target feature"]
        neg = base.copy()
        neg[35] = "true"
        expanded = base.copy()
        expanded[15] = "spoke expansion"
        (d / "P00001.0.tsv").write_text("\n".join("\t".join(r) for r in [base, neg, expanded]))
        result = intact(self.cache, "P00001", self.proteins)
        self.assertEqual(
            [r["reason"] for r in result], ["direct_interaction", "negative_assay", "association_or_complex"]
        )
        self.assertEqual(result[0]["query_feature"], "target feature")
        self.assertEqual(ranked(result), [])

    def test_synbip_sequence_and_affinity_are_parent_level(self):
        d = self.root / "synbip"
        fields = ["SBPID", "Protein Scaffold Name"]
        table(
            d / "General_Information_of_SBP.csv",
            fields,
            [dict(zip(fields, ["SBP1", "Peptide aptamer"])), dict(zip(fields, ["SBP2", "Fab"]))],
        )
        cols = ["SBP ID", "SBP Name", "Affinity", "Target Uniprot ID", "TargetID"]
        table(
            d / "Target_to_SBP_Mapping_with_Mechanism_of_Action.csv",
            cols,
            [
                dict(zip(cols, [key, "reported binder", "Kd: 7 nM", "P00001", "ST1"]))
                for key in ("SBP1", "SBP2")
            ],
        )
        (d / "Sequence_Information_of_SBPs.txt").write_text(
            f">SBP1|peptide binder\n{OTHER}\n>SBP2|Fab\n{SEQ}\n"
        )
        (row,) = synbip(self.cache, self.proteins)["P00001"]
        self.assertEqual(
            (row["partner_sequences"], row["affinity"], row["kd_nM"], row["status"]),
            (OTHER, "Kd: 7 nM", 7, "parent_only"),
        )

    def test_public_cli_works_without_antibody_indexes_and_reports_offline_misses(self):
        q = self.root / "q.fa"
        q.write_text(f">sample/one\n{SEQ}\n>sample/two\n{OTHER}\n")
        out = self.root / "out"
        cmd = [
            sys.executable,
            str(ROOT / "run_pipeline.py"),
            "--sources",
            "partners",
            "--partner-offline",
            "--partner-cache",
            str(self.root),
            "--partner-human-fasta",
            str(self.fasta),
            "--reference-dir",
            str(self.root / "empty"),
            "--query",
            str(q),
            "--outdir",
            str(out),
        ]
        subprocess.run(cmd + ["--partner-sources", "cellphonedb"], check=True, capture_output=True)
        self.assertFalse((out / "binding_candidates.csv").exists())
        rows = read_csv(out / "partner_evidence.csv")
        self.assertEqual({r["query_id"] for r in rows}, {"sample/one", "sample/two"})
        self.assertEqual({r["status"] for r in rows}, {"parent_only"})
        self.assertTrue((out / "sample_one/partner_candidates.csv").is_file())
        summary = read_csv(out / "partner_summary.csv")
        self.assertEqual({r["query_supported"] for r in summary}, {"False"})
        self.assertEqual({r["query_id"] for r in summary}, {"sample/one", "sample/two"})
        failed = subprocess.run(cmd + ["--partner-sources", "intact"], capture_output=True)
        self.assertNotEqual(failed.returncode, 0)
        self.assertEqual({r["status"] for r in read_csv(out / "partner_evidence.csv")}, {"source_error"})


@unittest.skipUnless(importlib.util.find_spec("biotite"), "Partner structure tests use Biotite")
class PartnerGeometry(unittest.TestCase):
    def test_covered_site_and_substituted_contact(self):
        from partner_structures import Chain, site_result

        p = Protein("P1", "TEST", "Test", SEQ)
        chain = Chain("P1", SEQ, "Test", SEQ, [(1, len(SEQ))], {n: n for n in range(1, len(SEQ) + 1)})
        q = SEQ[4:24]
        m = Match("fragment", {n + 4: n for n in range(len(q))}, 1)

        def score(sites, query=q, limit=0):
            return site_result({"A": sites}, {"A": chain}, query, p, m, limit, 0.8)

        self.assertEqual(score([5, 9])["status"], "site_conserved")
        self.assertEqual(score([3, 9])["site_covered"], 1)
        mutated = q[:4] + "Y" + q[5:]
        self.assertEqual(score([5, 9], mutated)["reason"], "site_difference_limit")
        self.assertEqual(score([5, 9, 11, 17, 20], mutated, 1)["status"], "site_similar")
        chain.modified_residues[9] = "PTR"
        self.assertEqual(score([5, 9])["reason"], "modified_target_contact")

    def test_real_biological_assemblies(self):
        from partner_structures import load_entry, structure_rows

        path = ROOT / "reference/partners"
        online = os.environ.get("AB_SEARCH_LIVE_PARTNERS") == "1"
        pdbs = ["2b0u", "3w11", "4asd"]
        if not online and (
            not (path / "human-reviewed-isoforms.fasta.gz").exists()
            or not all((path / "pdbe" / (p + ".cif.gz")).exists() for p in pdbs)
        ):
            self.skipTest("Set AB_SEARCH_LIVE_PARTNERS=1 to fetch the official PDBe test entries")
        cache = Cache(path, offline=not online)
        from partner_sources import UNIPROT

        proteins = cache.read("human-reviewed-isoforms.fasta.gz", UNIPROT, read_proteome)
        outputs = []
        for r, pdb_id in zip(read_query_records(ROOT / "examples/protein_partners.fasta"), pdbs):
            p, m, _, _ = resolve_query(r.sequence, proteins)
            outputs.append(
                structure_rows(load_entry(cache, pdb_id), pdb_id, proteins, p, r.sequence, m, 4.5, 0, 0.90)
            )
        fst, insr, kdr = outputs
        (activin,) = ranked(fst)
        self.assertEqual(
            (activin["partner_id"], activin["status"], activin["partner_stoichiometry"]),
            ("P08476", "exact_structure", "2"),
        )
        self.assertGreater(activin["site_covered"], 0)
        self.assertEqual(activin["site_size"], activin["site_covered"])
        (insulin,) = [r for r in insr if r.get("partner_id") == "P01308"]
        self.assertEqual(insulin["reason"], "additional_target_chain_outside_query")
        self.assertEqual(sorted(map(len, insulin["partner_sequences"].split("|"))), [21, 30])
        self.assertIn("731-746", insulin["complex_context"])
        self.assertEqual(ranked(insr), [])
        self.assertEqual([r["reason"] for r in kdr], ["target_domain_outside_query"])
        self.assertTrue(all(not r.get("partner_id") for r in kdr))


if __name__ == "__main__":
    unittest.main()
