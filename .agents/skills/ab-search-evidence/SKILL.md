---
name: ab-search-evidence
description: Use when explaining an ab-search candidate, investigating missing or conflicting hits, or integrating an evidence source into the antibody or human-partner outputs.
---

# Interpret a binding record

Start with the query record ID and the source information in `run_metadata.json`. For a multi-record FASTA, `query_manifest.tsv` identifies the per-query directory. An unavailable source or missing offline snapshot can leave an empty candidate table. Check the per-source error or review reason before interpreting that result.

| Evidence | How to trace it |
| --- | --- |
| Antibody structure | Follow `best_evidence_record` or `evidence_records` from `binding_candidates.csv` to `binding_evidence.csv` and the matching structure score. Inspect observed CDR contacts, mapping coverage, multichain context and assembly warnings. |
| Patent or literature | Check the POISK record ID, patent, reported chains and construct. PLAbDab `exact_chain_pair` is a paired-sequence reference. `target_mention_only` is a document-level lead. Conflicting reports remain visible in the evidence table. |
| Assay | Read the tested sequence, construct and assay context. The AVIDa-hIL6 label applies to a measured construct. Record raw affinity, qualifier and units together. A probability needs applicable calibration. |
| Human protein or peptide partner | Use `partner_summary.csv` for parent-level relationships, `partner_candidates.csv` for supported forms and `protein_binders.csv` for the grouped view. Trace the form key and assembly in `partner_evidence.csv`. A relationship to a precursor alone leaves binding to the submitted fragment unresolved. |

`binding_search.py` ranks exact input chain pairs. Its `evidence_grade` is an ordinal rank. Fv QC records chain provenance and warnings separately. Keep a recovered Fv linked to the original record and pair key. `protein_partners.py` groups supported forms by partner, deposited sequence, stoichiometry and modifications. Measurements retain their original record and construct.

When adding a source, exercise a supported record, an ambiguous or negative record, and an unavailable snapshot. Check the normal CLI outputs for one and multiple queries. [PLAbDab pairing methods](https://doi.org/10.1093/nar/gkad1056) · [PDBe/SIFTS mappings](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)
