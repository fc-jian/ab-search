---
name: anarci
description: Use for antibody domain numbering, scFv chain assignment, CDR extraction and IMGT Fv quality checks in ab-search.
---

# Antibody numbering in ab-search

Choose the entry point by the output being produced.

- **Structure scores.** Use `run_anarci_domains()` and the CDR definitions in `score_structures.py`. The default is Chothia. IMGT is also defined. Other schemes accepted by ANARCI can fall back to the scorer's Chothia intervals. Treat them as unsupported for interpreting structure scores.
- **Fv quality.** Use `FvChecker` in `fv_qc.py` for the `fv_*` columns. This path uses IMGT independently of the structure-scoring option. Raw full-chain sequences, numbered domains and a source-verified Fv have different meanings and grouping keys.

## Preserve the mapping

Pair each numbered domain with the metadata at the same sequence and domain index. A single scFv can contain a light domain followed by a heavy domain. For the raw ANARCI API, input sequence offsets are zero-based and the end is inclusive. Keep the insertion code with the scheme position and preserve the reported residue order. Gap characters consume no input residue.

Use the existing FASTA reader for the relevant command. An internal stop in a query must surface as an input error. Keep annotation failures in the relevant warning column and retain the source record.

Use full-chain SEQRES for numbering when available. Map sequence indices to the observed chain before selecting CDR atoms. The numbering engine, numbering scheme and CDR definition are separate choices. A backend change needs a comparison of domain types, spans, region assignments and QC warnings, including missing data and multi-domain inputs. Fv QC also uses ANARCI's germline tables.

See the [ANARCI API example](https://github.com/oxpig/ANARCI/blob/master/Example_scripts_and_sequences/anarci_API_example.py) and [AbNumber region definitions](https://abnumber.readthedocs.io/en/stable/#chain).
