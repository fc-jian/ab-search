---
name: sabdab-reference
description: Use when preparing a SAbDab reference snapshot, diagnosing missing chains or BLAST hits, or assessing SAbDab2 import compatibility in ab-search.
---

# SAbDab reference compatibility

## Reproduce an existing index

The current reader in `filter_seqs.py` expects a tab-delimited legacy summary with `pdb`, `Hchain`, `Lchain`, `antigen_chain` and `antigen_type`. Use PDB SEQRES headers compatible with that snapshot and preserve chain case. The default `--scope antigen` indexes protein antigen chains. `antibody-antigen` adds antibody chains and `pdb` selects all protein chains from the chosen entries.

`prepare_reference()` in `run_pipeline.py` records the summary and SEQRES digests and scope for the generated FASTA. It checks the FASTA digest when deciding to rebuild BLAST. Use a separate reference directory when comparing snapshots. Verify that expected PDB/chain pairs reach the FASTA and that ordinary queries resolve to the intended cached structures. A successful BLAST database build alone cannot detect incorrect joins.

## Assess a SAbDab2 source

SAbDab2 exports an instance-level CSV with fields including `INSTANCE`, `PDB_ID` and `SABDAB_ID`. It supplies IMGT-numbered structures and can split several domains on one deposited chain into separate logical chains. The legacy reader has no SAbDab2 adapter. Several current ID parsers also assume four-character PDB accessions.

For a proposed import, identify the accession and instance, distinguish SAbDab2 logical chains from deposited author chains, and map each domain or antigen segment through the retained PDB label chain and residue IDs. Record unresolved or ambiguous mappings. Use matched chain and type entries for mixed-component antigens. Compare a split-chain case, an extended accession and an assembly with repeated chains. When validating a merge, use the full heavy/light sequences for `merged_antibodies.csv` and the variable-domain numbering for `SABDAB_ID`.

The reference format and numbering definitions can change independently. Use the [SAbDab2 field descriptions](https://sabdab.opig.stats.ox.ac.uk/about), [upstream API schema](https://sabdab.opig.stats.ox.ac.uk/api/docs) and [RCSB identifier conventions](https://www.rcsb.org/docs/general-help/identifiers-in-pdb).
