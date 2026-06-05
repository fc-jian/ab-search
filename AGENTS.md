# AGENTS.md

## Project Purpose

This repository builds a SAbDab-backed antigen search and structure scoring pipeline.
Given a protein query FASTA, the pipeline:

1. Builds a BLAST database from SAbDab protein antigen chains.
2. Searches the query against those SAbDab chains with `blastp`.
3. Downloads matching antibody-antigen complex structures as cached `*.cif.gz` files.
4. Scores each candidate structure by antibody CDR-contacting antigen epitope residues.
5. Annotates antibody chains with ANARCI Chothia numbering by default and reports CDR-H/L regions.

## Environment

Expected tools and Python packages:

```bash
mamba install -c conda-forge -c bioconda biotite blast pandas biopython anarci
```

The current project environment already has `blastp`, `makeblastdb`, `pandas`,
`Bio`, `biotite`, and `anarci` available.

Use the local project skills when editing antibody or structure logic:

- `.agents/skills/anarci/SKILL.md` for ANARCI numbering and CDR extraction.
- `.agents/skills/biotite/SKILL.md` for mmCIF loading, chain selection, and distance calculations.

## Reference Inputs

Required files under `reference/`:

- `sabdab_summary_all.tsv`: SAbDab summary table.
- `pdb_seqres.fasta.gz`: RCSB PDB derived sequence FASTA.

Generated reference files:

- `reference/sabdab_seq.fasta`: SAbDab protein antigen-chain FASTA.
- `reference/blastdb.*`: BLAST protein database files.
- `reference/pdb_files/*.cif.gz`: downloaded structure cache.

## Main Command

Run the full test pipeline:

```bash
python run_pipeline.py \
  --query test_input.fasta \
  --outdir results/test_input \
  --threads 4
```

Run the multi-record test:

```bash
python run_pipeline.py \
  --query test_multi_input.fasta \
  --outdir results/test_multi_input \
  --threads 4
```

Important outputs:

- `results/test_input/blastp_hits.tsv`: BLAST hits.
- `results/test_input/download_manifest.tsv`: structure download/cache manifest.
- `results/test_input/structure_scores.csv`: final epitope and antibody annotation table.
- `results/test_input/merged_antibodies.csv`: exact heavy/light antibody sequence groups,
  including `support_pdb_ids` separated by `|` and `epitope_identity_range`.
- `results/test_input/structure_score_errors.tsv`: per-hit scoring failures.
- `results/test_input/logs/`: command logs.

For multi-record FASTA input, `run_pipeline.py` runs records sequentially and
creates one subdirectory per FASTA record ID. The output root also contains:

- `query_manifest.tsv`: record-to-subdirectory mapping.
- `merged_antibodies.csv`: antibody groups recomputed across every record's
  `structure_scores.csv`.

## Manual Pipeline

Prepare antigen-chain FASTA:

```bash
python filter_seqs.py \
  --fasta reference/pdb_seqres.fasta.gz \
  --sabdab reference/sabdab_summary_all.tsv \
  --output reference/sabdab_seq.fasta
```

Build BLAST database:

```bash
makeblastdb -in reference/sabdab_seq.fasta -dbtype prot -out reference/blastdb
```

Run BLAST:

```bash
mkdir -p results/test_input
blastp -query test_input.fasta -db reference/blastdb \
  -out results/test_input/blastp_hits.tsv \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs" \
  -max_target_seqs 10000 \
  -evalue 1e-5 \
  -num_threads 4
```

Download structures:

```bash
python download_pdb.py \
  --hits results/test_input/blastp_hits.tsv \
  --pdb-cache reference/pdb_files \
  --min-identity 0.7 \
  --manifest results/test_input/download_manifest.tsv
```

Score structures:

```bash
python score_structures.py \
  --query-fasta test_input.fasta \
  --hits results/test_input/blastp_hits.tsv \
  --sabdab reference/sabdab_summary_all.tsv \
  --pdb-seqres reference/pdb_seqres.fasta.gz \
  --pdb-cache reference/pdb_files \
  --output results/test_input/structure_scores.csv \
  --merged-antibody-output results/test_input/merged_antibodies.csv \
  --errors-output results/test_input/structure_score_errors.tsv \
  --min-identity 0.7 \
  --scheme chothia
```

## Script Notes

- `filter_seqs.py` defaults to `--scope antigen`, which keeps only protein antigen
  chains listed in SAbDab. Use `--scope antibody-antigen` to also include listed
  H/L antibody chains, or `--scope pdb` to include all protein chains from SAbDab
  PDB entries.
- `download_pdb.py` accepts `--min-identity 0.7` as 70 percent identity. Values
  above `1.0` are treated as already being percentages.
- `score_structures.py` uses Biotite with `model=1`, author chain/residue IDs,
  highest-occupancy alternate locations, CA distance threshold `8.0`, and minimum
  heavy-atom distance threshold `4.5`. It reads the first biological assembly via
  `pdbx.get_assembly()` by default; if no matching CDR-mediated ab-ag complex is
  found there, it falls back to the asymmetric unit and records the fallback in
  `structure_warnings`. When `--pdb-seqres` is supplied, antibody ANARCI annotation
  uses full SEQRES chain sequences while epitope geometry still uses observed
  structure residues.
- CDRs are Chothia by default (`--scheme chothia`). Epitope contacts are computed
  only against antibody CDR atoms.
- Multichain binding penalty: antibody CDR-contacting protein residues from non-target
  chains are included in the epitope denominator and treated as missing/mutated,
  unless the chain is a likely homodimer copy of the target antigen chain. The
  default homodimer exclusion threshold is `--homodimer-identity-threshold 0.95`.
- `merged_antibodies.csv` groups exact antibody sequences by
  `seq_heavy_full + seq_light_full` and reports antibody annotation columns plus
  `support_pdb_ids` and `epitope_identity_range`.
- ANARCI uses Chothia by default. Failed antibody annotations are reported in the
  `anarci_warnings` column rather than aborting the whole scoring run.

## Editing Guidelines

- Keep scripts standalone and runnable from the repository root.
- Preserve single-record output behavior. For multi-record FASTA input, keep each
  record's full pipeline output isolated in its own sanitized record-name subdirectory
  and regenerate the root-level merged antibody summary from all record score tables.
- Prefer cached reference data and structure files; do not redownload files that
  already exist and are non-empty.
- Keep output tables stable and append extra columns only after the documented
  core columns.
- Preserve chain IDs exactly as provided by SAbDab and PDB/mmCIF author fields.
- When changing structure scoring, test with `./test_input.fasta` and inspect
  `./results/test_input/structure_scores.csv`.
