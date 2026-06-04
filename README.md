# ab-search

SAbDab-backed antigen search and antibody-antigen epitope scoring pipeline.

Given a query protein FASTA, the pipeline searches SAbDab antigen chains with
BLAST, downloads matching antibody-antigen complex structures, and reports
structure-based epitope identity plus antibody Chothia CDR annotations.

## 1. Environment

Install the required command-line tools and Python packages:

```bash
mamba install -c conda-forge -c bioconda biotite blast pandas biopython anarci
```

Check the environment:

```bash
blastp -version
makeblastdb -version
python - <<'PY'
import pandas
import Bio
import biotite
import anarci
print("environment OK")
PY
```

## 2. Prepare Reference Data

Create the reference directory:

```bash
mkdir -p reference
```

Download the SAbDab summary table from:

```text
https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/search/?all=true#downloads
```

Save it as:

```text
reference/sabdab_summary_all.tsv
```

Download all PDB SEQRES sequences:

```bash
wget -O reference/pdb_seqres.fasta.gz \
  https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
```

Build the SAbDab antigen-chain FASTA:

```bash
python filter_seqs.py \
  --fasta reference/pdb_seqres.fasta.gz \
  --sabdab reference/sabdab_summary_all.tsv \
  --output reference/sabdab_seq.fasta
```

Build the BLAST database:

```bash
makeblastdb \
  -in reference/sabdab_seq.fasta \
  -dbtype prot \
  -out reference/blastdb
```

## 3. Run The Pipeline

Use the all-in-one runner:

```bash
python run_pipeline.py \
  --query test_input.fasta \
  --outdir results/test_input \
  --threads 4
```

For a different query:

```bash
python run_pipeline.py \
  --query path/to/query.fasta \
  --outdir results/query_name \
  --threads 4
```

The runner reuses existing reference files, BLAST DB files, and cached mmCIF
structures when they are already present.

## 4. Main Outputs

For the example command, outputs are written under `results/test_input/`:

- `blastp_hits.tsv`: BLAST tabular hits.
- `download_manifest.tsv`: downloaded or cached PDB mmCIF files.
- `structure_scores.csv`: per-hit structure scoring table.
- `merged_antibodies.csv`: exact heavy/light antibody sequence groups with
  `support_pdb_ids` separated by `|` and `epitope_identity_range` summarizing the
  per-structure identity range for that antibody group.
- `structure_score_errors.tsv`: per-hit scoring failures, if any.
- `logs/`: command logs.

## 5. Scoring Notes

`structure_scores.csv` reports epitope identity using antibody CDR-contacting
antigen residues:

- CDRs are defined with Chothia numbering by default (`--scheme chothia`).
- Only antigen residues contacting antibody CDR atoms are counted as epitope
  residues.
- Target antigen-chain epitope residues are aligned to the query sequence and
  scored as identical or mutated.
- If antibody CDRs also contact other non-antibody protein chains, those contact
  residues are included in the epitope denominator as missing multichain binding
  residues.
- Other chains are excluded from this multichain penalty when their sequence
  identity to the target antigen chain is at least `--homodimer-identity-threshold`
  (default `0.95`), treating them as target-chain homodimer copies.
- Structures are read from the first biological assembly by default. If no
  matching CDR-mediated antibody-antigen complex is found there, scoring falls
  back to the asymmetric unit and reports this in `structure_warnings`.

Useful scoring parameters:

```bash
python run_pipeline.py \
  --query test_input.fasta \
  --outdir results/test_input \
  --ca-distance-threshold 8.0 \
  --atom-distance-threshold 4.5 \
  --homodimer-identity-threshold 0.95 \
  --scheme chothia
```

## 6. Manual Structure Scoring

After BLAST and structure download, scoring can be run directly:

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
  --scheme chothia
```
