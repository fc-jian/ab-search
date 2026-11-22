# ab-search

Search a protein FASTA for published antibody sequences. The pipeline connects
SAbDab structures, POISK patent records, Oxford PLAbDab paired sequences and an
experimentally labelled human IL-6 panel. Each candidate links back to the
records that support it.

## Quick start

Python 3.10+ with NumPy, pandas and Biopython runs the patent and human IL-6
searches. These two indexes can be built and searched on Windows or Linux.

```bash
python prepare_patents.py --download
python prepare_hil6.py --download
python run_pipeline.py --sources patents avida --query examples/human_controls.fasta --outdir results/human_controls
python run_pipeline.py --sources avida --query examples/avida_hil6_wt.fasta --outdir results/il6_wt
```

To include structures, prepare the SAbDab references and tools described below.
Running without `--sources` searches all prepared indexes. Use
`--sources sabdab patents avida` to explicitly require all three backends.
The main output is `binding_candidates.csv`. A multi-query file produces a
combined candidate table with a separate row per query and antibody sequence pair.

## Optional protein and peptide partners

Add `--protein-partners` to an antibody search to also find published protein
and peptide binders of reviewed human proteins. Use `--sources partners` for
a standalone partner search. The input can be a full protein, mature peptide
or fragment. The query sequence is matched locally to the reviewed human
UniProt proteome including isoforms. Public services receive accessions and
PDB IDs. The structural backend uses Biotite and SciPy.

```bash
python run_pipeline.py --query examples/protein_partners.fasta --protein-partners --outdir results/combined
python run_pipeline.py --sources partners --query examples/protein_partners.fasta --outdir results/partners
python run_pipeline.py --sources partners --query examples/protein_partners.fasta --partner-pdb 2b0u 3w11 4asd --outdir results/partner_controls
```

| Source name for `--partner-sources` | Data used |
| --- | --- |
| `pdbe` | PDBe's best-structures API discovers entries across the experimental PDB. SIFTS in the updated mmCIF maps the deposited polypeptide residues in each biological assembly to UniProt. [PDBe API](https://www.ebi.ac.uk/pdbe/pdbe-rest-api) · [SIFTS mappings](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html) |
| `cellphonedb` | Version 5.0.0 protein–protein interaction records, searched in both directions with receptor/ligand roles and complex subunits. [Public release](https://github.com/ventolab/cellphonedb-data/tree/v5.0.0/data) |
| `gtopdb` | Endogenous pairings and human peptide interactions joined to the matching peptide release. Records can provide mature chain sequences, modifications, assay descriptions, original units and references. Set the environment variable `GTP_API_KEY` for authenticated downloads, or supply the three official matching CSVs in the cache's `gtopdb/` directory. [Downloads](https://www.guidetopharmacology.org/download.jsp) · [Programmatic access](https://blog.guidetopharmacology.org/2026/11/20/gtopdb-download-changes-coming-soon/) |
| `intact` | Paginated PSICQUIC MITAB 2.7 records for human participants. Rows retain the experiment, participant features, negative flag decision, and publication. Associations and expanded complexes remain identifiable. [IntAct](https://www.ebi.ac.uk/intact/) |
| `synbip` | Published target–binder records and original sequences for human targets and supported non-Ig scaffolds including monobodies, DARPins, affibodies and anticalins. Raw affinity and the source record accompany each hit. [SYNBIP downloads](https://synbip.idrblab.net/full-data-download) |

The default partner search uses all five sources. The source snapshots live
under `reference/partners/`. Use `--partner-cache DIR` for a separate snapshot
and `--partner-offline` to reproduce a run from that snapshot. Use
`--partner-human-fasta FILE` for a local reviewed human UniProt FASTA.
`--partner-pdb` restricts the discovered structures by ID. An optional
`--partner-max-structures N` limits the number per query, prioritizing
structures that span the input region. The run metadata records the number
discovered, selected IDs, source outcomes, URLs and file hashes. A source
failure is visible in the evidence and metadata. An all-source failure for a
resolved query returns a nonzero exit code after writing the audit.

| Partner output | Contents |
| --- | --- |
| `partner_summary.csv` | One row per query and partner ID including all sources and the number of distinct forms supported for the input. Start here to see both established parent-protein relationships and fragment-supported hits. |
| `partner_candidates.csv` | Structure-supported forms with deposited partner chain sequences, stoichiometry, supporting structures and contact coverage. |
| `partner_evidence.csv` | Source-level records, assembly chains, query-relative contact positions, construct context, assay fields and decisions. Reported point values for Kd, Ki and IC50 have separate nM columns with their original qualifier and raw affinity. Sequence basis identifies reference precursors and reported partner chains. |

Multi-record runs produce these tables per query and at the output root.
GtoPdb, CellPhoneDB, IntAct and SYNBIP can establish a reported relationship
even when its assayed construct is unspecified. Such records receive
`parent_only`. A structure requires actual heavy-atom contact (default 4.5 Å)
with a human polypeptide in a biological assembly. An exact deposited
construct can qualify directly. For transfer, every observed target contact
must map uniquely into the query and satisfy the site's substitution and
identity limits (default zero substitutions). Extra target segments, a
multi-protein interface or a modified target contact are recorded as context.
The peptide label in structural output means that each partner chain is at
most 50 residues. Chemical components and nucleic acids are excluded from
the partner search. Evidence grades describe published experiments and structures.

The three-structure command above exercises these controls. The FST query
matches the 288-residue follistatin chains bound to the activin A dimer
(2B0U). The INSR query spans residues 28–337 and insulin also contacts a
separate receptor αCT peptide (3W11). The mapped KDR kinase-domain construct
(4ASD) falls outside the KDR 141–320 input. To include the real structures in
the test suite on a fresh checkout, run
`AB_SEARCH_LIVE_PARTNERS=1 python -m unittest discover -s tests -v` in the
documented environment. Cached controls run offline on subsequent invocations.
[2B0U](https://www.rcsb.org/structure/2B0U) · [3W11](https://www.rcsb.org/structure/3W11) · [4ASD](https://www.rcsb.org/structure/4ASD)

### Sources and traceability

| Source | Searchable evidence | Reference preparation |
| --- | --- | --- |
| SAbDab + RCSB | Experimental antibody–antigen complexes and CDR contacts | Sections 1 and 2 below |
| POISK 1.1 | Antibody chains, target sequences, patent and record IDs, assay text, coordinate evidence and audit fields | `prepare_patents.py` reads the official CSV or a local CSV.GZ. [Versioned release](https://zenodo.org/records/21775623) |
| PLAbDab | Exact paired-chain references and potential leads from target names mentioned in a document | The same builder indexes the Oxford paired CSV.GZ. [Source and pairing methods](https://doi.org/10.1093/nar/gkad1056) |
| AVIDa-hIL6 | Assay labels for 38,599 VHH sequences tested across 31 human IL-6 constructs | `prepare_hil6.py` verifies both files from the pinned release. The data are distributed under CC BY-NC 4.0. [Dataset](https://huggingface.co/datasets/COGNANO/AVIDa-hIL6) · [Study](https://arxiv.org/abs/2306.03329) |

To index local patent exports, run
`python prepare_patents.py --poisk path/to/poisk.csv --plabdab path/to/paired_sequences.csv.gz`.
The builder stages a new SQLite file and replaces the old index after successful
import. POISK 1.1 contains 93,193 records. In the verified release, 81,795 have
a searchable antigen and at least one antibody chain. The Oxford snapshot used
here has 176,894 paired records. Source hashes and import counts are recorded
in the index and copied to each search's `run_metadata.json`.

### Fragment and homolog evidence

1. An exact antigen query can retrieve a reported patent pair or experimental complex.
   Admission of a POISK pair requires a confirmed audit without flags, a supported
   protein construct, and complete, unambiguous paired chains or an explicit VHH.
   Incomplete chains, contradictory results, complex-dependent binding and
   canonical-only target assignments remain visible in `binding_evidence.csv`.
2. For a fragment or homolog, every residue of the verified region or structural
   contact set must be covered by one unambiguous alignment. Individual BLAST HSPs
   are never combined to claim coverage. The search checks literal containment
   first. For an equal-length pair of at least eight residues with ≥90% identity,
   a unique, ungapped optimum can retain terminal substitutions. Local
   alignments need at least 12 aligned residue pairs and the configured overall
   identity (default 70%, including alignment gaps in the denominator).
3. Missing residues, unknown site letters and insertions within a mapped linear
   region prevent transfer. The default allows zero site substitutions. For a
   controlled relaxation, `--max-site-substitutions 1 --min-site-identity 0.90`
   applies both thresholds. A change outside the site can leave the candidate eligible.
   Conserved sequence alone leaves folding, accessibility and assay conditions uncertain.
4. Cross-construct ranking currently covers source-identified human, mouse, rat,
   chicken and bovine single-protein targets. Exact published target lookups remain
   available for other organisms. In the structural path an unassignable contact
   set or additional contacting protein chain prompts review.
5. POISK structural assertions require residue letters matching the patent-listed
   sequence and structure-determination evidence. A bare numeric range, domain
   name or competition result cannot localize a fragment by itself. The curated
   table anchors selected regions by unique peptide sequence and release-target
   SHA-256. It also records construct and terminal requirements. For example, the
   Aβ experiment refers to a processed N terminus, while the archived target
   sequence is its APP precursor. [Curated records](data/verified_patent_sites.tsv) · [Original Aβ experiment](https://patents.google.com/patent/US7189819B2/en)

The human PSCA controls exercise a real coordinate offset. The original mapping
used a 123-residue PSCA sequence. The 114-residue reference in POISK record
`008508` starts nine residues later. Mouse 1G8's mapped 40-residue region is
located by its sequence. `PSCA_middle` recovers mu1G8, while `PSCA_C_terminal`
retains only one residue of the region and produces an empty candidate table.
The patent localizes a binding region by fragment experiments. Individual
atomic contacts within that region remain unresolved. [Fragment mapping](https://patents.google.com/patent/US6790939B2/en) · [1G8 source patent](https://patents.google.com/patent/US8088908B2/en)

### Reading the scores

| Output | Interpretation |
| --- | --- |
| `binding_candidates.csv` | One exact chain pair per query with its source records. `evidence_grade` is ordinal, from 4 for an exact structure or positive assay, through 3 for an exact patent report, 2 for a conserved mapped site and 1 for a mapped site with permitted substitutions. |
| `binding_evidence.csv` | Per-record decisions, coverage, substitutions, raw affinity and assay context. A review reason records why sequence evidence could not establish a candidate. |
| `plabdab_references.csv` | `exact_chain_pair` links the published pair to a candidate. `target_mention_only` records additional pairs found by a target-name or gene-name token from an exact POISK target or SAbDab complex. Both retain the Oxford ID, pairing method and source title. Document mentions do not increase the evidence grade. |
| `run_metadata.json` | Sources, fingerprints, thresholds and the result paths for this run. SAbDab-specific BLAST limits and per-query structure error counts are included when active. |

`binding_probability` estimates the chance of a positive result for the
specified assay endpoint. `probability_status` identifies the applicable model
or explains a blank estimate. A missing probability is unknown.

**Human IL-6 panel.** The bundled procedure estimates the AVIDa binary
phage-display label for previously profiled VHHs and exact published constructs
(including their tags). Five folds split entire antigen constructs by a stable
hash. For each test fold, three other folds supply each VHH's response history,
and one separate fold calibrates positivity by response-frequency and assay-count
bins. The evaluated VHH is excluded from the calibration counts. Missing labels
stay missing. An estimate requires three prior observations and 30 calibration
measurements. The estimate is `(positive + 1) / (measured + 2)`.

Validation includes positive and negative labels before positive candidates are
selected. Of 573,891 measured pairs, 553,287 could be scored, including 15,257 of
20,980 positives. The pooled held-out Brier score was 0.00206 against 0.02684 for
calibration-fold prevalence on the same evaluated rows (lower is better).
`reference/avida_validation.csv` gives counts and scores per held-out construct.
This evaluates transfer across constructs within one experimental panel. Generalization
to a new VHH repertoire, another target, or another assay requires separate validation.

**Independent assays for another query.** Pass `--calibration-assays labels.csv`
to estimate the positive fraction of an evidence tier for one endpoint. Columns
are `query_sha256,antibody_key,evidence_tier,assay,study_id,bound`. Use the same
query and chain keys as the candidate table, one protocol and positivity threshold
per file, and independently measured 0/1 outcomes. Assign tiers using retrieval
evidence available before the independent assay. Use the original patent or PDB
identifier as the study ID where applicable. The estimator removes the candidate
itself and overlapping support documents, collapses duplicate antibody labels and
requires both outcomes across at least ten distinct antibodies and two studies.
It uses the same add-one estimate and exposes its sample count.

**Equilibrium occupancy.** With `--free-antibody-nm C`, a confirmed exact patent
record with an explicitly labelled point Kd from a supported binding assay can
also report `equilibrium_bound_fraction = C / (C + Kd)` under a 1:1 equilibrium
model. Units are converted to nM. This field describes occupancy at the chosen
concentration. IC50, bounds and affinities reported for a different construct
leave it blank.

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

Multi-record FASTA files are processed sequentially:

```bash
python run_pipeline.py \
  --query test_multi_input.fasta \
  --outdir results/test_multi_input \
  --threads 4
```

For multi-record input, each FASTA record ID gets its own output subdirectory:

```text
results/test_multi_input/
|-- M1R/
|-- M1R_copy/
|-- merged_antibodies.csv
`-- query_manifest.tsv
```

Record IDs are preserved as query IDs. Characters unsuitable for directory names
are replaced with `_`; the runner stops if two record IDs would map to the same
directory name.

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

For multi-record input:

- `<record_id>/`: the complete normal output set for one FASTA record.
- `query_manifest.tsv`: maps input record IDs to generated subdirectories and files.
- Top-level `merged_antibodies.csv`: exact antibody sequence groups recomputed from
  all records. Its `query_id_list`, `support_pdb_ids`, support count, and
  `epitope_identity_range` summarize the complete multi-query run.

Antibody chain descriptions:

- `structure_scores.csv` includes `heavy_chain_description` and
  `light_chain_description`, read from the corresponding PDB SEQRES FASTA header.
- `merged_antibodies.csv` includes `heavy_chain_description_list` and
  `light_chain_description_list`, merged with the same unique `|`-separated logic
  used for `antigen_name_list`.

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

## 7. Find Chain Names

Use `find_chain_name.py` to annotate plain one-sequence-per-line protein input
against any BLAST protein database:

```bash
python find_chain_name.py \
  --input sequences.txt \
  --blastdb reference/blastdb \
  --output chain_names.csv \
  --min-identity 0.7 \
  --threads 4
```

If `blastp` is not on `PATH`, pass it explicitly:

```bash
python find_chain_name.py \
  --input sequences.txt \
  --blastdb reference/blastdb \
  --output chain_names.csv \
  --blastp /home/jianfc/miniforge3/envs/bioinfo/bin/blastp
```

The output CSV contains one row per input line, including the best matching
subject ID, FASTA description, identity, query coverage, e-value, and bitscore.
