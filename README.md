# ab-search

Search a protein FASTA for published antibody sequences and, optionally, human protein or peptide binding partners. Results link each candidate to its source record and show whether the evidence applies to the submitted sequence, a mapped region, or the parent protein.

## Install

Use Python 3.10 or newer. The patent and assay indexes use NumPy, pandas and Biopython. For all backends, create a Conda environment on a platform with Bioconda packages, such as Linux or WSL.

```bash
mamba create -n ab-search -c conda-forge -c bioconda python=3.11 numpy pandas biopython scipy biotite blast hmmer anarci
conda activate ab-search
```

SAbDab scoring uses `blastp`, `makeblastdb` and `hmmscan` from this environment. Run commands from the repository root.

## Run a search

Supply your own `query.fasta` containing one or more protein sequences with FASTA headers. To search published patent records, build the local index once, then run the query.

```bash
python prepare_patents.py --download
python run_pipeline.py --query query.fasta --sources patents --outdir results/query
```

Omitting `--sources` selects the prepared antibody sources. Select `partners` explicitly, or add `--protein-partners` to an antibody search. Each multi-record query produces one directory per FASTA record, a `query_manifest.tsv` and combined tables at the output root.

## Choose a source

| Source | Scope and preparation |
| --- | --- |
| `patents` | Indexes POISK records and optional PLAbDab paired-chain references. For local exports, use `python prepare_patents.py --poisk poisk.csv --plabdab paired_sequences.csv.gz`. Supply `--curations regions.tsv` for independently verified patent regions. [POISK release](https://zenodo.org/records/21775623) · [PLAbDab](https://opig.stats.ox.ac.uk/webapps/plabdab/) |
| `avida` | Looks up measured VHHs for exact human IL-6 assay constructs. Build with `python prepare_hil6.py --download`. The upstream dataset is distributed under CC BY-NC 4.0. [AVIDa-hIL6 dataset and terms](https://huggingface.co/datasets/COGNANO/AVIDa-hIL6) |
| `sabdab` | Requires a compatible legacy SAbDab summary TSV and PDB SEQRES FASTA. The first run builds a protein BLAST index and caches matching experimental structures. [SAbDab2 format](https://sabdab.opig.stats.ox.ac.uk/about) · [RCSB SEQRES FASTA](https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz) |
| `partners` | Resolves the input to a reviewed human UniProt protein and looks up PDBe/SIFTS, CellPhoneDB, Guide to Pharmacology, IntAct and SYNBIP evidence. Restrict sources with `--partner-sources` and reuse a prepared snapshot with `--partner-cache DIR --partner-offline`. [PDBe/SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html) · [Guide to Pharmacology downloads](https://www.guidetopharmacology.org/download.jsp) |

Reference files and indexes are prepared on your machine under `reference/` by default. Each provider supplies its own data and usage terms.

### SAbDab structure scoring

Place your compatible summary at `reference/sabdab_summary_all.tsv` and the SEQRES file at `reference/pdb_seqres.fasta.gz`, or pass their paths explicitly.

```bash
python run_pipeline.py --query query.fasta --sources sabdab --sabdab-summary summary.tsv --pdb-seqres pdb_seqres.fasta.gz --threads 4
```

The reader expects legacy columns such as `pdb`, `Hchain`, `Lchain` and `antigen_chain` matching author chains in deposited structures. Current SAbDab2 summaries use an instance-level CSV, extended identifiers and logical chain IDs that can differ from deposited chain IDs. They require a separate mapping step before this reader can use them. [SAbDab2 format and identifiers](https://sabdab.opig.stats.ox.ac.uk/about)

The scorer numbers antibody domains from SEQRES when available, maps the CDRs to observed residues and uses the first biological assembly. If that assembly is unavailable or yields no matching CDR contacts, it tries the asymmetric unit and records the reason in `structure_warnings`. A contact meets either the Cα threshold (8 Å) or the minimum heavy-atom threshold (4.5 Å). The scorer defines Chothia and IMGT CDR ranges. Chothia is the default.

### Human protein and peptide partners

```bash
python run_pipeline.py --query query.fasta --sources partners --outdir results/partners
```

PDBe provides experimental assemblies and SIFTS residue mappings. Other sources can establish a relationship to the parent protein. Read the construct, partner form, original measurement and review reason to determine which conclusions apply to a fragment. For Guide to Pharmacology, use `GTP_API_KEY` when required by the provider or supply the matching official tables in the local cache.

## Read the results

| File | What to look for |
| --- | --- |
| `binding_candidates.csv` | One row per query and exact antibody chain pair. `evidence_tier` identifies the best supporting record. `evidence_grade` is an ordinal rank. |
| `binding_evidence.csv` | Source-level decisions, coverage, substitutions, reported affinity, conditions and review reasons. A probability is available only with applicable calibration. |
| `plabdab_references.csv` | `exact_chain_pair` identifies a matching pair. `target_mention_only` records a document to investigate. |
| `structure_scores.csv` and `merged_antibodies.csv` | Per-hit geometry and exact full-chain groups, with supporting PDB IDs and the observed epitope identity range. Download and scoring failures have separate manifests. |
| `protein_binders.csv` | One fragment-supported binder per query and partner ID. `supported_forms_json` retains the deposited sequences, stoichiometry, modifications and supporting records for each form. |
| `partner_candidates.csv`, `partner_summary.csv` and `partner_evidence.csv` | Supported partner forms, relationships at the parent-protein level, and the source records behind each decision. |
| `run_metadata.json` | Active sources, input and reference fingerprints, search parameters and source availability. Check this file when interpreting an empty result. |

Run `python run_pipeline.py --help` for search options. For a list containing one protein sequence per line, `python find_chain_name.py --help` describes how to look up the best chain description in an existing protein BLAST database.
