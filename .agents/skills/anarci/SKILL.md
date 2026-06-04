---
name: anarci
description: Use ANARCI for antibody variable-domain sequence annotation, chain-type classification, numbering with IMGT/Kabat/Chothia/Martin/AHo schemes, and CDR/framework extraction from VH/VL sequences.
---

# ANARCI antibody annotation skill

Use this skill when writing code or shell commands for antibody variable-domain annotation with **ANARCI**.

## Scope

Use ANARCI for:

- Numbering VH, Vκ, Vλ, and TCR variable-domain amino-acid sequences.
- Classifying chain type.
- Extracting CDR and framework regions.
- Batch-processing FASTA files.
- Producing TSV/CSV annotation tables.

Do **not** run ANARCI directly on nucleotide sequences. Translate first.

## Defaults

Unless the user specifies otherwise:

- Use `scheme="imgt"`.
- Treat input as antibody variable-domain amino-acid sequences.
- Accept mixed heavy/light FASTA files.
- Preserve input sequence IDs.
- Write failed sequences to a separate output.
- Use one output row per recognized domain.
- Report distances/positions/residue ranges using the selected numbering scheme.

## Installation check

```bash
ANARCI --help

python - <<'PY'
from anarci import anarci
print("ANARCI import OK")
PY
````

Typical environment:

```bash
conda create -n anarci python=3.10 -y
conda activate anarci
conda install -c bioconda hmmer -y
pip install ANARCI
```

If installation fails, use the project’s documented GitHub installation route.

## Command-line usage

```bash
ANARCI -i input.fasta -s imgt -o output_imgt
```

Common schemes:

```text
imgt
kabat
chothia
martin
aho
```

Use IMGT by default for repertoire and antibody-engineering workflows.

## Python API

```python
from anarci import anarci

records = [
    ("seq1", "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSA..."),
    ("seq2", "DIQMTQSPSSLSASVGDRVTITCRASQ..."),
]

numbering, alignment_details, hit_tables = anarci(
    records,
    scheme="imgt",
    output=False,
)
```

ANARCI returns:

* `numbering`: numbered residue assignments.
* `alignment_details`: domain and chain metadata.
* `hit_tables`: HMM hit details.

Downstream code should parse these defensively, because local ANARCI versions/forks may differ slightly.

## IMGT regions

For IMGT antibody variable domains:

```python
IMGT_REGIONS = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}
```

Do not mix CDR definitions across schemes. If using Kabat, Chothia, Martin, or AHo, use that scheme’s own boundaries.

## Minimal FASTA parser

```python
def read_fasta(path: str) -> list[tuple[str, str]]:
    records = []
    name = None
    chunks = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper().replace("*", "")))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        records.append((name, "".join(chunks).upper().replace("*", "")))

    return records
```

## Core helper functions

```python
from __future__ import annotations

from dataclasses import dataclass
from anarci import anarci


@dataclass
class AnarciDomain:
    sequence_id: str
    chain_type: str | None
    scheme: str
    numbered: list[tuple[tuple[int, str], str]]
    raw_detail: object


def run_anarci(
    records: list[tuple[str, str]],
    scheme: str = "imgt",
) -> tuple[list[AnarciDomain], list[str]]:
    numbering, alignment_details, hit_tables = anarci(
        records,
        scheme=scheme,
        output=False,
    )

    domains = []
    failed_ids = []

    for (seq_id, _), seq_numbering, seq_detail in zip(
        records,
        numbering,
        alignment_details,
    ):
        if seq_numbering is None:
            failed_ids.append(seq_id)
            continue

        for domain_numbering, domain_start, domain_end in seq_numbering:
            chain_type = None

            if isinstance(seq_detail, list) and seq_detail:
                detail0 = seq_detail[0]
                if isinstance(detail0, dict):
                    chain_type = detail0.get("chain_type") or detail0.get("chain")

            domains.append(
                AnarciDomain(
                    sequence_id=seq_id,
                    chain_type=chain_type,
                    scheme=scheme,
                    numbered=domain_numbering,
                    raw_detail=seq_detail,
                )
            )

    return domains, failed_ids


def extract_region(
    numbered: list[tuple[tuple[int, str], str]],
    start: int,
    end: int,
) -> str:
    residues = []

    for (position_number, insertion_code), aa in numbered:
        if start <= position_number <= end and aa != "-":
            residues.append(aa)

    return "".join(residues)


def extract_imgt_regions(
    numbered: list[tuple[tuple[int, str], str]]
) -> dict[str, str]:
    return {
        region: extract_region(numbered, start, end)
        for region, (start, end) in IMGT_REGIONS.items()
    }
```

## Make an annotation table

```python
import pandas as pd


def make_annotation_table(domains: list[AnarciDomain]) -> pd.DataFrame:
    rows = []

    for domain in domains:
        regions = extract_imgt_regions(domain.numbered)

        row = {
            "sequence_id": domain.sequence_id,
            "scheme": domain.scheme,
            "chain_type": domain.chain_type,
        }

        for region in ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]:
            seq = regions[region]
            key = region.lower()
            row[key] = seq
            row[f"{key}_len"] = len(seq)

        rows.append(row)

    return pd.DataFrame(rows)
```

## End-to-end script skeleton

```python
import argparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-fasta", required=True)
    parser.add_argument("-o", "--output-tsv", required=True)
    parser.add_argument("--scheme", default="imgt")
    parser.add_argument("--failed-output", default=None)
    args = parser.parse_args()

    records = read_fasta(args.input_fasta)
    domains, failed_ids = run_anarci(records, scheme=args.scheme)

    table = make_annotation_table(domains)
    table.to_csv(args.output_tsv, sep="\t", index=False)

    if args.failed_output:
        with open(args.failed_output, "w") as handle:
            for seq_id in failed_ids:
                handle.write(f"{seq_id}\n")

    print(f"Input sequences: {len(records)}")
    print(f"Numbered domains: {len(domains)}")
    print(f"Failed sequences: {len(failed_ids)}")


if __name__ == "__main__":
    main()
```

Example:

```bash
python annotate_antibodies.py \
  -i antibodies.fasta \
  -o antibodies_imgt.tsv \
  --scheme imgt \
  --failed-output failed_ids.txt
```

## Output columns

Recommended TSV columns:

```text
sequence_id
scheme
chain_type
fr1
cdr1
fr2
cdr2
fr3
cdr3
fr4
fr1_len
cdr1_len
fr2_len
cdr2_len
fr3_len
cdr3_len
fr4_len
```

Only add germline V/J, species, or e-value columns if they are actually available from the parsed ANARCI output or another tool. Do not infer detailed V(D)J assignment from numbering alone.

## Position formatting

ANARCI positions are usually represented as:

```python
((position_number, insertion_code), amino_acid)
```

Use this helper for position-level tables:

```python
def format_position(position: tuple[int, str]) -> str:
    number, insertion = position
    insertion = insertion.strip()
    return f"{number}{insertion}" if insertion else str(number)
```

Preserve insertion codes such as `111A`, `111B`, etc.

## Quality checks

For production scripts:

* Validate amino-acid characters.
* Flag sequences with stop codons.
* Record failed ANARCI calls.
* Require non-empty CDR3 for normal VH/VL annotation.
* Be cautious with partial sequences.
* Be cautious with full Ig constructs containing signal peptide, constant region, tags, or linkers.
* Emit one row per recognized domain if multiple domains are found.

## When to use AbNumber

Use `abnumber` for small Python analyses or notebooks that need convenient FR/CDR access:

```python
from abnumber import Chain

chain = Chain("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSW...", scheme="imgt")
print(chain.cdr1_seq)
print(chain.cdr2_seq)
print(chain.cdr3_seq)
```

Use raw ANARCI for batch pipelines, exact control of failed sequences, and compatibility with existing ANARCI workflows.

## Rules for Codex

When generating ANARCI code:

* Default to IMGT.
* Keep input IDs unchanged.
* Separate successful and failed sequences.
* Emit one row per domain.
* Include `scheme` in every output table.
* Do not assume every sequence is a valid antibody.
* Do not infer V/J germline genes from numbering alone.
* Preserve insertion codes in position-level outputs.
* Keep parsing defensive against ANARCI version differences.
