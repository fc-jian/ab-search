---
name: biotite
description: Use for mmCIF reading, author and label chain mapping, assembly copies and residue contacts in the ab-search structural backends.
---

# Structure coordinates in ab-search

## Select the backend

| Backend | Coordinate and identity conventions |
| --- | --- |
| SAbDab, `score_structures.py` | `load_structure_variants()` reads cached text `.cif.gz`, requests model 1 and occupancy-based altloc selection, and tries the first biological assembly. SAbDab matching uses the deposited author chain IDs. The scorer can fall back to the asymmetric unit, recording the reason in `structure_warnings`. |
| Protein partners, `partner_structures.py` | `structure_rows()` visits listed biological assemblies and uses label chain and residue IDs to join SIFTS. It retains `auth_asym_id` and `sym_id` for chain-instance reporting. A missing assembly produces a review record. |

The deposited author residue key includes the chain, residue number and insertion code. Within an assembly, retain the assembly ID and distinguish symmetry copies. The same chain label can occur in multiple copies. The SAbDab helper currently filters by author chain alone, so check for this ambiguity before reusing it for a new assembly analysis.

## Measure contacts

For antibody scores, number the available full-chain sequence and project CDRs onto observed coordinates. `contact_residues_for_chain()` counts a residue when the CDR Cα distance is at most 8 Å **or** the minimum heavy-atom distance is at most 4.5 Å by default. The multichain penalty includes other contacting protein chains, subject to the configured homodimer exclusion. The partners backend uses its own heavy-atom cutoff and checks mapped residues and construct context.

For changes to selection or assembly handling, check an insertion-coded pair, two copies bearing one chain ID and the missing-assembly behavior of the affected backend. Use the cached compressed input format in I/O checks.

See the [Biotite assembly API](https://www.biotite-python.org/latest/apidoc/biotite.structure.io.pdbx.get_assembly.html) and [RCSB identifier conventions](https://www.rcsb.org/docs/general-help/identifiers-in-pdb).
