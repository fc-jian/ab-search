---
name: biotite
description: Use Biotite for protein/RNA/DNA structure parsing, especially CIF/mmCIF/BCIF IO, AtomArray filtering, residue/atom selection, and residue-atom distance calculations.
---

# Biotite structure-analysis skill

Use this skill when writing Python code for protein/RNA/DNA structure parsing and geometric analysis with **Biotite**, especially for PDBx/mmCIF or BinaryCIF IO and residue/atom distance calculations.

## Scope

Prefer Biotite for:

* Reading `.cif`, `.mmcif`, `.pdbx`, and `.bcif` structure files.
* Converting PDBx/mmCIF content into `AtomArray` or `AtomArrayStack`.
* Filtering atoms by chain, residue ID, residue name, atom name, element, hetero flag, etc.
* Calculating distances between specific atoms, between all atoms in two residues, or between residue centroids/representative atoms.
* Writing modified structures back to mmCIF/PDBx when needed.

Do not use Biopython’s `Bio.PDB` unless the user explicitly asks for it or existing project code already depends on it.

## Imports

Use these imports by default:

```python
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
```

## Core Biotite concepts

Biotite represents structures as:

* `AtomArray`: one model, shape roughly `(n_atoms,)`.
* `AtomArrayStack`: multiple models/frames, shape roughly `(n_models, n_atoms)`.
* Annotation arrays: `array.chain_id`, `array.res_id`, `array.res_name`, `array.atom_name`, `array.element`, etc.
* Coordinates: `array.coord`, shape `(n_atoms, 3)` for `AtomArray`.

Use NumPy-style boolean masks for selections:

```python
mask = (
    (atoms.chain_id == "A")
    & (atoms.res_id == 145)
    & (atoms.atom_name == "CA")
)
selection = atoms[mask]
```

Single integer indexing returns an `Atom`; boolean/index-array slicing returns an `AtomArray`.

## Reading CIF/mmCIF/BCIF

For text mmCIF/PDBx:

```python
def load_cif(path: str | Path, model: int = 1) -> struc.AtomArray:
    path = Path(path)
    cif = pdbx.CIFFile.read(path)
    atoms = pdbx.get_structure(
        cif,
        model=model,
        altloc="occupancy",
        use_author_fields=True,
        include_bonds=False,
    )
    return atoms
```

For BinaryCIF:

```python
def load_bcif(path: str | Path, model: int = 1) -> struc.AtomArray:
    path = Path(path)
    bcif = pdbx.BinaryCIFFile.read(path)
    atoms = pdbx.get_structure(
        bcif,
        model=model,
        altloc="occupancy",
        use_author_fields=True,
        include_bonds=False,
    )
    return atoms
```

Default choices:

* Use `model=1` unless the user explicitly asks for all models.
* Use `altloc="occupancy"` for real structural analysis, because it chooses the highest-occupancy alternate location per residue.
* Use `use_author_fields=True` when the user refers to PDB-style chain IDs/residue numbering.
* Use `include_bonds=True` only when bond topology is needed. For distance calculations, coordinates are sufficient.

If the code must support both `.cif` and `.bcif`:

```python
def load_structure(path: str | Path, model: int = 1) -> struc.AtomArray:
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in {".cif", ".mmcif", ".pdbx"}:
        file = pdbx.CIFFile.read(path)
    elif suffix == ".bcif":
        file = pdbx.BinaryCIFFile.read(path)
    else:
        raise ValueError(f"Unsupported structure format: {path.suffix}")

    return pdbx.get_structure(
        file,
        model=model,
        altloc="occupancy",
        use_author_fields=True,
        include_bonds=False,
    )
```

## Basic cleanup

For protein-only residue atom analysis, usually remove non-amino-acid atoms:

```python
protein = atoms[struc.filter_amino_acids(atoms)]
```

To remove waters:

```python
no_water = atoms[atoms.res_name != "HOH"]
```

To remove hydrogens:

```python
heavy = atoms[atoms.element != "H"]
```

For most residue-distance tasks, use heavy atoms unless hydrogen-specific distances are requested.

## Atom selection helper

Write small explicit helper functions instead of burying masks inline.

```python
def select_atoms(
    atoms: struc.AtomArray,
    chain_id: str | None = None,
    res_id: int | None = None,
    atom_name: str | None = None,
    res_name: str | None = None,
    include_hydrogen: bool = True,
) -> struc.AtomArray:
    mask = np.ones(atoms.array_length(), dtype=bool)

    if chain_id is not None:
        mask &= atoms.chain_id == chain_id
    if res_id is not None:
        mask &= atoms.res_id == res_id
    if atom_name is not None:
        mask &= atoms.atom_name == atom_name
    if res_name is not None:
        mask &= atoms.res_name == res_name
    if not include_hydrogen:
        mask &= atoms.element != "H"

    return atoms[mask]
```

## Strict single-atom lookup

Use this when calculating a distance between named atoms such as `A:145:CA` and `B:32:CB`.

```python
def get_single_atom(
    atoms: struc.AtomArray,
    chain_id: str,
    res_id: int,
    atom_name: str,
) -> struc.Atom:
    sel = select_atoms(
        atoms,
        chain_id=chain_id,
        res_id=res_id,
        atom_name=atom_name,
    )

    if sel.array_length() == 0:
        raise ValueError(f"Atom not found: chain={chain_id}, res_id={res_id}, atom={atom_name}")
    if sel.array_length() > 1:
        raise ValueError(
            f"Ambiguous atom selection: chain={chain_id}, res_id={res_id}, "
            f"atom={atom_name}, matches={sel.array_length()}"
        )

    return sel[0]
```

## Distance between two named atoms

Prefer `struc.distance()` over manual Euclidean distance.

```python
def atom_distance(
    atoms: struc.AtomArray,
    atom1: tuple[str, int, str],
    atom2: tuple[str, int, str],
) -> float:
    chain1, res1, name1 = atom1
    chain2, res2, name2 = atom2

    a1 = get_single_atom(atoms, chain1, res1, name1)
    a2 = get_single_atom(atoms, chain2, res2, name2)

    return float(struc.distance(a1, a2))
```

Example:

```python
atoms = load_structure("input.cif")
d = atom_distance(atoms, ("A", 145, "CA"), ("B", 32, "CB"))
print(f"{d:.3f} Å")
```

## Distance between two residues: minimum heavy-atom distance

For residue-residue contact analysis, use the minimum pairwise heavy-atom distance by default.

```python
def residue_min_distance(
    atoms: struc.AtomArray,
    residue1: tuple[str, int],
    residue2: tuple[str, int],
    include_hydrogen: bool = False,
) -> float:
    chain1, res1 = residue1
    chain2, res2 = residue2

    r1 = select_atoms(
        atoms,
        chain_id=chain1,
        res_id=res1,
        include_hydrogen=include_hydrogen,
    )
    r2 = select_atoms(
        atoms,
        chain_id=chain2,
        res_id=res2,
        include_hydrogen=include_hydrogen,
    )

    if r1.array_length() == 0:
        raise ValueError(f"Residue not found: chain={chain1}, res_id={res1}")
    if r2.array_length() == 0:
        raise ValueError(f"Residue not found: chain={chain2}, res_id={res2}")

    # Broadcasting:
    # r1.coord[:, None, :] has shape (n1, 1, 3)
    # r2.coord[None, :, :] has shape (1, n2, 3)
    # result has shape (n1, n2)
    dist_matrix = struc.distance(r1.coord[:, None, :], r2.coord[None, :, :])
    return float(np.min(dist_matrix))
```

Example:

```python
atoms = load_structure("input.cif")
d_min = residue_min_distance(atoms, ("A", 145), ("B", 32))
print(f"Minimum heavy-atom distance: {d_min:.3f} Å")
```

## Distance matrix between residues

Use this when the user needs the actual atom-atom distance matrix.

```python
def residue_distance_matrix(
    atoms: struc.AtomArray,
    residue1: tuple[str, int],
    residue2: tuple[str, int],
    include_hydrogen: bool = False,
) -> tuple[np.ndarray, list[str], list[str]]:
    chain1, res1 = residue1
    chain2, res2 = residue2

    r1 = select_atoms(atoms, chain_id=chain1, res_id=res1, include_hydrogen=include_hydrogen)
    r2 = select_atoms(atoms, chain_id=chain2, res_id=res2, include_hydrogen=include_hydrogen)

    if r1.array_length() == 0:
        raise ValueError(f"Residue not found: chain={chain1}, res_id={res1}")
    if r2.array_length() == 0:
        raise ValueError(f"Residue not found: chain={chain2}, res_id={res2}")

    distances = struc.distance(r1.coord[:, None, :], r2.coord[None, :, :])

    labels1 = [f"{chain1}:{res1}:{name}" for name in r1.atom_name]
    labels2 = [f"{chain2}:{res2}:{name}" for name in r2.atom_name]

    return distances, labels1, labels2
```

## Closest atom pair between two residues

```python
def closest_atoms_between_residues(
    atoms: struc.AtomArray,
    residue1: tuple[str, int],
    residue2: tuple[str, int],
    include_hydrogen: bool = False,
) -> dict[str, object]:
    chain1, res1 = residue1
    chain2, res2 = residue2

    r1 = select_atoms(atoms, chain_id=chain1, res_id=res1, include_hydrogen=include_hydrogen)
    r2 = select_atoms(atoms, chain_id=chain2, res_id=res2, include_hydrogen=include_hydrogen)

    if r1.array_length() == 0:
        raise ValueError(f"Residue not found: chain={chain1}, res_id={res1}")
    if r2.array_length() == 0:
        raise ValueError(f"Residue not found: chain={chain2}, res_id={res2}")

    distances = struc.distance(r1.coord[:, None, :], r2.coord[None, :, :])
    i, j = np.unravel_index(np.argmin(distances), distances.shape)

    return {
        "chain1": chain1,
        "res_id1": res1,
        "res_name1": str(r1.res_name[i]),
        "atom_name1": str(r1.atom_name[i]),
        "chain2": chain2,
        "res_id2": res2,
        "res_name2": str(r2.res_name[j]),
        "atom_name2": str(r2.atom_name[j]),
        "distance_angstrom": float(distances[i, j]),
    }
```

## Contact search between two chains

For simple contact tables, calculate all residue pairs with minimum heavy-atom distance below a cutoff.

```python
def residue_contacts_between_chains(
    atoms: struc.AtomArray,
    chain1: str,
    chain2: str,
    cutoff: float = 4.0,
) -> list[dict[str, object]]:
    atoms = atoms[struc.filter_amino_acids(atoms)]
    atoms = atoms[atoms.element != "H"]

    a1 = atoms[atoms.chain_id == chain1]
    a2 = atoms[atoms.chain_id == chain2]

    if a1.array_length() == 0:
        raise ValueError(f"No atoms found for chain {chain1}")
    if a2.array_length() == 0:
        raise ValueError(f"No atoms found for chain {chain2}")

    residues1 = sorted(set(int(x) for x in a1.res_id))
    residues2 = sorted(set(int(x) for x in a2.res_id))

    contacts: list[dict[str, object]] = []

    for res1 in residues1:
        for res2 in residues2:
            d = residue_min_distance(atoms, (chain1, res1), (chain2, res2))
            if d <= cutoff:
                contacts.append(
                    {
                        "chain1": chain1,
                        "res_id1": res1,
                        "chain2": chain2,
                        "res_id2": res2,
                        "min_distance_angstrom": d,
                    }
                )

    return contacts
```

For large structures, avoid this O(n_residue_pairs) implementation and use a spatial-neighbor strategy such as `CellList` if performance becomes a bottleneck.

## Writing CIF/mmCIF

```python
def write_cif(atoms: struc.AtomArray, path: str | Path) -> None:
    path = Path(path)
    cif = pdbx.CIFFile()
    pdbx.set_structure(cif, atoms)
    cif.write(path)
```

If preserving non-coordinate mmCIF metadata is important, do not create a new empty `CIFFile`; read the original file, modify the structure, then call `set_structure()` on the original object.

## Common pitfalls

1. **Model handling**

   * `get_structure(file)` without `model` returns an `AtomArrayStack`.
   * Most residue-distance helper functions should use `model=1` and return an `AtomArray`.

2. **Author vs label residue numbering**

   * Use `use_author_fields=True` when the user gives residue IDs as seen in PDB/mmCIF viewers or publications.
   * Use `use_author_fields=False` only when the task explicitly requires canonical PDBx label IDs.

3. **Alternate conformations**

   * Use `altloc="occupancy"` for distance calculations.
   * Avoid `altloc="all"` unless the user wants all alternate locations; it can duplicate atoms and make selections ambiguous.

4. **Hydrogens**

   * Most deposited structures lack hydrogens.
   * For residue-residue distance or contact analysis, default to heavy atoms.

5. **Insertion codes**

   * PDB-style residue identity may require insertion codes in some structures.
   * If residue selection by `res_id` is ambiguous, inspect available annotation arrays and include insertion-code fields if present.

6. **Chain IDs**

   * With `use_author_fields=True`, chain IDs should match common PDB viewer chain labels.
   * With label fields, chain IDs may differ from author chain IDs.

7. **Units**

   * Coordinates and distances are in Ångström for PDB/mmCIF structures.

## Minimal script template

```python
from __future__ import annotations

from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx


def load_structure(path: str | Path, model: int = 1) -> struc.AtomArray:
    path = Path(path)

    if path.suffix.lower() in {".cif", ".mmcif", ".pdbx"}:
        file = pdbx.CIFFile.read(path)
    elif path.suffix.lower() == ".bcif":
        file = pdbx.BinaryCIFFile.read(path)
    else:
        raise ValueError(f"Unsupported structure format: {path.suffix}")

    return pdbx.get_structure(
        file,
        model=model,
        altloc="occupancy",
        use_author_fields=True,
        include_bonds=False,
    )


def select_atoms(
    atoms: struc.AtomArray,
    chain_id: str | None = None,
    res_id: int | None = None,
    atom_name: str | None = None,
    include_hydrogen: bool = True,
) -> struc.AtomArray:
    mask = np.ones(atoms.array_length(), dtype=bool)

    if chain_id is not None:
        mask &= atoms.chain_id == chain_id
    if res_id is not None:
        mask &= atoms.res_id == res_id
    if atom_name is not None:
        mask &= atoms.atom_name == atom_name
    if not include_hydrogen:
        mask &= atoms.element != "H"

    return atoms[mask]


def get_single_atom(
    atoms: struc.AtomArray,
    chain_id: str,
    res_id: int,
    atom_name: str,
) -> struc.Atom:
    sel = select_atoms(atoms, chain_id=chain_id, res_id=res_id, atom_name=atom_name)

    if sel.array_length() == 0:
        raise ValueError(f"Atom not found: chain={chain_id}, res_id={res_id}, atom={atom_name}")
    if sel.array_length() > 1:
        raise ValueError(
            f"Ambiguous atom selection: chain={chain_id}, res_id={res_id}, "
            f"atom={atom_name}, matches={sel.array_length()}"
        )

    return sel[0]


def atom_distance(
    atoms: struc.AtomArray,
    atom1: tuple[str, int, str],
    atom2: tuple[str, int, str],
) -> float:
    a1 = get_single_atom(atoms, *atom1)
    a2 = get_single_atom(atoms, *atom2)
    return float(struc.distance(a1, a2))


def residue_min_distance(
    atoms: struc.AtomArray,
    residue1: tuple[str, int],
    residue2: tuple[str, int],
) -> float:
    chain1, res1 = residue1
    chain2, res2 = residue2

    r1 = select_atoms(atoms, chain_id=chain1, res_id=res1, include_hydrogen=False)
    r2 = select_atoms(atoms, chain_id=chain2, res_id=res2, include_hydrogen=False)

    if r1.array_length() == 0:
        raise ValueError(f"Residue not found: chain={chain1}, res_id={res1}")
    if r2.array_length() == 0:
        raise ValueError(f"Residue not found: chain={chain2}, res_id={res2}")

    distances = struc.distance(r1.coord[:, None, :], r2.coord[None, :, :])
    return float(np.min(distances))


if __name__ == "__main__":
    atoms = load_structure("input.cif")

    ca_distance = atom_distance(atoms, ("A", 145, "CA"), ("B", 32, "CA"))
    print(f"CA-CA distance: {ca_distance:.3f} Å")

    min_distance = residue_min_distance(atoms, ("A", 145), ("B", 32))
    print(f"Residue minimum heavy-atom distance: {min_distance:.3f} Å")
```

## Coding rules for Codex

When generating Biotite code:

* Prefer explicit helper functions over one-off inline masks.
* Always validate empty selections and ambiguous single-atom selections.
* Use `model=1` unless the user requests multi-model analysis.
* Use `altloc="occupancy"` unless the user requests a different alternate-location policy.
* Use `use_author_fields=True` for residue numbers/chains supplied by humans.
* Default residue-residue distances to minimum heavy-atom distance.
* Report distances in Å.
* Add concise tests when editing a codebase: load a small fixture CIF, select known atoms, assert non-empty selections, and assert distance output is finite and positive.
