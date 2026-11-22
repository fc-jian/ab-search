"""Conservative transfer of an observed binding site to a fixed protein query."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import islice

from Bio.Align import PairwiseAligner


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
PROTEIN_ALPHABET = AMINO_ACIDS | frozenset("BJOUXZ")
# The initial cross-construct evidence panel covers explicitly identified host proteins.
HOSTS = {"human", "homo sapiens", "mouse", "mus musculus", "rat", "rattus norvegicus",
         "chicken", "gallus gallus", "bovine", "bos taurus"}


@dataclass(frozen=True)
class Match:
    kind: str
    mapping: dict[int, int]  # reference -> query, both zero-based
    identity: float
    ambiguous: bool = False


def sequence_match(query: str, reference: str, *, allow_local: bool, min_identity: float) -> Match | None:
    if not query or not reference:
        return None
    if query == reference:
        return Match("exact", dict(enumerate(range(len(query)))), 1.0)
    if query in reference or reference in query:
        shorter, longer = sorted((query, reference), key=len)
        start = longer.index(shorter)
        query_start, ref_start = (0, start) if shorter == query else (start, 0)
        mapping = {ref_start + i: query_start + i for i in range(len(shorter))}
        return Match("fragment" if shorter == query else "contains_reference", mapping, 1.0,
                     longer.find(shorter, start + 1) != -1)
    if not allow_local:
        return None

    # An unambiguous end-to-end alignment can retain a terminal substitution
    # that local alignment would clip, including on a short assayed peptide.
    aligner = PairwiseAligner(mode="local", match_score=2, mismatch_score=-2,
                              open_gap_score=-6, extend_gap_score=-1)
    if len(reference) == len(query) and len(query) >= 8:
        same = sum(r == q and r in AMINO_ACIDS for r, q in zip(reference, query))
        if same / len(query) >= max(0.90, min_identity):
            aligner.mode = "global"
            full = list(islice(aligner.align(reference, query), 2))
            if len(full) == 1 and full[0].shape[1] == len(query):
                return Match("full_length", dict(enumerate(range(len(query)))), same / len(query))
            aligner.mode = "local"
    # Internal gaps count in the local identity denominator. Tied placements
    # do not establish a unique site.
    alignments = list(islice(aligner.align(reference, query), 2))
    if not alignments:
        return None
    ref_indices, query_indices = alignments[0].indices
    mapping = {int(r): int(q) for r, q in zip(ref_indices, query_indices) if r >= 0 and q >= 0}
    matches = sum(reference[r] == query[q] and reference[r] in AMINO_ACIDS for r, q in mapping.items())
    identity = matches / len(ref_indices)
    if len(mapping) < 12 or identity < min_identity:
        return None
    return Match("local", mapping, identity, len(alignments) > 1)


def assess_site(query: str, reference: str, match: Match, positions: list[int],
                max_substitutions: int, min_site_identity: float, *, linear: bool) -> dict:
    """Positions are 1-based on the specified reference. Missing sites never pass."""
    metrics = {"site_size": len(positions)}
    if not positions:
        return {**metrics, "reason": "site_not_localized"}
    if match.ambiguous:
        return {**metrics, "reason": "ambiguous_alignment"}
    indices = [p - 1 for p in positions]
    pairs = [(p, match.mapping[p]) for p in indices if p in match.mapping]
    covered = len(pairs)
    metrics.update(site_covered=covered, site_coverage=covered / len(indices))
    if covered != len(indices):
        return {**metrics, "reason": "site_not_fully_covered"}
    if any(reference[p] not in AMINO_ACIDS or query[q] not in AMINO_ACIDS for p, q in pairs):
        return {**metrics, "reason": "unknown_site_residue"}
    if linear and pairs[-1][1] - pairs[0][1] + 1 != covered:
        return {**metrics, "reason": "insertion_in_mapped_region"}
    changes = sum(reference[p] != query[q] for p, q in pairs)
    identity = 1 - changes / covered
    metrics.update(site_substitutions=changes, site_identity=identity)
    if changes > max_substitutions or identity < min_site_identity:
        return {**metrics, "reason": "site_difference_limit"}
    return {**metrics, "tier": "site_similar" if changes else "site_conserved"}
