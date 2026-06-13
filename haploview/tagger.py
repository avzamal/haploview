"""Tag-SNP selection (Haploview's pairwise Tagger / de Bakker et al. 2005).

A tag SNP captures another SNP when their r^2 is >= a cutoff (default 0.8)
within a distance window (default 500 kb). The pairwise algorithm greedily
selects, at each step, the SNP that captures the most still-uncaptured SNPs,
until every SNP is captured (a SNP always captures itself). Each captured SNP is
then assigned its "best tag" -- the chosen tag it has the highest r^2 with.

This is a faithful transcription of ``edu/mit/wi/tagger/Tagger.findTags``
(PAIRWISE_ONLY mode). Tie-breaks among equally good candidate tags are resolved
deterministically by genomic order (the original used hash-map order, which is
not reproducible); the resulting tag set is an equivalent valid cover.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .ld import LDTable
from .model import Dataset

DEFAULT_RSQ_CUTOFF = 0.8


@dataclass
class TagResult:
    tags: List[int]                                   # chosen tag SNP indices, in order
    tag_captures: Dict[int, List[int]]                # tag -> all SNP indices it can capture
    best_tag: Dict[int, int]                          # snp -> its best tag index
    best_rsq: Dict[int, float]                        # snp -> r^2 with its best tag
    untagged: List[int]                               # SNPs not captured by any tag
    rsq_cutoff: float = DEFAULT_RSQ_CUTOFF
    mean_rsq: float = 0.0
    percent_captured: int = 0


def _rsq(ld: LDTable, i: int, j: int) -> float:
    if i == j:
        return 1.0
    pair = ld.get(i, j)
    return pair.rsq if pair is not None else 0.0


def select_tags(dataset: Dataset, ld: LDTable,
                rsq_cutoff: float = DEFAULT_RSQ_CUTOFF) -> TagResult:
    n = dataset.n_markers
    capture = list(range(n))

    # all_tagged[t] = every SNP t can capture (r^2 >= cutoff, within LD window),
    # including t itself; tagged[t] is the working set of still-uncaptured ones.
    all_tagged: Dict[int, set] = {}
    for t in range(n):
        s = {t}
        for c in capture:
            if c != t and _rsq(ld, t, c) >= rsq_cutoff:
                s.add(c)
        all_tagged[t] = s
    tagged = {t: set(s) for t, s in all_tagged.items()}

    potential = set(range(n))
    sites_to_capture = set(capture)
    tags: List[int] = []

    while sites_to_capture:
        if not potential:
            break
        # pick the tag covering the most uncaptured sites; deterministic tie-break
        best = max(potential, key=lambda t: (len(tagged[t]), -t))
        if len(tagged[best]) == 0:
            potential.discard(best)
            continue
        newly = set(tagged[best])
        tags.append(best)
        potential.discard(best)
        for pt in list(potential):
            tagged[pt] -= newly
            if len(tagged[pt]) == 0:
                potential.discard(pt)
        sites_to_capture -= newly
        sites_to_capture.discard(best)

    # map every captured SNP to the chosen tags that can capture it; pick the best
    tag_captures = {t: sorted(all_tagged[t]) for t in tags}
    best_tag: Dict[int, int] = {}
    best_rsq: Dict[int, float] = {}
    for snp in capture:
        candidates = [t for t in tags if snp in all_tagged[t]]
        if not candidates:
            continue
        bt = max(candidates, key=lambda t: _rsq(ld, snp, t))
        best_tag[snp] = bt
        best_rsq[snp] = _rsq(ld, snp, bt)

    untagged = sorted(sites_to_capture)
    captured = [s for s in capture if s in best_tag]
    mean = sum(best_rsq[s] for s in captured) / len(captured) if captured else 0.0
    pct = int(len(captured) / len(capture) * 100) if capture else 0

    return TagResult(tags=tags, tag_captures=tag_captures, best_tag=best_tag,
                     best_rsq=best_rsq, untagged=untagged, rsq_cutoff=rsq_cutoff,
                     mean_rsq=mean, percent_captured=pct)
