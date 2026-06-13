"""Haplotype block detection.

Faithful transcription of ``edu/mit/wi/haploview/FindBlocks.java`` for all three
methods: Gabriel et al. (2002) confidence-interval blocks (default), the
four-gamete test, and the solid-spine-of-LD method. All constants are the
originals.
"""

from __future__ import annotations

from typing import List

from .ld import LDTable
from .model import Dataset

CUT_HIGH_CI = 0.98
CUT_LOW_CI = 0.70
MAF_THRESH = 0.05
CUT_LOW_CI_VAR = [0, 0, 0.80, 0.50, 0.50]
MAX_DIST = [0, 0, 20000, 30000, 1000000]
REC_HIGH_CI = 0.90
INFORM_FRAC = 0.95
FOUR_GAMETE_CUTOFF = 0.01
SPINE_DP = 0.80


def _sorted_pairs(pairs):
    """Stable ascending sort by separation, then full reverse (as in Java)."""
    pairs = sorted(pairs, key=lambda p: p[2])  # stable on insertion order
    pairs.reverse()
    return pairs


def do_gabriel(dataset: Dataset, dprime: LDTable) -> List[List[int]]:
    n = dataset.n_markers
    skip = [dataset.markers[x].maf < MAF_THRESH for x in range(n)]

    strong_pairs = []
    for x in range(n - 1):
        for y in range(x + 1, n):
            pair = dprime.get(x, y)
            if pair is None:
                continue
            if skip[x] or skip[y]:
                continue
            if pair.lod < -90:
                continue
            if pair.high_ci < CUT_HIGH_CI or pair.low_ci < CUT_LOW_CI:
                continue
            sep = abs(dataset.markers[y].position - dataset.markers[x].position)
            strong_pairs.append((x, y, sep))

    strong_pairs = _sorted_pairs(strong_pairs)

    blocks: List[List[int]] = []
    used = [False] * (n + 1)
    for first, last, sep in strong_pairs:
        num_strong = num_rec = num_in_group = 0
        this_block: List[int] = []
        if used[first] or used[last]:
            continue
        for x in range(first, last + 1):
            if not skip[x]:
                num_in_group += 1
        if num_in_group < 4 and sep > MAX_DIST[num_in_group]:
            continue

        this_block.append(first)
        for y in range(first + 1, last + 1):
            if skip[y]:
                continue
            this_block.append(y)
            for x in range(first, y):
                if skip[x]:
                    continue
                pair = dprime.get(x, y)
                if pair is None:
                    continue
                if pair.lod < -90:
                    continue
                if pair.lod == 0 and pair.low_ci == 0 and pair.high_ci == 0:
                    continue
                if num_in_group < 5:
                    if pair.low_ci > CUT_LOW_CI_VAR[num_in_group] and pair.high_ci >= CUT_HIGH_CI:
                        num_strong += 1
                else:
                    if pair.low_ci > CUT_LOW_CI and pair.high_ci >= CUT_HIGH_CI:
                        num_strong += 1
                if pair.high_ci < REC_HIGH_CI:
                    num_rec += 1

        if num_in_group > 3:
            if num_strong + num_rec < 6:
                continue
        elif num_in_group > 2:
            if num_strong + num_rec < 3:
                continue
        else:
            if num_strong + num_rec < 1:
                continue

        if (num_strong / (num_strong + num_rec)) > INFORM_FRAC:
            _insert_sorted(blocks, this_block)
            for u in range(first, last + 1):
                used[u] = True
    return blocks


def do_four_gamete(dataset: Dataset, dprime: LDTable) -> List[List[int]]:
    n = dataset.n_markers
    strong_pairs = []
    for x in range(n - 1):
        for y in range(x + 1, n):
            pair = dprime.get(x, y)
            if pair is None:
                continue
            num_gam = sum(1 for f in pair.freqs if f > FOUR_GAMETE_CUTOFF + 1e-8)
            if num_gam > 3:
                continue
            strong_pairs.append((x, y, y - x - 1))

    strong_pairs = _sorted_pairs(strong_pairs)

    blocks: List[List[int]] = []
    used = [False] * (n + 1)
    for first, last, sep in strong_pairs:
        if used[first] or used[last]:
            continue
        is_block = True
        for y in range(first + 1, last + 1):
            for x in range(first, y):
                pair = dprime.get(x, y)
                if pair is None:
                    continue
                num_gam = sum(1 for f in pair.freqs if f > FOUR_GAMETE_CUTOFF + 1e-8)
                if num_gam > 3:
                    is_block = False
                    break
            if not is_block:
                break
        if is_block:
            _insert_sorted(blocks, list(range(first, last + 1)))
            for u in range(first, last + 1):
                used[u] = True
    return blocks


def do_spine(dataset: Dataset, dprime: LDTable) -> List[List[int]]:
    n = dataset.n_markers
    blocks: List[List[int]] = []
    # verticalExtent / horizontalExtent persist across markers in the original
    vertical_extent = 0
    horizontal_extent = 0
    i = 0
    while i < n:
        baddies = 0
        for j in range(i + 1, i + dprime.length(i)):
            pair = dprime.get(i, j)
            if pair is None:
                continue
            if pair.dprime < SPINE_DP:
                if baddies < 1:
                    baddies += 1
                else:
                    vertical_extent = j - 1
                    break
            vertical_extent = j
        m = vertical_extent
        while m > i:
            for k in range(i, m):
                pair = dprime.get(k, m)
                if pair is None:
                    continue
                if pair.dprime < SPINE_DP:
                    if baddies < 1:
                        baddies += 1
                    else:
                        break
                horizontal_extent = k + 1
            if horizontal_extent == m:
                blocks.append(list(range(i, m + 1)))
                i = m
            m -= 1
        i += 1
    return blocks


def _insert_sorted(blocks: List[List[int]], new_block: List[int]):
    """Insert a block keeping the list ordered by first marker index."""
    if not blocks:
        blocks.append(new_block)
        return
    for b in range(len(blocks)):
        if new_block[0] < blocks[b][0]:
            blocks.insert(b, new_block)
            return
    blocks.append(new_block)


METHODS = {
    "gabriel": do_gabriel,
    "4gam": do_four_gamete,
    "spine": do_spine,
}


def find_blocks(method: str, dataset: Dataset, dprime: LDTable) -> List[List[int]]:
    try:
        fn = METHODS[method]
    except KeyError:
        raise ValueError(f"unknown block method: {method!r}")
    return fn(dataset, dprime)
