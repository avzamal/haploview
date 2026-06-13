"""Within-block multi-marker haplotype frequency estimation (EM).

Transcribes the core of ``edu/mit/wi/haploview/EM.java`` (without the trio /
association machinery, which is irrelevant to unrelated VCF samples):

* a block is split into sub-segments of <=8 markers (base-8 partitioning);
* a standard phase EM (PSEUDOCOUNT 0.1, 20 iterations) is run per segment,
  keeping haplotypes with frequency > 0.001;
* segment haplotypes are ligated and a second EM produces final block
  haplotype frequencies;
* each individual is then assigned its most likely diplotype, used to genotype
  samples in the new haploblock VCF.

Haplotype alleles are coded 0 = major allele, 1 = minor allele at each locus,
packed into an integer bitmask (bit j = locus j).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Optional, Tuple

import numpy as np

from .model import Dataset

PSEUDOCOUNT = 0.1
MISSING_LIMIT = 4
MAX_ITER = 20
KEEP_THRESH = 0.001


def _block_sizes(n: int) -> List[int]:
    if n < 9:
        return [n]
    ones = n % 8
    eights = (n - ones) // 8
    if ones == 0:
        return [8] * eights
    sizes = [8] * (eights - 1)
    second_last = (8 + ones) // 2
    sizes.append(second_last)
    sizes.append(8 + ones - second_last)
    return sizes


def _locus_allele_pairs(code: int) -> List[Tuple[int, int]]:
    """Allowed (hap1, hap2) allele pairs at a locus given minor-dosage code."""
    if code == 0:
        return [(0, 0)]
    if code == 2:
        return [(1, 1)]
    if code == 1:
        return [(0, 1), (1, 0)]
    # missing: genotype unknown
    return [(0, 0), (0, 1), (1, 0), (1, 1)]


def _reconstructions(codes: List[int]) -> List[Tuple[int, int]]:
    """Distinct unordered (h1, h2) haplotype-mask reconstructions for one individual."""
    seen = set()
    out: List[Tuple[int, int]] = []
    for combo in product(*[_locus_allele_pairs(c) for c in codes]):
        h1 = h2 = 0
        for j, (x, y) in enumerate(combo):
            h1 |= x << j
            h2 |= y << j
        key = (h1, h2) if h1 <= h2 else (h2, h1)
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def _run_em(individuals: List[List[Tuple[int, int]]], num_poss: int) -> Dict[int, float]:
    """Phase EM over a list of individuals' reconstruction lists."""
    prob: Dict[int, float] = {}

    def get(h):
        return prob.get(h, PSEUDOCOUNT)

    total = num_poss * PSEUDOCOUNT
    for recs in individuals:
        if len(recs) == 1:
            h1, h2 = recs[0]
            prob[h1] = get(h1) + 1.0
            prob[h2] = get(h2) + 1.0
            total += 2.0
    for h in list(prob):
        prob[h] = prob[h] / total
    default = PSEUDOCOUNT / total

    def getn(p, h, d):
        return p.get(h, d)

    for _ in range(MAX_ITER):
        # E-step
        weighted = []
        for recs in individuals:
            ps = []
            s = 0.0
            for (h1, h2) in recs:
                p = getn(prob, h1, default) * getn(prob, h2, default)
                ps.append(p)
                s += p
            if s <= 0:
                s = 1.0
            weighted.append([(recs[k][0], recs[k][1], ps[k] / s) for k in range(len(recs))])
        # M-step
        newprob: Dict[int, float] = {}
        total = num_poss * 1e-10
        for recs in weighted:
            for (h1, h2, p) in recs:
                newprob[h1] = newprob.get(h1, 1e-10) + p
                newprob[h2] = newprob.get(h2, 1e-10) + p
                total += 2.0 * p
        for h in list(newprob):
            newprob[h] /= total
        prob = newprob
        default = 1e-10 / total
    return prob


@dataclass
class BlockHaplotypes:
    marker_indices: List[int]          # indices into the Dataset
    haplotypes: List[List[int]]        # each is allele list (0=major,1=minor) per locus
    frequencies: List[float]
    # per-sample assigned diplotype as (hap_index_a, hap_index_b) or None
    sample_diplotypes: List[Optional[Tuple[int, int]]]


def estimate_block_haplotypes(dataset: Dataset, marker_indices: List[int]) -> BlockHaplotypes:
    n = len(marker_indices)
    sizes = _block_sizes(n)
    block_len = n

    # minor-dosage codes per included marker, shape (n, n_samples)
    code_rows = np.vstack([dataset.minor_dosage_row(mi) for mi in marker_indices])

    # decide which individuals are usable (missing-data filter), per Haploview
    seg_bounds = []
    s = 0
    for sz in sizes:
        seg_bounds.append((s, s + sz))
        s += sz

    usable = []
    per_ind_codes = []
    for j in range(dataset.n_samples):
        col = code_rows[:, j]
        total_missing = int(np.count_nonzero(col == -1))
        too_many = False
        for (a, b) in seg_bounds:
            if int(np.count_nonzero(col[a:b] == -1)) >= MISSING_LIMIT:
                too_many = True
                break
        if not too_many and total_missing <= 1 + block_len // 3:
            usable.append(j)
            per_ind_codes.append([int(x) for x in col])

    # per-segment EM
    seg_haps: List[List[int]] = []        # kept hap masks per segment
    seg_recons: List[List[List[Tuple[int, int]]]] = []
    for (a, b) in seg_bounds:
        seg_codes = [codes[a:b] for codes in per_ind_codes]
        recons = [_reconstructions(c) for c in seg_codes]
        num_poss = 2 ** (b - a)
        prob = _run_em(recons, num_poss)
        kept = sorted(h for h, p in prob.items() if p > KEEP_THRESH)
        if not kept:
            kept = sorted(prob, key=lambda h: -prob[h])[:1]
        seg_haps.append(kept)
        seg_recons.append(recons)

    # ligation: build super-reconstructions consistent across segments
    super_inds: List[List[Tuple[int, int]]] = []
    seg_keepset = [set(k) for k in seg_haps]
    offsets = [a for (a, b) in seg_bounds]
    for idx in range(len(usable)):
        # for each segment, the kept reconstructions
        per_seg_choices = []
        ok = True
        for s_i, (a, b) in enumerate(seg_bounds):
            choices = [(h1, h2) for (h1, h2) in seg_recons[s_i][idx]
                       if h1 in seg_keepset[s_i] and h2 in seg_keepset[s_i]]
            if not choices:
                # fall back to all reconstructions if kept set excludes them
                choices = seg_recons[s_i][idx]
            per_seg_choices.append(choices)
        # combine segments into full-length haplotype-pair reconstructions
        combos = []
        seen = set()
        for combo in product(*per_seg_choices):
            full_h1 = 0
            full_h2 = 0
            for s_i, (h1, h2) in enumerate(combo):
                full_h1 |= h1 << offsets[s_i]
                full_h2 |= h2 << offsets[s_i]
            key = (full_h1, full_h2) if full_h1 <= full_h2 else (full_h2, full_h1)
            if key not in seen:
                seen.add(key)
                combos.append(key)
        super_inds.append(combos)

    num_poss_full = 2 ** n
    full_prob = _run_em(super_inds, num_poss_full)

    haps = sorted(h for h, p in full_prob.items() if p > KEEP_THRESH)
    if not haps:
        haps = sorted(full_prob, key=lambda h: -full_prob[h])[:1]
    freqs = [full_prob[h] for h in haps]
    hap_index = {h: i for i, h in enumerate(haps)}

    # assign each usable sample its most likely diplotype among kept haplotypes
    sample_dip: List[Optional[Tuple[int, int]]] = [None] * dataset.n_samples
    default = 1e-12
    for k, j in enumerate(usable):
        best = None
        best_p = -1.0
        for (h1, h2) in super_inds[k]:
            if h1 not in hap_index or h2 not in hap_index:
                continue
            p = full_prob.get(h1, default) * full_prob.get(h2, default)
            if p > best_p:
                best_p = p
                best = (hap_index[h1], hap_index[h2])
        sample_dip[j] = best

    hap_lists = [[(h >> j) & 1 for j in range(n)] for h in haps]
    return BlockHaplotypes(marker_indices, hap_lists, freqs, sample_dip)
