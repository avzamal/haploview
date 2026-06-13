"""Pairwise linkage disequilibrium.

Direct transcription of Haploview's ``HaploData.computeDPrime`` /
``count_haps`` / ``estimate_p`` (``edu/mit/wi/haploview/HaploData.java``):
a two-locus EM over unphased genotypes giving D', r^2, LOD and the
Gabriel-style 95% confidence interval on D'. Constants and rounding match the
original so output is identical at its printed precision.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .model import Dataset

LN10 = math.log(10.0)
TOLERANCE = 0.00000001

# indices into the four-haplotype arrays
AA, AB, BA, BB = 0, 1, 2, 3


def round_double(d: float, places: int) -> float:
    """Replicate Haploview's Util.roundDouble (Java Math.round, half-up)."""
    factor = 10 ** places
    return math.floor(d * factor + 0.5) / factor


@dataclass
class PairwiseLinkage:
    dprime: float
    lod: float
    rsq: float
    low_ci: float
    high_ci: float
    # haplotype frequencies in Haploview order: [AA, AB, BB, BA]
    freqs: List[float]


def _count_table(m1: np.ndarray, m2: np.ndarray):
    """Count obligate two-marker haplotypes + double-hets from minor dosages.

    Mirrors the per-chromosome accumulation in computeDPrime: homozygotes add
    two obligate haplotypes, single-marker hets add one haplotype to each of
    two cells, double-hets are ambiguous.
    """
    valid = (m1 != -1) & (m2 != -1)
    het1 = m1 == 1
    het2 = m2 == 1
    dh = int(np.count_nonzero(valid & het1 & het2))

    def c(mask):
        return int(np.count_nonzero(valid & mask))

    # neither marker het (both homozygous) -> two identical obligate haplotypes
    aa = 2 * c((~het1) & (~het2) & (m1 == 0) & (m2 == 0))
    ab = 2 * c((~het1) & (~het2) & (m1 == 0) & (m2 == 2))
    ba = 2 * c((~het1) & (~het2) & (m1 == 2) & (m2 == 0))
    bb = 2 * c((~het1) & (~het2) & (m1 == 2) & (m2 == 2))

    # marker1 het only: shares marker2 allele, splits over marker1
    e = c(het1 & (~het2) & (m2 == 0)); aa += e; ba += e
    e = c(het1 & (~het2) & (m2 == 2)); ab += e; bb += e
    # marker2 het only: shares marker1 allele, splits over marker2
    e = c((~het1) & het2 & (m1 == 0)); aa += e; ab += e
    e = c((~het1) & het2 & (m1 == 2)); ba += e; bb += e

    return aa, ab, ba, bb, dh


def compute_dprime(dataset: Dataset, pos1: int, pos2: int) -> Optional[PairwiseLinkage]:
    """LD statistics for the marker pair (pos1, pos2). None if monomorphic."""
    if dataset.markers[pos1].maf == 0 or dataset.markers[pos2].maf == 0:
        return None

    m1 = dataset.minor_dosage_row(pos1)
    m2 = dataset.minor_dosage_row(pos2)
    k_aa, k_ab, k_ba, k_bb, dh = _count_table(m1, m2)

    # monomorphic check: a margin is empty and there are no double hets
    r1 = k_aa + k_ab
    r2 = k_ba + k_bb
    c1 = k_aa + k_ba
    c2 = k_ab + k_bb
    if (r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0) and dh == 0:
        return PairwiseLinkage(1.0, 0.0, 0.0, 0.0, 0.0, [])

    known = [float(k_aa), float(k_ab), float(k_ba), float(k_bb)]
    total_chroms = known[AA] + known[AB] + known[BA] + known[BB] + 2.0 * dh
    pA1 = (known[AA] + known[AB] + dh) / total_chroms
    pB1 = 1.0 - pA1
    pA2 = (known[AA] + known[BA] + dh) / total_chroms
    pB2 = 1.0 - pA2

    prob = [0.0, 0.0, 0.0, 0.0]
    num_haps = [0.0, 0.0, 0.0, 0.0]

    def count_haps(em_round: int):
        num_haps[AA] = known[AA]
        num_haps[AB] = known[AB]
        num_haps[BA] = known[BA]
        num_haps[BB] = known[BB]
        if em_round > 0:
            denom = (prob[AA] * prob[BB]) + (prob[AB] * prob[BA])
            num_haps[AA] += dh * (prob[AA] * prob[BB]) / denom
            num_haps[BB] += dh * (prob[AA] * prob[BB]) / denom
            num_haps[AB] += dh * (prob[AB] * prob[BA]) / denom
            num_haps[BA] += dh * (prob[AB] * prob[BA]) / denom

    def estimate_p(const_prob: float):
        total = num_haps[AA] + num_haps[AB] + num_haps[BA] + num_haps[BB] + 4.0 * const_prob
        for x in range(4):
            prob[x] = (num_haps[x] + const_prob) / total
            if prob[x] < 1e-10:
                prob[x] = 1e-10

    # initial conditions: const_prob = 0.1 then one count/estimate step
    for x in range(4):
        prob[x] = 0.1
    count_haps(0)
    estimate_p(0.1)

    # EM iterations with const_prob = 0
    count = 1
    loglike = -999999999.0
    while count < 1000:
        oldloglike = loglike
        count_haps(count)
        loglike = (known[AA] * math.log(prob[AA]) + known[AB] * math.log(prob[AB])
                   + known[BA] * math.log(prob[BA]) + known[BB] * math.log(prob[BB])) / LN10 \
            + (dh * math.log(prob[AA] * prob[BB] + prob[AB] * prob[BA])) / LN10
        if abs(loglike - oldloglike) < TOLERANCE:
            break
        estimate_p(0.0)
        count += 1

    loglike1 = (known[AA] * math.log(prob[AA]) + known[AB] * math.log(prob[AB])
                + known[BA] * math.log(prob[BA]) + known[BB] * math.log(prob[BB])
                + dh * math.log(prob[AA] * prob[BB] + prob[AB] * prob[BA])) / LN10
    loglike0 = (known[AA] * math.log(pA1 * pA2) + known[AB] * math.log(pA1 * pB2)
                + known[BA] * math.log(pB1 * pA2) + known[BB] * math.log(pB1 * pB2)
                + dh * math.log(2 * pA1 * pA2 * pB1 * pB2)) / LN10

    num = prob[AA] * prob[BB] - prob[AB] * prob[BA]
    if num < 0:
        # flip to report positive D'
        prob[AA], prob[AB] = prob[AB], prob[AA]
        prob[BB], prob[BA] = prob[BA], prob[BB]
        pA2, pB2 = pB2, pA2
        num_haps[AA], num_haps[AB] = num_haps[AB], num_haps[AA]
        num_haps[BB], num_haps[BA] = num_haps[BA], num_haps[BB]
        known[AA], known[AB] = known[AB], known[AA]
        known[BB], known[BA] = known[BA], known[BB]
        num = prob[AA] * prob[BB] - prob[AB] * prob[BA]

    denom1 = (prob[AA] + prob[BA]) * (prob[BA] + prob[BB])
    denom2 = (prob[AA] + prob[AB]) * (prob[AB] + prob[BB])
    denom = denom1 if denom1 < denom2 else denom2
    dprime = num / denom
    rsq = num * num / (pA1 * pB1 * pA2 * pB2)

    # confidence bounds (Gabriel et al. 2002): posterior over D' on a flat prior
    lsurface = [0.0] * 101
    for i in range(101):
        dpr = i * 0.01
        tmpAA = dpr * denom + pA1 * pA2
        tmpAB = pA1 - tmpAA
        tmpBA = pA2 - tmpAA
        tmpBB = pB1 - tmpBA
        if i == 100:
            tmpAA = max(tmpAA, 1e-10)
            tmpAB = max(tmpAB, 1e-10)
            tmpBA = max(tmpBA, 1e-10)
            tmpBB = max(tmpBB, 1e-10)
        lsurface[i] = (known[AA] * math.log(tmpAA) + known[AB] * math.log(tmpAB)
                       + known[BA] * math.log(tmpBA) + known[BB] * math.log(tmpBB)
                       + dh * math.log(tmpAA * tmpBB + tmpAB * tmpBA)) / LN10

    total_prob = 0.0
    for i in range(101):
        lsurface[i] -= loglike1
        lsurface[i] = math.pow(10.0, lsurface[i])
        total_prob += lsurface[i]

    low_i = 0
    sum_prob = 0.0
    for i in range(101):
        sum_prob += lsurface[i]
        if sum_prob > 0.05 * total_prob and sum_prob - lsurface[i] < 0.05 * total_prob:
            low_i = i - 1
            break

    high_i = 0
    sum_prob = 0.0
    for i in range(100, -1, -1):
        sum_prob += lsurface[i]
        if sum_prob > 0.05 * total_prob and sum_prob - lsurface[i] < 0.05 * total_prob:
            high_i = i + 1
            break
    if high_i > 100:
        high_i = 100

    freqs = [prob[AA], prob[AB], prob[BB], prob[BA]]
    return PairwiseLinkage(
        round_double(dprime, 3),
        round_double(loglike1 - loglike0, 2),
        round_double(rsq, 3),
        low_i / 100.0,
        high_i / 100.0,
        freqs,
    )


class LDTable:
    """Pairwise LD over a marker set, windowed by physical distance.

    Replicates ``generateDPrimeTable`` / ``DPrimeTable``: a pair is computed
    only when the markers are within ``max_dist_bp`` (0 = no limit). Lookups
    beyond the window return ``None`` just like the original.
    """

    def __init__(self, dataset: Dataset, max_dist_bp: int = 500_000):
        self.dataset = dataset
        self.max_dist_bp = max_dist_bp
        n = dataset.n_markers
        self._rows: List[List[Optional[PairwiseLinkage]]] = [[] for _ in range(n)]
        self._window_end = list(range(n))  # last in-window marker index per row
        for x in range(n - 1):
            px = dataset.markers[x].position
            row: List[Optional[PairwiseLinkage]] = []
            last = x
            for y in range(x + 1, n):
                sep = dataset.markers[y].position - px
                if max_dist_bp > 0 and sep > max_dist_bp:
                    break
                row.append(compute_dprime(dataset, x, y))
                last = y
            self._rows[x] = row
            self._window_end[x] = last

    def get(self, x: int, y: int) -> Optional[PairwiseLinkage]:
        if x > y:
            x, y = y, x
        if x >= self.dataset.n_markers - 1:
            return None
        idx = y - x - 1
        row = self._rows[x]
        if 0 <= idx < len(row):
            return row[idx]
        return None

    def length(self, x: int) -> int:
        """Number of markers from x to the last in-window marker, inclusive."""
        return self._window_end[x] - x + 1
