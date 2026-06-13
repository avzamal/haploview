"""Marker quality-control: Hardy-Weinberg exact test, call rate, MAF.

The HWE p-value is a direct transcription of Haploview's ``hwCalculate``
(``edu/mit/wi/pedfile/CheckData.java``), itself the Wigginton & Abecasis (2005)
exact test, so results match the original tool.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

from .model import Dataset, Marker


def hwe_exact_p(obs_het: int, obs_hom1: int, obs_hom2: int) -> float:
    """Exact two-sided Hardy-Weinberg p-value (Wigginton & Abecasis 2005).

    Arguments are the heterozygote count and the two homozygote counts; the
    routine self-corrects which homozygote is the rare one.
    """
    obs_hom_r = obs_hom1
    obs_hom_c = obs_hom2
    diplotypes = obs_hom_r + obs_het + obs_hom_c
    if diplotypes == 0:
        return 1.0
    rare = obs_hom_r * 2 + obs_het
    hets = obs_het
    if rare > diplotypes:
        rare = 2 * diplotypes - rare
    if hets > rare:
        # numerically impossible; treat as monomorphic-ish
        return 1.0

    tail_probs = [0.0] * (rare + 1)

    mid = int(rare * (2 * diplotypes - rare) / (2 * diplotypes))
    if ((rare & 1) ^ (mid & 1)) != 0:
        mid += 1

    het = mid
    hom_r = (rare - mid) // 2
    hom_c = diplotypes - het - hom_r

    tail_probs[mid] = 1.0
    total = tail_probs[mid]

    het = mid
    while het > 1:
        tail_probs[het - 2] = (tail_probs[het] * het * (het - 1.0)) / (
            4.0 * (hom_r + 1.0) * (hom_c + 1.0))
        total += tail_probs[het - 2]
        hom_r += 1
        hom_c += 1
        het -= 2

    het = mid
    hom_r = (rare - mid) // 2
    hom_c = diplotypes - het - hom_r
    while het <= rare - 2:
        tail_probs[het + 2] = (tail_probs[het] * 4.0 * hom_r * hom_c) / (
            (het + 2.0) * (het + 1.0))
        total += tail_probs[het + 2]
        hom_r -= 1
        hom_c -= 1
        het += 2

    if total <= 0:
        return 1.0
    tail_probs = [p / total for p in tail_probs]

    top = tail_probs[hets]
    for i in range(hets + 1, rare + 1):
        top += tail_probs[i]
    other = tail_probs[hets]
    for i in range(hets - 1, -1, -1):
        other += tail_probs[i]

    if top > 0.5 and other > 0.5:
        return 1.0
    return min(top, other) * 2.0


@dataclass
class QCParams:
    """Marker-exclusion thresholds (defaults match Haploview's CLI)."""

    min_geno_percent: float = 0.75   # -minGeno
    hwe_cutoff: float = 0.001        # -hwcutoff (0 disables)
    min_maf: float = 0.0             # -minMAF (0 disables)
    missing_threshold: float = 0.5   # -missingcutoff: drop indivs > this missing


@dataclass
class MarkerStatus:
    index: int
    name: str
    kept: bool
    reason: str


def filter_markers(dataset: Dataset, params: QCParams):
    """Return (kept_indices, statuses) applying Haploview-style QC.

    A marker is dropped if monomorphic, if its genotyping rate is below
    ``min_geno_percent``, if its HWE p-value is below ``hwe_cutoff``, or if its
    MAF is below ``min_maf``.
    """
    kept: List[int] = []
    statuses: List[MarkerStatus] = []
    for i, mk in enumerate(dataset.markers):
        geno_rate = 1.0 - mk.missing_fraction
        reason = "OK"
        keep = True
        if mk.maf <= 0.0:
            keep, reason = False, "monomorphic"
        elif geno_rate < params.min_geno_percent:
            keep, reason = False, "genopercent"
        elif params.min_maf > 0 and mk.maf < params.min_maf:
            keep, reason = False, "minMAF"
        elif params.hwe_cutoff > 0 and mk.hwe_p < params.hwe_cutoff:
            keep, reason = False, "HWE"
        statuses.append(MarkerStatus(i, mk.name, keep, reason))
        if keep:
            kept.append(i)
    return kept, statuses
