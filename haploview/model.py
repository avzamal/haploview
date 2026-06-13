"""Core data model: markers and a genotype matrix.

Genotypes are stored as an ``int8`` matrix of shape ``(n_markers, n_samples)``.
Each entry is the dosage (0, 1, 2) of a marker's ``a2`` allele, or ``-1`` for
missing data. This is enough to reconstruct everything the LD/EM math needs:
het iff dosage == 1, the two homozygotes are dosage 0 and 2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

MISSING = -1


def _round3(x: float) -> float:
    """Java Math.round half-up to 3 decimals (Haploview's Util.roundDouble)."""
    return math.floor(x * 1000 + 0.5) / 1000


@dataclass
class Marker:
    """A single biallelic marker.

    ``a1``/``a2`` are the allele characters; the genotype dosage counts copies
    of ``a2``. ``major``/``minor`` are assigned from observed frequencies, the
    way Haploview does (major = more frequent allele).
    """

    name: str
    chrom: str
    position: int
    a1: str
    a2: str
    # Filled in by Dataset.compute_stats()
    maf: float = 0.0
    missing_fraction: float = 0.0
    hwe_p: float = 1.0
    minor_is_a2: bool = True
    # genotype counts of (a2a2 hom, het, a1a1 hom) among non-missing
    counts: tuple = (0, 0, 0)

    @property
    def major_allele(self) -> str:
        return self.a2 if not self.minor_is_a2 else self.a1

    @property
    def minor_allele(self) -> str:
        return self.a2 if self.minor_is_a2 else self.a1


class Dataset:
    """Markers plus their genotype matrix and sample identifiers."""

    def __init__(self, markers: List[Marker], genotypes: np.ndarray,
                 samples: List[str]):
        if genotypes.shape[0] != len(markers):
            raise ValueError("genotype rows must equal number of markers")
        if genotypes.shape[1] != len(samples):
            raise ValueError("genotype cols must equal number of samples")
        self.markers = markers
        self.genotypes = genotypes.astype(np.int8, copy=False)
        self.samples = samples
        self.compute_stats()

    @property
    def n_markers(self) -> int:
        return len(self.markers)

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    def compute_stats(self) -> None:
        """Populate per-marker MAF, missingness, HWE p-value and major/minor."""
        from .qc import hwe_exact_p

        for i, mk in enumerate(self.markers):
            row = self.genotypes[i]
            nonmiss = row[row != MISSING]
            n = nonmiss.size
            n_hom_a2 = int(np.count_nonzero(nonmiss == 2))
            n_het = int(np.count_nonzero(nonmiss == 1))
            n_hom_a1 = int(np.count_nonzero(nonmiss == 0))
            mk.counts = (n_hom_a2, n_het, n_hom_a1)
            total_alleles = 2 * n
            if total_alleles == 0:
                mk.maf = 0.0
                mk.missing_fraction = 1.0
                mk.minor_is_a2 = True
                mk.hwe_p = 1.0
                continue
            count_a2 = 2 * n_hom_a2 + n_het
            p_a2 = count_a2 / total_alleles
            p_a1 = 1.0 - p_a2
            # Haploview stores the analysis MAF rounded to 3 decimals; the block
            # MAF threshold (0.05) is compared against this rounded value.
            mk.maf = _round3(min(p_a1, p_a2))
            mk.minor_is_a2 = p_a2 <= p_a1
            mk.missing_fraction = float(
                1.0 - n / self.n_samples) if self.n_samples else 0.0
            # HWE exact test on the rarer-allele homozygote convention
            n_minor_hom = n_hom_a2 if mk.minor_is_a2 else n_hom_a1
            n_major_hom = n_hom_a1 if mk.minor_is_a2 else n_hom_a2
            mk.hwe_p = hwe_exact_p(n_het, n_minor_hom, n_major_hom)

    def minor_dosage_row(self, i: int) -> np.ndarray:
        """Return marker ``i`` recoded as minor-allele dosage (missing = -1).

        0 = homozygous major, 1 = het, 2 = homozygous minor.
        """
        row = self.genotypes[i]
        mk = self.markers[i]
        if mk.minor_is_a2:
            return row
        out = row.copy()
        mask = out != MISSING
        out[mask] = 2 - out[mask]
        return out

    def subset(self, indices: List[int]) -> "Dataset":
        markers = [self.markers[i] for i in indices]
        genos = self.genotypes[indices, :]
        return Dataset(markers, genos, list(self.samples))

    def drop_high_missing_samples(self, missing_threshold: float = 0.5) -> "Dataset":
        """Drop individuals genotyped at fewer than (1 - threshold) of markers.

        Matches Haploview's per-individual filter (PedFile.check): an individual
        with ``genoPC < 1 - missingThreshold`` is removed before marker stats are
        computed. Returns a new Dataset (or self if nothing is dropped).
        """
        if self.n_markers == 0:
            return self
        present = (self.genotypes != MISSING).mean(axis=0)
        keep = [j for j in range(self.n_samples)
                if present[j] >= (1.0 - missing_threshold)]
        if len(keep) == self.n_samples:
            return self
        genos = self.genotypes[:, keep]
        samples = [self.samples[j] for j in keep]
        # rebuild markers fresh so per-marker stats are recomputed on survivors
        markers = [Marker(name=m.name, chrom=m.chrom, position=m.position,
                          a1=m.a1, a2=m.a2) for m in self.markers]
        return Dataset(markers, genos, samples)
