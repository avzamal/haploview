"""Read genotypes from a VCF file (plain or gzipped).

Only biallelic SNPs are loaded; indels, multiallelic sites and symbolic alleles
are skipped (Haploview operates on biallelic markers). Genotypes are taken from
the ``GT`` subfield; phasing (``|`` vs ``/``) is ignored because the LD/EM math
treats input as unphased, exactly like the original tool.
"""

from __future__ import annotations

import gzip
from typing import List, Optional

import numpy as np

from ..model import Dataset, Marker, MISSING

_BASES = {"A", "C", "G", "T"}


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _parse_gt(field: str, gt_index: int) -> int:
    """Return ALT-allele dosage (0/1/2) or MISSING from a sample column."""
    if field == "." or field == "":
        return MISSING
    parts = field.split(":")
    if gt_index >= len(parts):
        return MISSING
    gt = parts[gt_index]
    if gt in (".", "./.", ".|.", ""):
        return MISSING
    alleles = gt.replace("|", "/").split("/")
    dosage = 0
    for a in alleles:
        if a == ".":
            return MISSING
        if a == "0":
            continue
        elif a == "1":
            dosage += 1
        else:
            # allele index >1 shouldn't occur for a biallelic record
            return MISSING
    return dosage


def read_vcf(path: str, chrom: Optional[str] = None,
             max_markers: Optional[int] = None) -> Dataset:
    samples: List[str] = []
    markers: List[Marker] = []
    rows: List[np.ndarray] = []

    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                continue
            if not samples:
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            c, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            if chrom is not None and c != chrom:
                # markers are position-sorted within a contig; if we've already
                # collected this contig, stop once we pass it
                if markers and markers[-1].chrom == chrom:
                    break
                continue
            if ref not in _BASES or alt not in _BASES:
                continue  # not a simple biallelic SNP
            fmt = cols[8].split(":")
            try:
                gt_index = fmt.index("GT")
            except ValueError:
                continue
            dosage = np.fromiter(
                (_parse_gt(g, gt_index) for g in cols[9:]),
                dtype=np.int8, count=len(samples))
            name = vid if vid not in (".", "") else f"{c}:{pos}"
            markers.append(Marker(name=name, chrom=c, position=int(pos),
                                  a1=ref, a2=alt))
            rows.append(dosage)
            if max_markers is not None and len(markers) >= max_markers:
                break

    if not markers:
        raise ValueError(f"no biallelic SNPs read from {path}")
    genotypes = np.vstack(rows)
    return Dataset(markers, genotypes, samples)
