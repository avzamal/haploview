"""Format auto-detection and a single dispatch entry point."""

from __future__ import annotations

import gzip
from typing import Optional

from ..model import Dataset


def _peek(path: str) -> str:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.strip():
                return line
    return ""


def detect_format(path: str) -> str:
    low = path.lower()
    if low.endswith((".vcf", ".vcf.gz")):
        return "vcf"
    if low.endswith((".hmp", ".hmp.txt", ".hapmap")):
        return "hapmap"
    if low.endswith(".ped"):
        return "linkage"
    first = _peek(path)
    if first.startswith("##fileformat=VCF") or first.startswith("#CHROM"):
        return "vcf"
    if first.split()[:1] == ["rs#"] or first.lower().startswith("rs#"):
        return "hapmap"
    return "linkage"


def read_any(path: str, fmt: Optional[str] = None,
             info: Optional[str] = None, chrom: Optional[str] = None,
             max_markers: Optional[int] = None) -> Dataset:
    if fmt in (None, "auto"):
        fmt = detect_format(path)
    if fmt == "vcf":
        from .vcf_reader import read_vcf
        return read_vcf(path, chrom=chrom, max_markers=max_markers)
    if fmt == "hapmap":
        from .hapmap_reader import read_hapmap
        return read_hapmap(path)
    if fmt == "linkage":
        from .linkage_reader import read_linkage
        return read_linkage(path, info_path=info)
    raise ValueError(f"unknown format: {fmt!r}")
