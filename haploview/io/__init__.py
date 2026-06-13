"""Genotype input readers and format conversion."""

from .vcf_reader import read_vcf
from .hapmap_reader import read_hapmap
from .linkage_reader import read_linkage
from .detect import detect_format, read_any

__all__ = [
    "read_vcf",
    "read_hapmap",
    "read_linkage",
    "detect_format",
    "read_any",
]
