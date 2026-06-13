"""haploview: a modern reimplementation of the Haploview LD / block tool."""

__version__ = "0.1.0"

from .model import Dataset, Marker
from .ld import LDTable, compute_dprime, PairwiseLinkage
from .blocks import find_blocks
from .haplotype_em import estimate_block_haplotypes

__all__ = [
    "Dataset",
    "Marker",
    "LDTable",
    "compute_dprime",
    "PairwiseLinkage",
    "find_blocks",
    "estimate_block_haplotypes",
    "__version__",
]
