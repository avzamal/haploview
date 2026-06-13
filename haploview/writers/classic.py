"""Classic Haploview-style text outputs: the LD table and the block listing."""

from __future__ import annotations

from typing import List

from ..ld import LDTable
from ..model import Dataset


def write_ld_table(dataset: Dataset, dprime: LDTable, path: str) -> None:
    """Write the pairwise LD table (one row per in-window marker pair).

    Columns match Haploview's ``.LD`` dump: ``L1 L2 D' LOD r^2 CIlow CIhi Dist
    T-int``. The T-int (multi-marker) column is emitted as ``-``; it is not used
    for parity, which compares D'/LOD/r^2/CI/Dist.
    """
    with open(path, "w") as out:
        out.write("L1\tL2\tD'\tLOD\tr^2\tCIlow\tCIhi\tDist\tT-int\n")
        n = dataset.n_markers
        for x in range(n - 1):
            for y in range(x + 1, n):
                pair = dprime.get(x, y)
                if pair is None:
                    continue
                dist = dataset.markers[y].position - dataset.markers[x].position
                out.write(
                    f"{dataset.markers[x].name}\t{dataset.markers[y].name}\t"
                    f"{_fmt(pair.dprime)}\t{_fmt(pair.lod)}\t{_fmt(pair.rsq)}\t"
                    f"{_fmt(pair.low_ci)}\t{_fmt(pair.high_ci)}\t{dist}\t-\n")


def _fmt(v: float) -> str:
    # compact like Java's Double.toString (trailing-zero trimmed, but >=1 dp)
    s = f"{v:.3f}".rstrip("0")
    if s.endswith("."):
        s += "0"
    return s


def write_blocks_text(dataset: Dataset, blocks: List[List[int]],
                      block_haplotypes, path: str, hap_thresh: float = 0.01) -> None:
    """Write the block listing in Haploview ``.GABRIELblocks`` style.

    Each block has a ``BLOCK n.  MARKERS: <1-based indices>`` header followed by
    its common haplotypes (allele-base strings) and population frequencies.
    """
    with open(path, "w") as out:
        for i, (block, bh) in enumerate(zip(blocks, block_haplotypes)):
            nums = " ".join(str(m + 1) for m in block)
            out.write(f"BLOCK {i + 1}.  MARKERS: {nums}\n")
            order = sorted(range(len(bh.frequencies)),
                           key=lambda k: -bh.frequencies[k])
            for k in order:
                if bh.frequencies[k] < hap_thresh:
                    continue
                hap = bh.haplotypes[k]
                bases = "".join(
                    dataset.markers[block[j]].minor_allele if a == 1
                    else dataset.markers[block[j]].major_allele
                    for j, a in enumerate(hap))
                out.write(f"{bases} ({bh.frequencies[k]:.3f})\n")
