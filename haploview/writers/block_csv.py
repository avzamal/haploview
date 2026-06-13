"""CSV outputs describing block membership and block haplotypes."""

from __future__ import annotations

import csv
from typing import List

from ..model import Dataset


def write_block_csv(dataset: Dataset, blocks: List[List[int]], block_haplotypes,
                    method: str, membership_path: str, summary_path: str) -> None:
    """Write two CSVs: one (block, SNP) membership row per SNP, and a per-block
    haplotype summary."""
    with open(membership_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block_id", "method", "chrom", "block_start", "block_end",
                    "block_span_bp", "n_snps", "snp_index_in_block", "snp_id",
                    "snp_pos", "ref_allele", "alt_allele", "maf"])
        for i, block in enumerate(blocks):
            chrom = dataset.markers[block[0]].chrom
            start = dataset.markers[block[0]].position
            end = dataset.markers[block[-1]].position
            for j, m in enumerate(block):
                mk = dataset.markers[m]
                w.writerow([f"BLOCK{i + 1}", method, chrom, start, end,
                            end - start, len(block), j + 1, mk.name, mk.position,
                            mk.major_allele, mk.minor_allele, f"{mk.maf:.4f}"])

    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block_id", "method", "chrom", "block_start", "block_end",
                    "n_snps", "snp_ids", "haplotype", "frequency"])
        for i, (block, bh) in enumerate(zip(blocks, block_haplotypes)):
            chrom = dataset.markers[block[0]].chrom
            start = dataset.markers[block[0]].position
            end = dataset.markers[block[-1]].position
            snp_ids = ";".join(dataset.markers[m].name for m in block)
            order = sorted(range(len(bh.frequencies)),
                           key=lambda k: -bh.frequencies[k])
            for k in order:
                hap = bh.haplotypes[k]
                bases = "".join(
                    dataset.markers[block[j]].minor_allele if a == 1
                    else dataset.markers[block[j]].major_allele
                    for j, a in enumerate(hap))
                w.writerow([f"BLOCK{i + 1}", method, chrom, start, end,
                            len(block), snp_ids, bases, f"{bh.frequencies[k]:.4f}"])
