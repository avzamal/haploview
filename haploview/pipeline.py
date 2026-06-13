"""End-to-end analysis pipeline shared by the CLI and tests."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from .blocks import find_blocks
from .haplotype_em import BlockHaplotypes, estimate_block_haplotypes
from .io import read_any
from .ld import LDTable
from .model import Dataset
from .qc import QCParams, filter_markers
from .tagger import TagResult, select_tags
from .writers import (write_block_csv, write_block_vcf, write_blocks_text,
                      write_ld_table, write_tag_csv, write_tag_vcf)


@dataclass
class AnalysisResult:
    dataset: Dataset                       # post-QC dataset actually analysed
    ld: LDTable
    blocks: List[List[int]]
    block_haplotypes: List[BlockHaplotypes]
    dropped: list = field(default_factory=list)
    tags: Optional[TagResult] = None


def run_analysis(input_path: str, fmt: Optional[str] = None,
                 info: Optional[str] = None, method: str = "gabriel",
                 max_distance_kb: int = 500, qc: Optional[QCParams] = None,
                 chrom: Optional[str] = None, max_markers: Optional[int] = None,
                 tag: bool = True, tag_rsq: float = 0.8) -> AnalysisResult:
    qc = qc or QCParams()
    raw = read_any(input_path, fmt=fmt, info=info, chrom=chrom,
                   max_markers=max_markers)
    # drop individuals with too much missing data before computing marker stats
    raw = raw.drop_high_missing_samples(qc.missing_threshold)
    kept, statuses = filter_markers(raw, qc)
    dataset = raw.subset(kept)
    ld = LDTable(dataset, max_dist_bp=max_distance_kb * 1000)
    blocks = find_blocks(method, dataset, ld)
    block_haps = [estimate_block_haplotypes(dataset, b) for b in blocks]
    dropped = [s for s in statuses if not s.kept]
    tags = select_tags(dataset, ld, rsq_cutoff=tag_rsq) if tag else None
    return AnalysisResult(dataset, ld, blocks, block_haps, dropped, tags)


def write_outputs(result: AnalysisResult, prefix: str, method: str,
                  hap_thresh: float = 0.01) -> List[str]:
    paths = []
    ld_path = f"{prefix}.LD"
    write_ld_table(result.dataset, result.ld, ld_path)
    paths.append(ld_path)

    blocks_path = f"{prefix}.blocks"
    write_blocks_text(result.dataset, result.blocks, result.block_haplotypes,
                      blocks_path, hap_thresh=hap_thresh)
    paths.append(blocks_path)

    membership = f"{prefix}.blocks.csv"
    summary = f"{prefix}.block_haplotypes.csv"
    write_block_csv(result.dataset, result.blocks, result.block_haplotypes,
                    method, membership, summary)
    paths += [membership, summary]

    vcf_path = f"{prefix}.blocks.vcf"
    write_block_vcf(result.dataset, result.blocks, result.block_haplotypes,
                    method, vcf_path)
    paths.append(vcf_path)

    if result.tags is not None:
        tag_vcf = f"{prefix}.tagsnps.vcf"
        write_tag_vcf(result.dataset, result.tags.tags, tag_vcf)
        tag_csv = f"{prefix}.tags.csv"
        tag_summary = f"{prefix}.tags_summary.csv"
        write_tag_csv(result.dataset, result.tags, result.blocks,
                      tag_csv, tag_summary)
        paths += [tag_vcf, tag_csv, tag_summary]
    return paths
