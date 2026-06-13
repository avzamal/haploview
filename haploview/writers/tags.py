"""Write tag-SNP outputs: a subset VCF and a block<->tag<->SNP correspondence CSV."""

from __future__ import annotations

import csv
from typing import Dict, List

from ..model import Dataset, MISSING
from ..tagger import TagResult

VCF_VERSION = "VCFv4.2"


def _block_of(blocks: List[List[int]]) -> Dict[int, str]:
    out: Dict[int, str] = {}
    for i, block in enumerate(blocks):
        for m in block:
            out[m] = f"BLOCK{i + 1}"
    return out


def write_tag_vcf(dataset: Dataset, tags: List[int], path: str,
                  source: str = "haploview-py") -> None:
    """Write a standard biallelic VCF containing only the tag SNPs."""
    contigs = []
    seen = set()
    for t in tags:
        c = dataset.markers[t].chrom
        if c not in seen:
            seen.add(c)
            contigs.append(c)

    genos = dataset.genotypes
    gt = {0: "0/0", 1: "0/1", 2: "1/1", MISSING: "./."}
    with open(path, "w") as out:
        out.write(f"##fileformat={VCF_VERSION}\n")
        out.write(f"##source={source};content=tagSNPs\n")
        for c in contigs:
            out.write(f"##contig=<ID={c}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(dataset.samples) + "\n")
        for t in sorted(tags, key=lambda i: (dataset.markers[i].chrom,
                                             dataset.markers[i].position)):
            mk = dataset.markers[t]
            calls = "\t".join(gt[int(genos[t, j])] for j in range(dataset.n_samples))
            out.write(f"{mk.chrom}\t{mk.position}\t{mk.name}\t{mk.a1}\t{mk.a2}\t.\t"
                      f"PASS\tTAG\tGT\t{calls}\n")


def write_tag_csv(dataset: Dataset, result: TagResult, blocks: List[List[int]],
                  membership_path: str, summary_path: str) -> None:
    """Write the SNP->tag correspondence and the per-tag summary, both annotated
    with haploblock membership."""
    block_of = _block_of(blocks)
    tag_set = set(result.tags)

    with open(membership_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["snp_id", "chrom", "pos", "block_id", "is_tag",
                    "best_tag", "best_tag_pos", "best_tag_block",
                    "r2_with_tag", "captured"])
        for i, mk in enumerate(dataset.markers):
            bt = result.best_tag.get(i)
            captured = bt is not None
            w.writerow([
                mk.name, mk.chrom, mk.position, block_of.get(i, "-"),
                "yes" if i in tag_set else "no",
                dataset.markers[bt].name if captured else "-",
                dataset.markers[bt].position if captured else "-",
                block_of.get(bt, "-") if captured else "-",
                f"{result.best_rsq[i]:.3f}" if captured else "-",
                "yes" if captured else "no",
            ])

    with open(summary_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tag_id", "chrom", "pos", "block_id", "n_captured",
                    "captured_snp_ids", "captured_snp_blocks"])
        for t in result.tags:
            mk = dataset.markers[t]
            captured = result.tag_captures[t]
            ids = ";".join(dataset.markers[c].name for c in captured)
            blks = ";".join(sorted(set(block_of.get(c, "-") for c in captured)))
            w.writerow([mk.name, mk.chrom, mk.position, block_of.get(t, "-"),
                        len(captured), ids, blks])
