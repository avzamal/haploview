"""End-to-end pipeline / output checks."""

import csv
import os

from haploview.pipeline import run_analysis, write_outputs
from haploview.qc import QCParams

DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def test_end_to_end_outputs(tmp_path):
    result = run_analysis(os.path.join(DATA, "sample.vcf"), method="gabriel")
    assert len(result.blocks) == 2

    prefix = str(tmp_path / "out")
    paths = write_outputs(result, prefix, "gabriel")
    for p in paths:
        assert os.path.exists(p)

    # block VCF: one record per block, multiallelic where expected
    vcf_lines = [l for l in open(prefix + ".blocks.vcf") if not l.startswith("#")]
    assert len(vcf_lines) == 2
    for line in vcf_lines:
        cols = line.split("\t")
        ref = cols[3]
        alts = cols[4].split(",")
        # REF and ALT alleles are distinct haplotype strings of equal length
        for a in alts:
            assert a == "." or len(a) == len(ref)

    # CSV membership row count equals total SNPs across blocks
    with open(prefix + ".blocks.csv") as fh:
        rows = list(csv.DictReader(fh))
    total_snps = sum(len(b) for b in result.blocks)
    assert len(rows) == total_snps


def test_block_haplotype_freqs_sum_to_one():
    result = run_analysis(os.path.join(DATA, "sample.vcf"), method="gabriel")
    for bh in result.block_haplotypes:
        assert abs(sum(bh.frequencies) - 1.0) < 0.05


def test_methods_all_run():
    for method in ("gabriel", "4gam", "spine"):
        result = run_analysis(os.path.join(DATA, "sample.vcf"), method=method)
        assert isinstance(result.blocks, list)
