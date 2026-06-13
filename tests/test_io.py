"""Cross-format reading equivalence."""

import os

import numpy as np

from haploview.io import read_vcf, read_hapmap, read_linkage

DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def _maf_vector(ds):
    return [round(m.maf, 3) for m in ds.markers]


def test_formats_agree_on_genotypes_and_maf():
    vcf = read_vcf(os.path.join(DATA, "sample.vcf"))
    hmp = read_hapmap(os.path.join(DATA, "sample.hmp"))
    ped = read_linkage(os.path.join(DATA, "sample.ped"),
                       os.path.join(DATA, "sample.info"))

    assert vcf.n_markers == hmp.n_markers == ped.n_markers
    assert vcf.n_samples == hmp.n_samples == ped.n_samples
    # The dosage reference allele can differ (VCF keys to ALT, linkage to the
    # minor allele), but the canonical minor-allele dosage is invariant.
    for i in range(vcf.n_markers):
        v = vcf.minor_dosage_row(i)
        h = hmp.minor_dosage_row(i)
        p = ped.minor_dosage_row(i)
        assert np.array_equal(v, h)
        assert np.array_equal(v, p)
    assert _maf_vector(vcf) == _maf_vector(hmp) == _maf_vector(ped)


def test_vcf_marker_metadata():
    vcf = read_vcf(os.path.join(DATA, "sample.vcf"))
    assert vcf.markers[0].name == "rs1"
    assert vcf.markers[0].chrom == "1"
    assert vcf.markers[0].position == 1000
