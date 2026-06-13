"""Tag-SNP selection (pairwise Tagger) checks."""

import numpy as np

from conftest import make_dataset
from haploview.ld import LDTable
from haploview.tagger import select_tags


def test_perfect_ld_one_tag():
    # five SNPs in perfect LD -> a single tag captures all of them
    row = [0, 0, 1, 1, 2, 2, 0, 2, 1, 0] * 4
    ds = make_dataset([list(row) for _ in range(5)])
    ld = LDTable(ds, max_dist_bp=500_000)
    res = select_tags(ds, ld, rsq_cutoff=0.8)
    assert len(res.tags) == 1
    # every SNP captured at r^2 == 1
    assert len(res.best_tag) == ds.n_markers
    assert all(abs(r - 1.0) < 1e-9 for r in res.best_rsq.values())


def test_independent_snps_each_its_own_tag():
    rng = np.random.default_rng(11)
    rows = [rng.integers(0, 3, size=300).tolist() for _ in range(4)]
    ds = make_dataset(rows)
    ld = LDTable(ds, max_dist_bp=500_000)
    res = select_tags(ds, ld, rsq_cutoff=0.8)
    assert len(res.tags) == ds.n_markers


def test_every_snp_is_captured():
    # mix of a correlated cluster and a lone SNP
    base = ([0] * 30 + [2] * 30)
    rows = [list(base), list(base), [0, 1, 2] * 20]
    ds = make_dataset(rows)
    ld = LDTable(ds, max_dist_bp=500_000)
    res = select_tags(ds, ld, rsq_cutoff=0.8)
    assert res.untagged == []
    for snp, tag in res.best_tag.items():
        # a SNP's best tag really does capture it at >= cutoff
        from haploview.tagger import _rsq
        assert _rsq(ld, snp, tag) >= 0.8 - 1e-9
