"""Block-finding behaviour on constructed LD patterns."""

import numpy as np

from conftest import make_dataset
from haploview.blocks import find_blocks
from haploview.ld import LDTable


def _strong_block_dataset():
    # 4 SNPs in tight LD (two haplotypes), plus a 5th independent SNP far away
    rng = np.random.default_rng(7)
    n = 200
    base = rng.integers(0, 2, size=n)          # 0/1 haplotype indicator
    other = rng.integers(0, 2, size=n)
    def hom(ind):  # homozygous dosage from a haplotype indicator
        return (ind * 2).tolist()
    rows = [hom(base), hom(base), hom(base), hom(base), hom(other)]
    positions = [1000, 1200, 1400, 1600, 400000]
    return make_dataset(rows, positions=positions)


def test_gabriel_finds_block():
    ds = _strong_block_dataset()
    ld = LDTable(ds, max_dist_bp=500_000)
    blocks = find_blocks("gabriel", ds, ld)
    assert len(blocks) == 1
    assert blocks[0] == [0, 1, 2, 3]


def test_four_gamete_and_spine_run():
    ds = _strong_block_dataset()
    ld = LDTable(ds, max_dist_bp=500_000)
    for method in ("4gam", "spine"):
        blocks = find_blocks(method, ds, ld)
        # the four correlated SNPs should be grouped together
        assert any(set([0, 1, 2, 3]).issubset(set(b)) for b in blocks)


def test_recombinant_splits_four_gamete():
    # introduce all four gametes between two SNPs -> they cannot be one block
    rng = np.random.default_rng(3)
    n = 300
    a = rng.integers(0, 2, size=n)
    b = rng.integers(0, 2, size=n)  # independent -> 4 gametes present
    rows = [(a * 2).tolist(), (b * 2).tolist()]
    ds = make_dataset(rows, positions=[1000, 1200])
    ld = LDTable(ds, max_dist_bp=500_000)
    blocks = find_blocks("4gam", ds, ld)
    assert blocks == []
