"""Within-block haplotype EM checks."""

from conftest import make_dataset
from haploview.haplotype_em import estimate_block_haplotypes, _block_sizes


def test_block_sizes_partition():
    assert _block_sizes(4) == [4]
    assert _block_sizes(8) == [8]
    assert _block_sizes(16) == [8, 8]
    # >8 with remainder splits the last segment
    assert sum(_block_sizes(10)) == 10
    assert all(s <= 8 for s in _block_sizes(20))


def test_two_haplotype_block_frequencies():
    # 3 SNPs forming two haplotypes 000 and 111 at 60/40 from homozygous indivs
    # 6 individuals 000/000, 4 individuals 111/111 -> freqs 0.6 / 0.4
    rows = [[0] * 6 + [2] * 4, [0] * 6 + [2] * 4, [0] * 6 + [2] * 4]
    ds = make_dataset(rows)
    bh = estimate_block_haplotypes(ds, [0, 1, 2])
    freqs = sorted(bh.frequencies, reverse=True)
    assert abs(freqs[0] - 0.6) < 0.02
    assert abs(freqs[1] - 0.4) < 0.02


def test_phase_resolution_with_double_het():
    # individuals all 0/0 except some double hets; EM should still converge
    rows = [[0, 0, 1, 1, 2, 0, 1, 0, 2, 1],
            [0, 0, 1, 1, 2, 0, 1, 0, 2, 1]]
    ds = make_dataset(rows)
    bh = estimate_block_haplotypes(ds, [0, 1])
    assert abs(sum(bh.frequencies) - 1.0) < 1e-6
