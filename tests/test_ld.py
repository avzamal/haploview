"""LD math checks against hand-computable cases."""

import numpy as np

from conftest import make_dataset
from haploview.ld import compute_dprime


def test_perfect_ld():
    # two markers in perfect LD: genotypes identical -> D'=1, r^2=1
    row = [0, 0, 1, 1, 2, 2, 0, 2, 1, 0] * 5
    ds = make_dataset([row, list(row)])
    pair = compute_dprime(ds, 0, 1)
    assert pair.dprime == 1.0
    assert pair.rsq == 1.0
    assert pair.lod > 0


def test_independent_markers():
    # independent markers -> low r^2
    rng = np.random.default_rng(0)
    a = rng.integers(0, 3, size=400).tolist()
    b = rng.integers(0, 3, size=400).tolist()
    ds = make_dataset([a, b])
    pair = compute_dprime(ds, 0, 1)
    assert pair.rsq < 0.1


def test_monomorphic_returns_none():
    ds = make_dataset([[0] * 20, [0, 1, 2] * 6 + [0, 1]])
    assert compute_dprime(ds, 0, 1) is None


def test_dprime_in_range_and_ci_ordered():
    rng = np.random.default_rng(1)
    a = rng.integers(0, 3, size=200).tolist()
    # b correlated with a but noisy
    b = [min(2, max(0, x + rng.integers(-1, 2))) for x in a]
    ds = make_dataset([a, b])
    p = compute_dprime(ds, 0, 1)
    assert 0.0 <= p.dprime <= 1.0
    assert 0.0 <= p.low_ci <= p.high_ci <= 1.0
