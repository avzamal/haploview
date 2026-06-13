"""Read genotypes from a classic Haploview linkage ``.ped`` + ``.info`` pair.

``.ped``: family, individual, father, mother, sex, affection, then two allele
columns per marker. Alleles are ``ACGT`` or ``1-4`` (1=A,2=C,3=G,4=T); ``0`` is
missing. ``.info``: marker name and base-pair position, one marker per line, in
the same order as the genotype columns.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np

from ..model import Dataset, Marker, MISSING

_NUM2BASE = {"1": "A", "2": "C", "3": "G", "4": "T"}


def _read_info(info_path: str):
    names, positions = [], []
    with open(info_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            toks = line.split()
            names.append(toks[0])
            positions.append(int(toks[1]) if len(toks) > 1 else len(positions))
    return names, positions


def read_linkage(ped_path: str, info_path: Optional[str] = None) -> Dataset:
    if info_path is None:
        raise ValueError("linkage format requires an .info file (--info)")
    names, positions = _read_info(info_path)
    n_markers = len(names)

    samples: List[str] = []
    # allele tokens per marker per individual: list over inds of list[(t1,t2)]
    raw: List[List[tuple]] = []
    with open(ped_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 6 + 2 * n_markers:
                continue
            fam, iid = toks[0], toks[1]
            samples.append(iid if iid not in samples else f"{fam}:{iid}")
            alleles = toks[6:6 + 2 * n_markers]
            pairs = [(alleles[2 * m], alleles[2 * m + 1]) for m in range(n_markers)]
            raw.append(pairs)

    n_ind = len(samples)
    # determine a1/a2 per marker from observed allele tokens
    markers: List[Marker] = []
    genotypes = np.full((n_markers, n_ind), MISSING, dtype=np.int8)
    for m in range(n_markers):
        counts = {}
        for ind in range(n_ind):
            for t in raw[ind][m]:
                if t != "0":
                    counts[t] = counts.get(t, 0) + 1
        observed = sorted(counts, key=lambda k: (-counts[k], k))
        a1 = observed[0] if observed else "0"
        a2 = observed[1] if len(observed) > 1 else (a1 if a1 != "0" else "0")
        numeric = all(t in _NUM2BASE for t in (a1, a2) if t != "0")
        b1 = _NUM2BASE.get(a1, a1) if numeric else a1
        b2 = _NUM2BASE.get(a2, a2) if numeric else a2
        markers.append(Marker(name=names[m], chrom="0", position=positions[m],
                              a1=b1, a2=b2))
        for ind in range(n_ind):
            t1, t2 = raw[ind][m]
            if t1 == "0" or t2 == "0":
                continue
            genotypes[m, ind] = (1 if t1 == a2 else 0) + (1 if t2 == a2 else 0)

    return Dataset(markers, genotypes, samples)
