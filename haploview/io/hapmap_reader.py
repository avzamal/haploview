"""Read genotypes from a HapMap genotype file (.hmp).

HapMap layout: 11 leading metadata columns
``rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode``
followed by one column per sample, each a two-character genotype (e.g. ``AG``;
``NN`` for missing).
"""

from __future__ import annotations

import gzip
from typing import List

import numpy as np

from ..model import Dataset, Marker, MISSING

_META_COLS = 11


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_hapmap(path: str) -> Dataset:
    samples: List[str] = []
    markers: List[Marker] = []
    rows: List[np.ndarray] = []

    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            toks = line.split()
            if toks[0] in ("rs#", "rs.", "rsID") or toks[0].lower().startswith("rs#"):
                samples = toks[_META_COLS:]
                continue
            if not samples:
                # tolerate header without rs# marker
                if toks[0].lower() in ("strand", "chrom"):
                    continue
            name = toks[0]
            alleles = toks[1]
            if "/" in alleles:
                a1, a2 = alleles.split("/")[:2]
            else:
                a1, a2 = alleles[0], alleles[-1]
            chrom = toks[2]
            pos = int(toks[3])
            geno_cols = toks[_META_COLS:]
            dosage = np.empty(len(geno_cols), dtype=np.int8)
            for i, g in enumerate(geno_cols):
                g = g.upper()
                if "N" in g or "-" in g or len(g) < 2:
                    dosage[i] = MISSING
                else:
                    dosage[i] = (1 if g[0] == a2 else 0) + (1 if g[1] == a2 else 0)
            markers.append(Marker(name=name, chrom=chrom, position=pos,
                                  a1=a1, a2=a2))
            rows.append(dosage)
            if not samples:
                samples = [f"s{i+1}" for i in range(len(geno_cols))]

    if not markers:
        raise ValueError(f"no markers read from {path}")
    genotypes = np.vstack(rows)
    return Dataset(markers, genotypes, samples)
