import numpy as np

from haploview.model import Dataset, Marker


def make_dataset(geno_rows, positions=None, names=None, alleles=None,
                 chrom="1"):
    """Build a Dataset from a list of per-marker dosage rows (count of a2).

    Each row is a list of 0/1/2/-1 over individuals.
    """
    n = len(geno_rows)
    positions = positions or [1000 * (i + 1) for i in range(n)]
    names = names or [f"rs{i+1}" for i in range(n)]
    alleles = alleles or [("A", "G")] * n
    markers = [Marker(name=names[i], chrom=chrom, position=positions[i],
                      a1=alleles[i][0], a2=alleles[i][1]) for i in range(n)]
    genos = np.array(geno_rows, dtype=np.int8)
    samples = [f"IND{j}" for j in range(genos.shape[1])]
    return Dataset(markers, genos, samples)
