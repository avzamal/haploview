"""Convert a loaded :class:`Dataset` into legacy Haploview input formats.

Used both as a user-facing ``haploview convert`` subcommand and to drive the
parity test, where the original Java Haploview must read the *same* genotypes
that the Python tool read from VCF.
"""

from __future__ import annotations

from ..model import Dataset, MISSING


def dataset_to_linkage(dataset: Dataset, ped_path: str, info_path: str) -> None:
    """Write a Dataset as linkage ``.ped`` + ``.info`` (unrelated founders)."""
    with open(info_path, "w") as info:
        for mk in dataset.markers:
            info.write(f"{mk.name}\t{mk.position}\n")

    genos = dataset.genotypes
    with open(ped_path, "w") as ped:
        for j, sample in enumerate(dataset.samples):
            fields = [sample, sample, "0", "0", "0", "0"]
            for i, mk in enumerate(dataset.markers):
                d = int(genos[i, j])
                if d == MISSING:
                    fields += ["0", "0"]
                elif d == 0:
                    fields += [mk.a1, mk.a1]
                elif d == 1:
                    fields += [mk.a1, mk.a2]
                else:
                    fields += [mk.a2, mk.a2]
            ped.write(" ".join(fields) + "\n")


def dataset_to_hapmap(dataset: Dataset, path: str) -> None:
    """Write a Dataset in HapMap genotype format."""
    header = ["rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center",
              "protLSID", "assayLSID", "panelLSID", "QCcode"] + list(dataset.samples)
    genos = dataset.genotypes
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for i, mk in enumerate(dataset.markers):
            meta = [mk.name, f"{mk.a1}/{mk.a2}", mk.chrom, str(mk.position),
                    "+", "NA", "NA", "NA", "NA", "NA", "NA"]
            calls = []
            for j in range(dataset.n_samples):
                d = int(genos[i, j])
                if d == MISSING:
                    calls.append("NN")
                elif d == 0:
                    calls.append(mk.a1 + mk.a1)
                elif d == 1:
                    calls.append(mk.a1 + mk.a2)
                else:
                    calls.append(mk.a2 + mk.a2)
            out.write("\t".join(meta + calls) + "\n")
