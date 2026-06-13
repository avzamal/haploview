"""Write haplotype blocks as a VCF where each block is one multiallelic record.

This is a new, machine-readable representation of Haploview's blocks: every
block becomes a single VCF record whose REF/ALT alleles are the block's common
haplotypes (REF = most frequent), with member SNPs and haplotype frequencies in
INFO and each sample genotyped by its most likely block diplotype.

The allele strings are the concatenated per-SNP bases of a haplotype. They are a
compact synthetic encoding of the haplotype, not a literal reference span.
"""

from __future__ import annotations

from typing import List

from ..model import Dataset

VCF_VERSION = "VCFv4.2"


def _hap_bases(dataset: Dataset, block: List[int], hap: List[int]) -> str:
    return "".join(
        dataset.markers[block[j]].minor_allele if a == 1
        else dataset.markers[block[j]].major_allele
        for j, a in enumerate(hap))


def write_block_vcf(dataset: Dataset, blocks: List[List[int]], block_haplotypes,
                    method: str, path: str, source: str = "haploview-py") -> None:
    contigs = []
    seen = set()
    for block in blocks:
        c = dataset.markers[block[0]].chrom
        if c not in seen:
            seen.add(c)
            contigs.append(c)

    with open(path, "w") as out:
        out.write(f"##fileformat={VCF_VERSION}\n")
        out.write(f"##source={source};method={method}\n")
        out.write('##ALT=<ID=HAP,Description="A haplotype block allele '
                  '(concatenated SNP bases)">\n')
        for c in contigs:
            out.write(f"##contig=<ID={c}>\n")
        out.write('##INFO=<ID=BLOCKID,Number=1,Type=String,Description="Block identifier">\n')
        out.write('##INFO=<ID=METHOD,Number=1,Type=String,Description="Block-definition method">\n')
        out.write('##INFO=<ID=NSNP,Number=1,Type=Integer,Description="Number of SNPs in block">\n')
        out.write('##INFO=<ID=START,Number=1,Type=Integer,Description="Block start position">\n')
        out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Block end position">\n')
        out.write('##INFO=<ID=SNPIDS,Number=.,Type=String,Description="Member SNP identifiers">\n')
        out.write('##INFO=<ID=SNPPOS,Number=.,Type=Integer,Description="Member SNP positions">\n')
        out.write('##INFO=<ID=HAPS,Number=.,Type=String,Description="Common haplotype allele strings (REF first)">\n')
        out.write('##INFO=<ID=HAPFREQ,Number=.,Type=Float,Description="Haplotype frequencies aligned with HAPS">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Most likely block diplotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(dataset.samples) + "\n")

        for i, (block, bh) in enumerate(zip(blocks, block_haplotypes)):
            chrom = dataset.markers[block[0]].chrom
            start = dataset.markers[block[0]].position
            end = dataset.markers[block[-1]].position
            # order haplotypes by frequency descending; REF = most common
            order = sorted(range(len(bh.haplotypes)), key=lambda k: -bh.frequencies[k])
            remap = {old: new for new, old in enumerate(order)}
            allele_strings = [_hap_bases(dataset, block, bh.haplotypes[k]) for k in order]
            freqs_ordered = [bh.frequencies[k] for k in order]
            ref = allele_strings[0]
            alts = allele_strings[1:]
            alt_field = ",".join(alts) if alts else "."

            snp_ids = ",".join(dataset.markers[m].name for m in block)
            snp_pos = ",".join(str(dataset.markers[m].position) for m in block)
            info = (f"BLOCKID=BLOCK{i + 1};METHOD={method};NSNP={len(block)};"
                    f"START={start};END={end};SNPIDS={snp_ids};SNPPOS={snp_pos};"
                    f"HAPS={','.join(allele_strings)};"
                    f"HAPFREQ={','.join(f'{f:.4f}' for f in freqs_ordered)}")

            gts = []
            for j in range(dataset.n_samples):
                dip = bh.sample_diplotypes[j]
                if dip is None or dip[0] not in remap or dip[1] not in remap:
                    gts.append("./.")
                else:
                    a, b = sorted((remap[dip[0]], remap[dip[1]]))
                    gts.append(f"{a}|{b}")

            out.write(f"{chrom}\t{start}\tBLOCK{i + 1}\t{ref}\t{alt_field}\t.\t"
                      f"PASS\t{info}\tGT\t" + "\t".join(gts) + "\n")
