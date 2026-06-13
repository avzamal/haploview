# haploview (modern reimplementation)

A faithful, modern reimplementation of the classic
[Haploview](https://www.broadinstitute.org/haploview) linkage-disequilibrium and
haplotype-block analysis tool (Barrett *et al.*, *Bioinformatics* 2005), with:

- **modern genotype inputs**: VCF (plain or `.gz`) and HapMap, plus the classic
  linkage `.ped`/`.info` format;
- the **classic outputs** (pairwise D′/r²/LOD/CI table, block listing, per-block
  haplotypes with population frequencies);
- two **new machine-readable outputs**:
  1. a **VCF where each haploblock is a single multiallelic record** whose
     REF/ALT alleles are the block's common haplotypes, with member SNPs and
     haplotype frequencies in `INFO` and each sample genotyped by its most
     likely block diplotype;
  2. a **CSV describing which SNPs belong to each block** (plus a per-block
     haplotype-frequency summary CSV).

The analytical core is transcribed directly from the original Haploview Java
source (`edu/mit/wi/haploview/HaploData.java`, `FindBlocks.java`, `EM.java` and
`edu/mit/wi/pedfile/CheckData.java`), so results match the original tool — see
[Parity](#parity-with-the-original) below.

## Install

```bash
pip install -e .          # installs the `haploview` command (needs numpy)
```

## Usage

```bash
# Analyse a VCF with the default Gabriel block method
haploview analyze -i genotypes.vcf -o results

# Choose a method and thresholds
haploview analyze -i genotypes.vcf.gz -m 4gam --max-distance 200 -o results
haploview analyze -i genotypes.hmp -m all -o results        # all three methods

# Linkage input needs the marker .info file
haploview analyze -i data.ped --info data.info -o results

# Convert genotypes to legacy formats
haploview convert -i genotypes.vcf --to linkage -o data
haploview convert -i genotypes.vcf --to hapmap  -o data
```

Outputs written under the `-o` prefix:

| file | contents |
|------|----------|
| `<prefix>.LD` | pairwise table `L1 L2 D' LOD r² CIlow CIhi Dist T-int` |
| `<prefix>.blocks` | block listing with common haplotypes and frequencies |
| `<prefix>.blocks.vcf` | **new** multiallelic VCF, one record per haploblock |
| `<prefix>.blocks.csv` | **new** one row per (block, SNP) membership |
| `<prefix>.block_haplotypes.csv` | per-block haplotype frequencies |
| `<prefix>.tagsnps.vcf` | **new** input VCF subset to just the tag SNPs |
| `<prefix>.tags.csv` | **new** per-SNP: its block, best tag, r², and whether it is a tag |
| `<prefix>.tags_summary.csv` | **new** per-tag: the SNPs (and blocks) it captures |

### Tag SNPs

A **tag SNP** is a marker chosen so that genotyping it captures other markers it
is in strong LD with: Haploview's Tagger greedily picks a minimal set of SNPs
such that every SNP is either a tag or has **r² ≥ 0.8** (default, within 500 kb)
with a chosen tag (de Bakker *et al.*, *Nat Genet* 2005). Each captured SNP is
assigned its *best tag* (the chosen tag it correlates with most strongly).

This tool runs the pairwise Tagger by default (disable with `--no-tag`, change
the cutoff with `--tag-rsq`) and writes the tag-SNP subset VCF plus the
correspondence CSVs, which annotate every tag/captured SNP with its haploblock.

The greedy selection is a transcription of `edu/mit/wi/tagger/Tagger.java`
(pairwise mode). Because the original breaks ties among equally good candidate
tags using hash-map order (not reproducible), this tool breaks them
deterministically; the resulting tag set is an equivalent valid cover. On the
real 300-marker subset below, 112 of 117 capture groups are identical to the
original, both capture 100% of SNPs, and mean r² matches (0.96 vs 0.956).

### Block methods and key thresholds

All three Haploview methods are implemented with the original constants:

- **`gabriel`** – Gabriel *et al.* (2002) confidence-interval method (default):
  strong-LD pairs require `CIhigh ≥ 0.98` and `CIlow ≥ 0.70`; a block needs
  >95% of informative pairs in strong LD (`recHighCI = 0.90`), with the
  small-block CI variants and distance caps (20/30 kb) of the original.
- **`4gam`** – four-gamete test (a pair is compatible if ≤3 of its four
  haplotypes exceed frequency `0.01`).
- **`spine`** – solid spine of LD (`D′ ≥ 0.80`, tolerating one bad marker).

### Marker / individual QC (Haploview defaults)

- individuals genotyped at < 50% of markers are dropped (`--`, `missingThreshold`);
- markers below 75% genotyping (`--min-geno`), failing Hardy-Weinberg at
  p < 0.001 (`--hwe-cutoff`, exact test), or below `--min-maf` are excluded;
- the analysis MAF is rounded to 3 decimals before the 0.05 block threshold,
  matching the original.

## Parity with the original

The repository ships a harness that runs the **original Haploview 4.1 jar** and
this tool on the *same* genotypes and diffs their output:

```bash
bash scripts/get_haploview.sh                 # download the original jar
bash scripts/fetch_sample_vcf.sh              # carve a subset from a real VCF
HWE=0 GENO=0.5 bash scripts/validate_against_haploview.sh
```

On a 300-marker / 601-sample subset of a real GATK VCF
([Springer 12864-2021-07768 additional file 24](https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-021-07768-y/MediaObjects/12864_2021_7768_MOESM24_ESM.vcf))
this produces, against the original tool:

```
LD pairs: 1792 compared, 0 only-haploview, 0 only-mine, 0 value mismatches
Blocks: haploview=38 mine=38 IDENTICAL
PARITY: PASS -- LD table and blocks are identical.
```

`tests/test_against_haploview.py` runs the same check on the committed sample
data as part of `pytest` (skipped automatically if Java or the jar are absent).

## Tests

```bash
pytest -q
```

## References

- Barrett JC, Fry B, Maller J, Daly MJ. *Haploview: analysis and visualization
  of LD and haplotype maps.* Bioinformatics. 2005;21(2):263–5.
- Gabriel SB *et al.* *The structure of haplotype blocks in the human genome.*
  Science. 2002;296(5576):2225–9.
- Wigginton JE, Cutler DJ, Abecasis GR. *A note on exact tests of
  Hardy-Weinberg equilibrium.* Am J Hum Genet. 2005;76(5):887–93.
- Original source: github.com/BroadInstituteSoftware/Haploview,
  github.com/jazzywhit/Haploview (`edu/mit/wi/haploview/`).
