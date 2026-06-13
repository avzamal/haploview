#!/usr/bin/env bash
# End-to-end parity check: run the original Haploview jar and this tool on the
# same genotypes (a subset of the Springer VCF) and diff their LD table + blocks.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

VCF="${VCF:-data/chr01_subset.vcf}"
JAR="${JAR:-third_party/Haploview4.1.jar}"
OUTDIR="${OUTDIR:-parity_out}"
MAXDIST="${MAXDIST:-500}"
HWE="${HWE:-0.001}"
GENO="${GENO:-0.75}"
mkdir -p "$OUTDIR"
# clear stale outputs: Haploview omits the blocks file when there are no blocks,
# so a previous run's file must not linger and contaminate the comparison
rm -f "$OUTDIR"/hv.* "$OUTDIR"/mine.*

# 0. inputs
[ -s "$VCF" ] || bash scripts/fetch_sample_vcf.sh "$VCF" "${NMARK:-300}"
[ -s "$JAR" ] || bash scripts/get_haploview.sh "$JAR"

# 1. convert the SAME genotypes to linkage for the original tool
haploview convert -i "$VCF" --to linkage -o "$OUTDIR/subset"

# 2. original Haploview (nogui), identical thresholds
java -Djava.awt.headless=true -jar "$JAR" -n \
  -pedfile "$OUTDIR/subset.ped" -info "$OUTDIR/subset.info" \
  -dprime -blockoutput GAB \
  -maxdistance "$MAXDIST" -hwcutoff "$HWE" -minGeno "$GENO" \
  -out "$OUTDIR/hv"

# 3. this tool on the VCF with the same thresholds
haploview analyze -i "$VCF" -m gabriel \
  --max-distance "$MAXDIST" --hwe-cutoff "$HWE" --min-geno "$GENO" \
  -o "$OUTDIR/mine"

# 4. compare
python scripts/compare_outputs.py \
  --hv-ld "$OUTDIR/hv.LD" --my-ld "$OUTDIR/mine.LD" \
  --hv-blocks "$OUTDIR/hv.GABRIELblocks" \
  --my-blocks-csv "$OUTDIR/mine.blocks.csv" \
  --info "$OUTDIR/subset.info"
