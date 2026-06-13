#!/usr/bin/env bash
# Stream the 220 MB Springer VCF and keep the header plus the first N biallelic
# SNP records on contig "01", writing a small committed-size subset. Because the
# file is position-sorted by contig and "01" is first, we exit as soon as we
# leave the contig or hit N markers -- the full file is never downloaded.
set -euo pipefail

URL="https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-021-07768-y/MediaObjects/12864_2021_7768_MOESM24_ESM.vcf"
OUT="${1:-data/chr01_subset.vcf}"
N="${2:-300}"
CONTIG="${3:-01}"

mkdir -p "$(dirname "$OUT")"

curl -s --max-time 600 "$URL" | awk -v n="$N" -v contig="$CONTIG" '
  /^#/ { print; next }
  {
    if ($1 == contig) {
      # keep only simple biallelic SNPs (single-base REF and ALT)
      if (length($4) == 1 && length($5) == 1) { print; c++ }
      if (c >= n) exit
    } else if (seen) {
      exit
    }
  }
  $1 == contig { seen = 1 }
' > "$OUT"

echo "wrote $OUT ($(grep -vc '^#' "$OUT") markers)"
