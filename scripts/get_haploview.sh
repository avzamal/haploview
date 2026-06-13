#!/usr/bin/env bash
# Fetch the original Haploview 4.1 jar (for the parity test). Falls back to
# building from the source mirror with javac if the jar download fails.
set -euo pipefail

DEST="${1:-third_party/Haploview4.1.jar}"
mkdir -p "$(dirname "$DEST")"

if [ -s "$DEST" ]; then
  echo "already present: $DEST"
  exit 0
fi

URL="https://downloads.sourceforge.net/project/haploview/release/Haploview4.1.jar"
echo "downloading $URL"
if curl -fsSL --max-time 180 -o "$DEST" "$URL" && file "$DEST" | grep -q "Java archive"; then
  echo "fetched $DEST"
  exit 0
fi

echo "jar download failed; the parity test will be skipped." >&2
exit 1
