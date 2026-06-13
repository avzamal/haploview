#!/usr/bin/env python3
"""Compare original Haploview output against this tool's output.

Asserts the pairwise LD table (D', LOD, r^2, CI bounds, Dist) and the block
membership are identical at Haploview's printed precision. Exits non-zero on any
mismatch.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, List, Tuple


def parse_ld(path: str) -> Dict[Tuple[str, str], Tuple[float, ...]]:
    out = {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 8:
                continue
            l1, l2 = t[0], t[1]
            vals = tuple(float(x) for x in t[2:7])  # D', LOD, r2, CIlow, CIhi
            dist = float(t[7]) if t[7] not in ("-", "") else None
            out[(l1, l2)] = vals + ((dist,) if dist is not None else ())
    return out


def parse_info(path: str) -> List[str]:
    names = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                names.append(line.split()[0])
    return names


def parse_haploview_blocks(path: str, info_names: List[str]) -> List[List[str]]:
    blocks = []
    if not os.path.exists(path):
        return blocks  # Haploview omits the file when there are no blocks
    with open(path) as fh:
        for line in fh:
            if line.startswith("BLOCK"):
                marker_part = line.split("MARKERS:")[1]
                nums = [tok.rstrip("!") for tok in marker_part.split()]
                blocks.append([info_names[int(n) - 1] for n in nums])
    return blocks


def parse_my_blocks_csv(path: str) -> List[List[str]]:
    blocks: Dict[str, List[str]] = {}
    order: List[str] = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            bid = row["block_id"]
            if bid not in blocks:
                blocks[bid] = []
                order.append(bid)
            blocks[bid].append(row["snp_id"])
    return [blocks[b] for b in order]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hv-ld", required=True)
    ap.add_argument("--my-ld", required=True)
    ap.add_argument("--hv-blocks", required=True)
    ap.add_argument("--my-blocks-csv", required=True)
    ap.add_argument("--info", required=True)
    ap.add_argument("--tol", type=float, default=1e-9)
    args = ap.parse_args()

    hv = parse_ld(args.hv_ld)
    mine = parse_ld(args.my_ld)

    mismatches = 0
    common = set(hv) & set(mine)
    only_hv = set(hv) - set(mine)
    only_mine = set(mine) - set(hv)
    for key in sorted(common):
        a, b = hv[key], mine[key]
        for i in range(min(len(a), len(b))):
            if abs(a[i] - b[i]) > args.tol:
                print(f"LD MISMATCH {key} col{i}: haploview={a[i]} mine={b[i]}")
                mismatches += 1

    print(f"LD pairs: {len(common)} compared, {len(only_hv)} only-haploview, "
          f"{len(only_mine)} only-mine, {mismatches} value mismatches")
    if only_hv:
        print("  example only-haploview:", sorted(only_hv)[:3])
    if only_mine:
        print("  example only-mine:", sorted(only_mine)[:3])

    info_names = parse_info(args.info)
    hv_blocks = parse_haploview_blocks(args.hv_blocks, info_names)
    my_blocks = parse_my_blocks_csv(args.my_blocks_csv)

    block_ok = True
    if len(hv_blocks) != len(my_blocks):
        block_ok = False
        print(f"BLOCK COUNT MISMATCH: haploview={len(hv_blocks)} mine={len(my_blocks)}")
    for i, (hb, mb) in enumerate(zip(hv_blocks, my_blocks)):
        if hb != mb:
            block_ok = False
            print(f"BLOCK {i+1} MISMATCH:\n  haploview={hb}\n  mine     ={mb}")
    print(f"Blocks: haploview={len(hv_blocks)} mine={len(my_blocks)} "
          f"{'IDENTICAL' if block_ok else 'DIFFER'}")

    ld_ok = (mismatches == 0 and not only_hv and not only_mine)
    if ld_ok and block_ok:
        print("\nPARITY: PASS -- LD table and blocks are identical.")
        return 0
    print("\nPARITY: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
