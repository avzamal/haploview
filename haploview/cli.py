"""Command-line interface for haploview."""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional

from .io import read_any
from .io.convert import dataset_to_hapmap, dataset_to_linkage
from .pipeline import run_analysis, write_outputs
from .qc import QCParams


def _add_analyze_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--input", "-i", required=True, help="genotype file")
    p.add_argument("--format", "-f", default="auto",
                   choices=["auto", "vcf", "hapmap", "linkage"],
                   help="input format (default: auto-detect)")
    p.add_argument("--info", help=".info marker file (required for linkage)")
    p.add_argument("--method", "-m", default="gabriel",
                   choices=["gabriel", "4gam", "spine", "all"],
                   help="block-definition method (default: gabriel)")
    p.add_argument("--out", "-o", required=True, help="output prefix")
    p.add_argument("--max-distance", type=int, default=500,
                   help="max LD comparison distance in kb (default: 500; 0 = unlimited)")
    p.add_argument("--min-maf", type=float, default=0.0,
                   help="drop markers below this MAF (default: 0)")
    p.add_argument("--min-geno", type=float, default=0.75,
                   help="minimum genotyping rate to keep a marker (default: 0.75)")
    p.add_argument("--hwe-cutoff", type=float, default=0.001,
                   help="drop markers with HWE p below this (default: 0.001; 0 disables)")
    p.add_argument("--hap-thresh", type=float, default=0.01,
                   help="min frequency for a haplotype to be reported (default: 0.01)")
    p.add_argument("--chrom", help="restrict to this chromosome/contig (VCF)")
    p.add_argument("--max-markers", type=int, help="cap number of markers read")


def _run_analyze(args) -> int:
    qc = QCParams(min_geno_percent=args.min_geno, hwe_cutoff=args.hwe_cutoff,
                  min_maf=args.min_maf)
    methods = ["gabriel", "4gam", "spine"] if args.method == "all" else [args.method]
    for method in methods:
        result = run_analysis(
            args.input, fmt=args.format, info=args.info, method=method,
            max_distance_kb=args.max_distance, qc=qc, chrom=args.chrom,
            max_markers=args.max_markers)
        prefix = args.out if len(methods) == 1 else f"{args.out}.{method}"
        paths = write_outputs(result, prefix, method, hap_thresh=args.hap_thresh)
        print(f"[{method}] {result.dataset.n_markers} markers, "
              f"{len(result.dropped)} dropped, {len(result.blocks)} blocks")
        for p in paths:
            print(f"  wrote {p}")
    return 0


def _run_convert(args) -> int:
    ds = read_any(args.input, fmt=args.format, info=args.info)
    if args.to == "linkage":
        dataset_to_linkage(ds, args.out + ".ped", args.out + ".info")
        print(f"wrote {args.out}.ped and {args.out}.info")
    else:
        dataset_to_hapmap(ds, args.out + ".hmp")
        print(f"wrote {args.out}.hmp")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="haploview",
        description="Modern reimplementation of Haploview LD / haplotype-block analysis.")
    sub = parser.add_subparsers(dest="command")

    analyze = sub.add_parser("analyze", help="run LD and block analysis")
    _add_analyze_args(analyze)

    convert = sub.add_parser("convert", help="convert genotypes to legacy formats")
    convert.add_argument("--input", "-i", required=True)
    convert.add_argument("--format", "-f", default="auto",
                         choices=["auto", "vcf", "hapmap", "linkage"])
    convert.add_argument("--info")
    convert.add_argument("--to", choices=["linkage", "hapmap"], default="linkage")
    convert.add_argument("--out", "-o", required=True, help="output prefix")

    # default to "analyze" when no subcommand is given
    args = parser.parse_args(argv)
    if args.command is None:
        analyze_args = analyze.parse_args(argv)
        return _run_analyze(analyze_args)
    if args.command == "analyze":
        return _run_analyze(args)
    if args.command == "convert":
        return _run_convert(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
