#!/usr/bin/env python3
"""Measure a munged sumstat's allele orientation and positions against gnomAD.

Draws random windows across the genome, reads both the sumstat and the gnomAD sites file
there with tabix (remote `gs://` paths work with `GCS_OAUTH_TOKEN`, which this sets from
gcloud when unset), and classifies every sumstat variant in the windows:

- `match`: `chr:pos:ref:alt` is in gnomAD -- position and orientation agree;
- `swapped`: gnomAD has it as `chr:pos:alt:ref` -- the munge has ref and alt the wrong
  way round, so beta and af are oriented wrong too;
- `absent`: nothing at that position with those alleles -- a variant gnomAD does not
  carry, a differently represented indel, or a wrong position.

For the matches it reports the correlation of the sumstat `af` with the gnomAD AF of the
chosen population and draws the AF-AF plot (SNPs and indels separately). The check is
where the munge scripts that cannot stream gnomAD themselves -- a hundred-million-row
GWAS Catalog file on a small machine -- get the plot the other external GWAS munges draw
inline. It is a sample, not a census: a wrong-strand or wrong-build file shows up in a
few hundred windows; a handful of mis-oriented variants does not.
"""

import argparse
import os
import random
import subprocess
import sys
from collections import Counter

GNOMAD_DEFAULT = "gs://daly-genetics-results/gnomad/gnomad.genomes.exomes.v4.1.1.sites.vep115.tsv.bgz"

# GRCh38 primary assembly lengths, chromosome X coded 23 as in the munged files and gnomAD
CHROMOSOME_LENGTHS = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259, 6: 170805979,
    7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422, 11: 135086622, 12: 133275309,
    13: 114364328, 14: 107043718, 15: 101991189, 16: 90338345, 17: 83257441, 18: 80373285,
    19: 58617616, 20: 64444167, 21: 46709983, 22: 50818468, 23: 156040895,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sumstat", required=True, help="munged, tabix-indexed sumstat (local or gs://)")
    parser.add_argument("--gnomad", default=GNOMAD_DEFAULT, help="gnomAD sites file (default: the daly v4.1.1 copy)")
    parser.add_argument("--af-col", default="AF_nfe", help="gnomAD AF column to compare against (default: AF_nfe)")
    parser.add_argument("--windows", type=int, default=200, help="number of random windows (default: 200)")
    parser.add_argument("--window-size", type=int, default=100_000, help="window width in bp (default: 100000)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--chromosomes", help="comma-separated chromosomes to draw windows from (default: all, X as 23)")
    parser.add_argument("--plot", help="AF-AF plot PNG path (default: <sumstat basename>.af_af.png in the cwd)")
    return parser.parse_args()


def tabix(path: str, region: str) -> list[list[str]]:
    out = subprocess.run(["tabix", path, region], capture_output=True, text=True, check=True).stdout
    return [line.split("\t") for line in out.splitlines() if line]


def header(path: str) -> list[str]:
    out = subprocess.run(["tabix", "-H", path], capture_output=True, text=True, check=True).stdout
    return out.splitlines()[-1].lstrip("#").split("\t")


def random_windows(n: int, width: int, seed: int, chromosomes: list[int] | None) -> list[tuple[int, int, int]]:
    lengths = {c: CHROMOSOME_LENGTHS[c] for c in (chromosomes or CHROMOSOME_LENGTHS)}
    rng = random.Random(seed)
    total = sum(lengths.values())
    windows = []
    for _ in range(n):
        point = rng.randrange(total)
        for chrom, length in lengths.items():
            if point < length:
                start = min(point, length - width)
                windows.append((chrom, start + 1, start + width))
                break
            point -= length
    return windows


def main() -> None:
    args = parse_args()
    if "GCS_OAUTH_TOKEN" not in os.environ and (args.sumstat.startswith("gs://") or args.gnomad.startswith("gs://")):
        os.environ["GCS_OAUTH_TOKEN"] = subprocess.run(
            ["gcloud", "auth", "print-access-token"], capture_output=True, text=True, check=True
        ).stdout.strip()

    sum_cols = header(args.sumstat)
    gn_cols = header(args.gnomad)
    s = {c: sum_cols.index(c) for c in ("chr", "pos", "ref", "alt", "af")}
    g = {c: gn_cols.index(c) for c in ("chr", "pos", "ref", "alt", args.af_col, "genome_or_exome")}

    counts = Counter()
    pairs = []  # (sumstat af, gnomAD af, is_snp)
    chromosomes = [int(c) for c in args.chromosomes.split(",")] if args.chromosomes else None
    for i, (chrom, start, end) in enumerate(random_windows(args.windows, args.window_size, args.seed, chromosomes), 1):
        region = f"{chrom}:{start}-{end}"
        rows = tabix(args.sumstat, region)
        if not rows:
            counts["empty windows"] += 1
            continue
        gnomad = {}
        for r in tabix(args.gnomad, region):
            key = (r[g["pos"]], r[g["ref"]], r[g["alt"]])
            # a variant in both genomes and exomes keeps the genomes row, the larger sample
            if key not in gnomad or r[g["genome_or_exome"]] == "g":
                gnomad[key] = r[g[args.af_col]]
        for r in rows:
            pos, ref, alt = r[s["pos"]], r[s["ref"]], r[s["alt"]]
            if (pos, ref, alt) in gnomad:
                counts["match"] += 1
                try:
                    pairs.append((float(r[s["af"]]), float(gnomad[(pos, ref, alt)]), len(ref) == 1 and len(alt) == 1))
                except ValueError:
                    counts["match without AF"] += 1
            elif (pos, alt, ref) in gnomad:
                counts["swapped"] += 1
            else:
                counts["absent"] += 1
        if i % 20 == 0:
            print(f"  {i} windows: {dict(counts)}", file=sys.stderr, flush=True)

    total = counts["match"] + counts["swapped"] + counts["absent"]
    print(f"windows: {args.windows} x {args.window_size} bp, sumstat variants seen: {total}")
    for key in ("match", "swapped", "absent", "match without AF", "empty windows"):
        if counts[key]:
            share = f" ({counts[key] / total:.1%})" if key in ("match", "swapped", "absent") else ""
            print(f"  {key}: {counts[key]}{share}")

    if pairs:
        import numpy as np

        a = np.array([p[0] for p in pairs])
        b = np.array([p[1] for p in pairs])
        snp = np.array([p[2] for p in pairs])
        print(f"  AF correlation vs gnomAD {args.af_col}: r = {np.corrcoef(a, b)[0, 1]:.4f} on {len(pairs)} matches")
        print(f"  |af - gnomAD| > 0.2: {(abs(a - b) > 0.2).sum()} ({(abs(a - b) > 0.2).mean():.2%})")
        # a mis-oriented row sits on the anti-diagonal; count the common ones there
        common = (b > 0.05) & (b < 0.95)
        anti = common & (abs(a - (1 - b)) < 0.05) & (abs(a - b) > 0.2)
        print(f"  common variants on the anti-diagonal (af ~ 1 - gnomAD): {anti.sum()} of {common.sum()}")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plot = args.plot or os.path.basename(args.sumstat).replace(".tsv.gz", "") + ".af_af.png"
        fig, axes = plt.subplots(1, 2, figsize=(14, 7))
        for ax, mask, label in ((axes[0], snp, "SNPs"), (axes[1], ~snp, "indels")):
            ax.scatter(a[mask], b[mask], s=2, alpha=0.3, rasterized=True)
            ax.plot([0, 1], [0, 1], "k--", linewidth=0.5, alpha=0.5)
            ax.set_xlabel("sumstat af (alt)")
            ax.set_ylabel(f"gnomAD {args.af_col}")
            ax.set_title(f"{label} (n={mask.sum():,})")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_aspect("equal")
        fig.suptitle(os.path.basename(args.sumstat))
        fig.tight_layout()
        fig.savefig(plot, dpi=120)
        print(f"  plot: {plot}")


if __name__ == "__main__":
    main()
