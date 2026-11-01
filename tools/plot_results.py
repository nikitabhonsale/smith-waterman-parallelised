#!/usr/bin/env python3
"""Plot the benchmark CSV written by sw_bench.

Usage:
    python3 tools/plot_results.py bench_results.csv [--outdir plots]

Produces three PNG files in --outdir:
    speedup.png         speedup over the sequential version vs. thread count (one panel per implementation)
    efficiency.png      parallel efficiency (speedup / threads) vs. thread count
    runtime_vs_size.png median fill time vs. problem size, sequential and each parallel version at max threads

Requires matplotlib (pip install matplotlib).
"""
import argparse
import csv
import os
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Okabe-Ito colorblind-safe palette, always assigned in this fixed order.
COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442", "#000000"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]


def load(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "size": int(r["size"]), "impl": r["impl"], "threads": int(r["threads"]),
                "median_ms": float(r["median_ms"]), "speedup": float(r["speedup_vs_seq"]),
                "efficiency": float(r["efficiency"]), "verified": int(r["verified"]),
            })
    return rows


def style(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11)
    ax.grid(True, linewidth=0.5, alpha=0.3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def per_size_lines(ax, rows, impl, value, sizes, thread_ticks):
    for k, n in enumerate(sizes):
        pts = sorted((r["threads"], r[value]) for r in rows if r["impl"] == impl and r["size"] == n)
        if pts:
            ax.plot(*zip(*pts), color=COLORS[k % len(COLORS)], marker=MARKERS[k % len(MARKERS)],
                    linewidth=2, markersize=6, label=f"N = {n}")
    ax.set_xticks(thread_ticks)


def plot_threads(rows, value, ylabel, fname, outdir, ideal):
    sizes = sorted({r["size"] for r in rows})
    ticks = sorted({r["threads"] for r in rows if r["impl"] != "sequential"})
    impls = [i for i in ("wavefront", "tiled") if any(r["impl"] == i for r in rows)]
    fig, axes = plt.subplots(1, len(impls), figsize=(5.5 * len(impls), 4.2), sharey=True, squeeze=False)
    for ax, impl in zip(axes[0], impls):
        per_size_lines(ax, rows, impl, value, sizes, ticks)
        if ideal == "linear":
            ax.plot(ticks, ticks, color="gray", linestyle="--", linewidth=1, label="ideal")
        elif ideal == "one":
            ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, label="ideal")
        style(ax, "OpenMP threads", ylabel, impl)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, fname), dpi=150)
    plt.close(fig)


def plot_runtime(rows, outdir):
    sizes = sorted({r["size"] for r in rows})
    max_t = max(r["threads"] for r in rows)
    series = [("sequential", 1, "sequential (1 thread)")]
    series += [(i, max_t, f"{i} ({max_t} threads)") for i in ("wavefront", "tiled")
               if any(r["impl"] == i for r in rows)]
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for k, (impl, t, label) in enumerate(series):
        pts = sorted((r["size"], r["median_ms"]) for r in rows if r["impl"] == impl and r["threads"] == t)
        if pts:
            ax.plot(*zip(*pts), color=COLORS[k], marker=MARKERS[k], linewidth=2, markersize=6, label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(sizes)
    ax.set_xticklabels([str(s) for s in sizes])
    style(ax, "Sequence length N (N x N matrix)", "Median fill time (ms)", "Runtime vs. problem size")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "runtime_vs_size.png"), dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv")
    p.add_argument("--outdir", default="plots")
    args = p.parse_args()

    rows = load(args.csv)
    if not rows:
        raise SystemExit("CSV has no data rows")
    if any(not r["verified"] for r in rows):
        print("warning: some rows are marked unverified (result differed from the sequential one)")
    os.makedirs(args.outdir, exist_ok=True)

    plot_threads(rows, "speedup", "Speedup over sequential", "speedup.png", args.outdir, "linear")
    plot_threads(rows, "efficiency", "Efficiency (speedup / threads)", "efficiency.png", args.outdir, "one")
    plot_runtime(rows, args.outdir)
    print(f"plots written to {args.outdir}/")


if __name__ == "__main__":
    main()
