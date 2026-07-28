#!/usr/bin/env python3
"""
plot_from_run_outputs.py

Standalone plotter + recap table for the outputs produced by your rate script:
  - run_outputs/run_points.csv  (required)
  - run_outputs/run_points.json (optional, not used)

It produces:
  - run_outputs/rates_overview.png  (Ri(E) per spectrum, log-log)
  - run_outputs/phit_overview.png   (P(hit)(E) per spectrum, semilogx)

It also prints a recap table to terminal:
  - Rate per spectrum file (source_file + particle)
  - Rate summed by particle type
  - TOTAL rate (sum of all files)

Usage:
  python3 plot_from_run_outputs.py
  python3 plot_from_run_outputs.py --csv run_outputs/run_points.csv --outdir run_outputs
"""

import argparse
import csv
import math
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")  # headless backend (no DISPLAY needed)

import matplotlib.pyplot as plt

DEFAULT_OUTDIR = Path("run_outputs")
DEFAULT_CSV = DEFAULT_OUTDIR / "run_points.csv"

def read_points_csv(csv_path: Path):
    rows = []
    with open(csv_path, newline="") as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            def ffloat(key, default=float("nan")):
                try:
                    return float(r.get(key, default))
                except (ValueError, TypeError):
                    return default

            def fint(key, default=0):
                try:
                    return int(float(r.get(key, default)))
                except (ValueError, TypeError):
                    return default

            rows.append({
                "source_file": (r.get("source_file", "") or "").strip(),
                "particle": (r.get("particle", "") or "").strip(),
                "index": fint("index", 0),
                "E_MeV": ffloat("E_MeV"),
                "R_i_per_s": ffloat("R_i_per_s", 0.0),
                "P_hit": ffloat("P_hit"),
                "N_gen": fint("N_gen", 0),
                "N_hit": fint("N_hit", 0),
                "Phi_bin": ffloat("Phi_bin", float("nan")),
                "dE_MeV": ffloat("dE_MeV", float("nan")),
                "flux_differential": ffloat("flux_differential", float("nan")),
            })
    return rows

def group_rows(rows):
    groups = defaultdict(list)
    for r in rows:
        key = (r["source_file"], r["particle"])
        groups[key].append(r)
    for k in list(groups.keys()):
        groups[k] = sorted(groups[k], key=lambda x: x["E_MeV"] if x["E_MeV"] == x["E_MeV"] else float("inf"))
    return groups

def plot_rates_overview(groups, outpath: Path):
    fig = plt.figure()
    ax = plt.gca()

    any_plotted = False
    for (src, part), rows in groups.items():
        E = [x["E_MeV"] for x in rows if x["E_MeV"] > 0 and x["R_i_per_s"] > 0]
        R = [x["R_i_per_s"] for x in rows if x["E_MeV"] > 0 and x["R_i_per_s"] > 0]
        if len(E) >= 2:
            ax.loglog(E, R, label=f"{Path(src).stem} [{part}]")
            any_plotted = True

    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("Rate contribution per point Ri (1/s)")
    ax.set_title("Rate contributions vs energy (log-log)")
    ax.grid(True, which="both")
    if any_plotted:
        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(str(outpath), dpi=200)
    plt.close(fig)

def plot_phit_overview(groups, outpath: Path):
    fig = plt.figure()
    ax = plt.gca()

    any_plotted = False
    for (src, part), rows in groups.items():
        E = [x["E_MeV"] for x in rows if x["E_MeV"] > 0 and (x["P_hit"] == x["P_hit"])]
        P = [x["P_hit"] for x in rows if x["E_MeV"] > 0 and (x["P_hit"] == x["P_hit"])]
        if len(E) >= 2:
            ax.semilogx(E, P, label=f"{Path(src).stem} [{part}]")
            any_plotted = True

    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("P(hit)  (EdepGas>0)")
    ax.set_title("Hit probability vs energy")
    ax.grid(True, which="both")
    if any_plotted:
        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(str(outpath), dpi=200)
    plt.close(fig)

def compute_rates(rows):
    """
    Your CSV is point-by-point. The per-spectrum rate is just:
      rate(source_file,particle) = sum_i R_i_per_s
    Then we sum by particle and overall.
    """
    rate_by_component = defaultdict(float)  # key=(source_file, particle)
    rate_by_particle = defaultdict(float)   # key=particle
    total_rate = 0.0

    for r in rows:
        Ri = r["R_i_per_s"]
        if Ri != Ri:  # NaN
            continue
        key = (r["source_file"], r["particle"])
        rate_by_component[key] += Ri
        rate_by_particle[r["particle"]] += Ri
        total_rate += Ri

    return rate_by_component, rate_by_particle, total_rate

def print_recap(rate_by_component, rate_by_particle, total_rate):
    # Component table
    comps = sorted(rate_by_component.items(), key=lambda kv: kv[1], reverse=True)

    print("\n================= RATE RECAP =================")
    print("Per component (spectrum file):")
    print(f"{'source_file':30s} {'particle':10s} {'rate [1/s]':>14s}")
    print("-" * 58)
    for (src, part), rate in comps:
        print(f"{src:30s} {part:10s} {rate:14.6g}")
    print("-" * 58)

    # Particle-summed table
    parts = sorted(rate_by_particle.items(), key=lambda kv: kv[1], reverse=True)

    print("\nSummed by particle:")
    print(f"{'particle':10s} {'rate [1/s]':>14s} {'fraction':>12s}")
    print("-" * 40)
    if total_rate > 0:
        for part, rate in parts:
            frac = rate / total_rate
            print(f"{part:10s} {rate:14.6g} {frac:12.4f}")
    else:
        for part, rate in parts:
            print(f"{part:10s} {rate:14.6g} {'-':>12s}")
    print("-" * 40)

    print(f"\nTOTAL RATE: {total_rate:.6g} 1/s")
    print("=============================================\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=str, default=str(DEFAULT_CSV), help="Path to run_points.csv")
    ap.add_argument("--outdir", type=str, default=str(DEFAULT_OUTDIR), help="Output directory for plots")
    args = ap.parse_args()

    csv_path = Path(args.csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not csv_path.exists():
        raise FileNotFoundError(f"Missing CSV: {csv_path}")

    rows = read_points_csv(csv_path)
    if not rows:
        raise RuntimeError("CSV contained no readable rows.")

    groups = group_rows(rows)

    # Plots
    out_rates = outdir / "rates_overview.png"
    out_phit = outdir / "phit_overview.png"
    plot_rates_overview(groups, out_rates)
    plot_phit_overview(groups, out_phit)

    # Recap tables
    rate_by_component, rate_by_particle, total_rate = compute_rates(rows)
    print_recap(rate_by_component, rate_by_particle, total_rate)

    print(f"Saved: {out_rates}")
    print(f"Saved: {out_phit}")

if __name__ == "__main__":
    main()
