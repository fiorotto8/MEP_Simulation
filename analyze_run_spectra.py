#!/usr/bin/env python3
"""
Refactored from Analysis_S.ipynb.

Computes flux-density histograms from Geant4 ROOT outputs and produces plots:
1) Entire gas volume (total + per-source)
2) Timepix-only hits (total + per-source)
3) Timepix-only hits with Gaussian smearing (total + per-source)

Requires:
- numpy, awkward, uproot, matplotlib, scipy, tqdm
- your local `analyze_run.py` (imported as ar)
"""

from __future__ import annotations

import argparse
import os
import sys
import csv
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import awkward as ak
import uproot
import matplotlib.pyplot as plt
import scipy.stats as st
from tqdm import tqdm

# Make sure analyze_run.py is importable (same policy as the notebook: sys.path.append("../"))
THIS_DIR = Path(__file__).resolve().parent
sys.path.append(str(THIS_DIR.parent))
import analyze_run as ar  # noqa: E402

def load_theta_limits_from_manifest(manifest_csv: Path) -> dict[str, tuple[float, float]]:
    out = {}
    with manifest_csv.open() as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            fn = row["filename"].strip()
            tmin = float(row["minTheta"])
            tmax = float(row["maxTheta"])
            out[fn] = (tmin, tmax)               # full name
            out[Path(fn).stem] = (tmin, tmax)    # stem
    return out

def compute_f_from_theta(R_mm: float, tmin_deg: float, tmax_deg: float) -> float:
    R_cm = R_mm / 10.0
    tmin = np.deg2rad(tmin_deg)
    tmax = np.deg2rad(tmax_deg)
    tmin = np.clip(tmin, 0.0, np.pi)
    tmax = np.clip(tmax, 0.0, np.pi)
    if tmax <= tmin:
        tmin, tmax = 0.0, np.pi
    A_patch = 2.0 * np.pi * (R_cm**2) * (np.cos(tmin) - np.cos(tmax))
    return np.pi * A_patch

def compute_integrated_rate(y: np.ndarray, sy: np.ndarray, bin_width_keV: float) -> tuple[float, float]:
    """
    Integrated rate R = sum_i y_i * ΔE.
    Uncertainty σ_R from quadrature: sqrt(sum_i (σ_i * ΔE)^2).

    Assumes `y` is a flux density per keV (or per bin_width_keV unit) and `sy` its 1σ uncertainty.
    """
    y = np.asarray(y, dtype=float)
    sy = np.asarray(sy, dtype=float)
    bw = float(bin_width_keV)
    rate = float(np.sum(y * bw))
    srate = float(np.sqrt(np.sum((sy * bw) ** 2)))
    return rate, srate

def save_integrated_rates_csv_per_source(
    out_csv: Path,
    directory: Path,
    energy_min_keV: float,
    energy_max_keV: float,
    bin_width_keV: float,
    s_smearing: float,
    rates_by_case: dict[str, dict[str, tuple[float, float]]],
) -> None:
    """
    Save integrated rates *per source* for multiple cases to CSV.

    `rates_by_case` maps:
        case -> { source_name -> (rate, uncertainty) }

    A "Total" source row can be included in each case dict.
    """
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "case",
        "source",
        "integrated_rate_1_over_s",
        "sigma_1_over_s",
        "directory",
        "energy_min_keV",
        "energy_max_keV",
        "bin_width_keV",
        "s_smearing",
    ]

    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for case, per_source in rates_by_case.items():
            # stable ordering: Total first if present, then alphabetical
            sources = list(per_source.keys())
            ordered = []
            if "Total" in per_source:
                ordered.append("Total")
            ordered += sorted([s for s in sources if s != "Total"])

            for src in ordered:
                r, sr = per_source[src]
                w.writerow(
                    {
                        "case": case,
                        "source": src,
                        "integrated_rate_1_over_s": f"{r:.10g}",
                        "sigma_1_over_s": f"{sr:.10g}",
                        "directory": str(directory),
                        "energy_min_keV": f"{energy_min_keV:.10g}",
                        "energy_max_keV": f"{energy_max_keV:.10g}",
                        "bin_width_keV": f"{bin_width_keV:.10g}",
                        "s_smearing": f"{s_smearing:.10g}",
                    }
                )


def split_by_n(n, vec):
    """
    Split `vec` into consecutive chunks whose lengths are given by `n`.

    Example:
        n = [2, 3, 1], vec = [10,11,12,13,14,15]
        -> [ [10,11], [12,13,14], [15] ]

    Returns a list of numpy arrays (or list slices if vec isn't array-like).
    """
    n = np.asarray(n, dtype=int)
    if n.ndim != 1:
        raise ValueError("`n` must be a 1D array of chunk sizes.")
    if np.any(n < 0):
        raise ValueError("All chunk sizes must be non-negative.")
    if np.sum(n) != len(vec):
        raise ValueError(f"sum(n)={np.sum(n)} does not match len(vec)={len(vec)}")

    out = []
    start = 0
    for size in n:
        out.append(vec[start : start + size])
        start += size
    return out


def compute_histograms(
    directory: Path,
    energy_min_keV: float,
    energy_max_keV: float,
    bin_width_keV: float,
    s_smearing: float,
    spectra_directory: Path,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray],
           np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray],
           np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray]:
    """
    Returns:
      bins,
      (yh_tot, syh_tot, yh_source, syh_source),
      (yh_tot_timepix, syh_tot_timepix, yh_source_timepix, syh_source_timepix),
      (yh_tot_timepix_sm, syh_tot_timepix_sm, yh_source_timepix_sm, syh_source_timepix_sm)
    """
    spectra_tables = ar.load_spectra_table(spectra_directory)

    files = ar.find_root_files(directory)
    if len(files) == 0:
        raise FileNotFoundError(f"No ROOT files found under: {directory}")

    sources = np.unique([info.source for info in files])
    paths: Dict[str, list[Path]] = {s: [] for s in sources}
    print("DEBUG sources:", sources)

    for info in files:
        paths[info.source].append(info.path)

    bins = np.arange(energy_min_keV, energy_max_keV + bin_width_keV, bin_width_keV)

    timepix_side_mm = 14.08

    # totals
    yh_tot = np.zeros(len(bins) - 1)
    syh_tot = np.zeros(len(bins) - 1)
    yh_source: Dict[str, np.ndarray] = {}
    syh_source: Dict[str, np.ndarray] = {}

    # timepix-only
    yh_tot_tp = np.zeros(len(bins) - 1)
    syh_tot_tp = np.zeros(len(bins) - 1)
    yh_source_tp: Dict[str, np.ndarray] = {}
    syh_source_tp: Dict[str, np.ndarray] = {}

    # timepix + smearing
    yh_tot_tp_sm = np.zeros(len(bins) - 1)
    syh_tot_tp_sm = np.zeros(len(bins) - 1)
    yh_source_tp_sm: Dict[str, np.ndarray] = {}
    syh_source_tp_sm: Dict[str, np.ndarray] = {}

    # timepix + smearing + GAGG
    yh_tot_tp_sm_g = np.zeros(len(bins) - 1)
    syh_tot_tp_sm_g = np.zeros(len(bins) - 1)
    yh_source_tp_sm_g: Dict[str, np.ndarray] = {}
    syh_source_tp_sm_g: Dict[str, np.ndarray] = {}

    # --- NEW: normalization from manifest (theta-limited sphere patch) ---
    manifest_csv = spectra_directory / "manifest.csv"
    if not manifest_csv.exists():
        raise FileNotFoundError(f"Missing manifest.csv in spectra_dir: {manifest_csv}")

    theta_map = load_theta_limits_from_manifest(manifest_csv)

    # sphere radius in mm: keep consistent with your generator
    # ideally pass as argument; for now use the same constant as in generation
    sphere_radius_mm = 130.0

    f_by_source: Dict[str, float] = {}
    missing = []
    for s in sources:
        if s not in theta_map:
            missing.append(s)
            continue
        tmin_deg, tmax_deg = theta_map[s]
        f_by_source[s] = compute_f_from_theta(sphere_radius_mm, tmin_deg, tmax_deg)

    if missing:
        raise KeyError(
            "These sources are present in ROOT outputs but missing in manifest.csv mapping: "
            + ", ".join(map(str, missing))
        )


    for s in tqdm(sources, desc="Sources"):
        yh = np.zeros(len(bins) - 1)
        syh = np.zeros(len(bins) - 1)

        yh_tp = np.zeros(len(bins) - 1)
        syh_tp = np.zeros(len(bins) - 1)

        yh_tp_sm = np.zeros(len(bins) - 1)
        syh_tp_sm = np.zeros(len(bins) - 1)

        yh_tp_sm_g = np.zeros(len(bins) - 1)
        syh_tp_sm_g = np.zeros(len(bins) - 1)

        for thisp in paths[s]:
            with uproot.open(thisp) as thisfile:
                # skip files without the expected tree
                if "Events;1" not in thisfile.keys():
                    continue

                thistree = thisfile["Events"]
                enes = thistree["edepGasTotal_keV"].array()  # awkward array
                nshots = thisfile["RunInfo"]["nBeamOnRequested"].array()[0]
                nhdep = np.array(thisfile["Events;1"]["nHitsEdep"].array())

                # hit coordinates split per-event
                hits_x = split_by_n(nhdep, thisfile["Hits;1"]["x_mm"].array())
                hits_y = split_by_n(nhdep, thisfile["Hits;1"]["y_mm"].array())
                hits_z = split_by_n(nhdep, thisfile["Hits;1"]["z_mm"].array())

                GAGGHit = thistree["GAGGHit"].array() #np.zeros(len(hits_x))#thistree["GAGGHit"].array()

                on_timepix = np.zeros(len(hits_x), dtype=int)
                half = timepix_side_mm / 2.0

                on_gagg = np.zeros(len(hits_x), dtype=int)

                for ev in range(len(hits_x)):
                    x = hits_x[ev]
                    y = hits_y[ev]
                    # z currently unused in selection, kept for parity with notebook
                    _z = hits_z[ev]
                    if np.any((x >= -half) & (x <= half) & (y >= -half) & (y <= half)):
                        on_timepix[ev] = 1
                    if GAGGHit[ev] != 0:
                        on_gagg[ev] = 1
                    

            if int(nshots) <= 0:
                continue

            if len(on_timepix) != len(enes):
                raise RuntimeError(
                    f"Number of events with hits ({len(on_timepix)}) != number of energy depositions ({len(enes)}) "
                    f"for file {thisp}"
                )

            # source spectrum integral (as in notebook)
            if s not in spectra_tables:
                raise KeyError(f"Source '{s}' not found in spectra table loaded from {spectra_directory}")

            Egrid, Fgrid = spectra_tables[s]

            edges = np.empty(Egrid.size + 1)
            edges[1:-1] = np.sqrt(Egrid[:-1] * Egrid[1:])  # geometric midpoints (log grid)

            r0 = Egrid[0] / Egrid[1]
            rN = Egrid[-1] / Egrid[-2]
            edges[0] = Egrid[0] * np.sqrt(r0)
            edges[-1] = Egrid[-1] * np.sqrt(rN)

            dE = edges[1:] - edges[:-1]
            flux_int = np.sum(Fgrid * dE)

            # histogram (all events)
            y_tmp, _ = np.histogram(enes, bins=bins)
            y_tmp = np.asarray(y_tmp)

            # normalization (same as notebook)
            #f = 4.0 * np.pi * (ar.SPHERE_RADIUS_MM / 10.0) ** 2 * np.pi
            f = f_by_source[s]
                        
            pint = (y_tmp + 1) / (nshots + 2)
            yh += pint * flux_int * f / bin_width_keV
            syh += (pint * (1.0 - pint) / (nshots + 3)) * (flux_int * f / bin_width_keV) ** 2

            # timepix-only
            mask_tp = (on_timepix == 1)
            y_tmp_tp, _ = np.histogram(enes[mask_tp], bins=bins)
            y_tmp_tp = np.asarray(y_tmp_tp)
            pint_tp = (y_tmp_tp + 1) / (nshots + 2)
            yh_tp += pint_tp * flux_int * f / bin_width_keV
            syh_tp += (pint_tp * (1.0 - pint_tp) / (nshots + 3)) * (flux_int * f / bin_width_keV) ** 2

            # timepix + smearing
            sm_enes = ak.to_numpy(enes)
            sm_enes = st.norm.rvs(loc=sm_enes, scale=s_smearing * sm_enes)
            
            y_tmp_tp_sm, _ = np.histogram(sm_enes[mask_tp], bins=bins)
            y_tmp_tp_sm = np.asarray(y_tmp_tp_sm)
            pint_tp_sm = (y_tmp_tp_sm + 1) / (nshots + 2)
            yh_tp_sm += pint_tp_sm * flux_int * f / bin_width_keV
            syh_tp_sm += (pint_tp_sm * (1.0 - pint_tp_sm) / (nshots + 3)) * (flux_int * f / bin_width_keV) ** 2

            # timepix + smearing + GAGG
            mask_tp_g = (on_timepix == 1) & (on_gagg == 0)

            y_tmp_tp_g, _ = np.histogram(sm_enes[mask_tp_g], bins=bins)
            y_tmp_tp_g = np.asarray(y_tmp_tp_g)
            pint_tp_g = (y_tmp_tp_g + 1) / (nshots + 2)
            yh_tp_sm_g += pint_tp_g * flux_int * f / bin_width_keV
            syh_tp_sm_g += (pint_tp_g * (1.0 - pint_tp_g) / (nshots + 3)) * (flux_int * f / bin_width_keV) ** 2

        yh_source[s] = yh
        syh_source[s] = np.sqrt(syh)
        yh_tot += yh
        syh_tot += syh

        yh_source_tp[s] = yh_tp
        syh_source_tp[s] = np.sqrt(syh_tp)
        yh_tot_tp += yh_tp
        syh_tot_tp += syh_tp

        yh_source_tp_sm[s] = yh_tp_sm
        syh_source_tp_sm[s] = np.sqrt(syh_tp_sm)
        yh_tot_tp_sm += yh_tp_sm
        syh_tot_tp_sm += syh_tp_sm

        yh_source_tp_sm_g[s] = yh_tp_sm_g
        syh_source_tp_sm_g[s] = np.sqrt(syh_tp_sm_g)
        yh_tot_tp_sm_g += yh_tp_sm_g
        syh_tot_tp_sm_g += syh_tp_sm_g

    # totals' sigma
    syh_tot = np.sqrt(syh_tot)
    syh_tot_tp = np.sqrt(syh_tot_tp)
    syh_tot_tp_sm = np.sqrt(syh_tot_tp_sm)
    syh_tot_tp_sm_g = np.sqrt(syh_tot_tp_sm_g)

    return (
        bins,
        yh_tot, syh_tot, yh_source, syh_source,
        yh_tot_tp, syh_tot_tp, yh_source_tp, syh_source_tp,
        yh_tot_tp_sm, syh_tot_tp_sm, yh_source_tp_sm, syh_source_tp_sm,
        yh_tot_tp_sm_g, syh_tot_tp_sm_g, yh_source_tp_sm_g, syh_source_tp_sm_g,
        sources,
    )


def plot_with_band(
    outpath: Path,
    title: str,
    bins: np.ndarray,
    y_tot: np.ndarray,
    sy_tot: np.ndarray,
    y_by_source: Dict[str, np.ndarray],
    sy_by_source: Dict[str, np.ndarray],
    sources,
    bin_width_keV: float,
    y_label: str = r"Flux density in detector [s$^{-1}$ keV$^{-1}$]",
    y_min: float = 1e-2,
):
    integrated_flux = float(np.sum(y_tot * bin_width_keV))
    sintegrated_flux = float(np.sqrt(np.sum((sy_tot**2) * (bin_width_keV**2))))

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.stairs(y_tot, bins, label="Total", linewidth=2)

    # 1σ band (pad to draw the last bin with step='post')
    y_low = np.maximum(y_tot - sy_tot, 1e-30)
    y_high = np.maximum(y_tot + sy_tot, 1e-30)
    y_low_p = np.r_[y_low, y_low[-1]]
    y_high_p = np.r_[y_high, y_high[-1]]

    ax.fill_between(bins, y_low_p, y_high_p, step="post", alpha=0.25, label=r"$\pm 1\sigma$")

    for s in sources:
        if s in y_by_source:
            ax.stairs(y_by_source[s], bins, label=str(s), linewidth=1)
            y_low = np.maximum(y_by_source[s] - sy_by_source[s], 1e-30)
            y_high = np.maximum(y_by_source[s] + sy_by_source[s], 1e-30)
            y_low_p = np.r_[y_low, y_low[-1]]
            y_high_p = np.r_[y_high, y_high[-1]]

            ax.fill_between(bins, y_low_p, y_high_p, step="post", alpha=0.25, label=None)

    ax.set_title(title)
    ax.set_yscale("log")
    ax.set_ylabel(y_label)
    ax.set_ylim(bottom=1e-4)
    ax.legend(loc="lower right")

    # place text at ~top-left; avoid clipping
    ax.text(
        0.02, 0.95,
        f"Integrated flux = ({integrated_flux:.2f} ± {sintegrated_flux:.2f}) 1/S",
        transform=ax.transAxes,
        fontsize=12,
        va="top",
        ha="left",
    )

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Refactored Analysis_S notebook runner.")
    p.add_argument("--directory", type=Path, help="Path containing run output ROOT files.", default=Path("run_outputs/"))
    p.add_argument("--energy-min", type=float, required=False, help="Energy cut min (keV).", default=0.0)
    p.add_argument("--energy-max", type=float, required=False, help="Energy cut max (keV).", default=100.0)
    p.add_argument("--bin-width", type=float, required=False, help="Histogram bin width (keV).", default=1.0)
    p.add_argument("--s-smearing", type=float, required=False, help="Gaussian relative smearing (e.g. 0.15 for 15 percent).", default=0.15)
    # Not requested, but useful and keeps notebook behavior reproducible:
    p.add_argument(
        "--spectra-dir",
        type=Path,
        default=(THIS_DIR.parent / "gpd3d/spectra"),
        help="Directory containing spectra tables (default: ../spectra relative to script).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    directory: Path = args.directory
    if not directory.exists():
        raise FileNotFoundError(f"Directory does not exist: {directory}")

    bins, yh_tot, syh_tot, yh_source, syh_source, \
        yh_tot_tp, syh_tot_tp, yh_source_tp, syh_source_tp, \
        yh_tot_tp_sm, syh_tot_tp_sm, yh_source_tp_sm, syh_source_tp_sm, \
        yh_tot_tp_sm_g, syh_tot_tp_sm_g, yh_source_tp_sm_g, syh_source_tp_sm_g, \
        sources = compute_histograms(
            directory=directory,
            energy_min_keV=args.energy_min,
            energy_max_keV=args.energy_max,
            bin_width_keV=args.bin_width,
            s_smearing=args.s_smearing,
            spectra_directory=args.spectra_dir,
        )

    plot_dir: Path = args.directory / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    plot_with_band(
        outpath=plot_dir / "flux_gas_volume_{0:.1f}_{1:.1f}_{2:.2f}.png".format(args.energy_min, args.energy_max, args.bin_width),
        title="Simulated energy deposition in gas - entire gas volume",
        bins=bins,
        y_tot=yh_tot,
        sy_tot=syh_tot,
        y_by_source=yh_source,
        sy_by_source=syh_source,
        sources=sources,
        bin_width_keV=args.bin_width,
    )

    plot_with_band(
        outpath=plot_dir / "flux_timepix_{0:.1f}_{1:.1f}_{2:.2f}.png".format(args.energy_min, args.energy_max, args.bin_width),
        title="Simulated energy deposition in gas - hits on Timepix",
        bins=bins,
        y_tot=yh_tot_tp,
        sy_tot=syh_tot_tp,
        y_by_source=yh_source_tp,
        sy_by_source=syh_source_tp,
        sources=sources,
        bin_width_keV=args.bin_width,
    )

    plot_with_band(
        outpath=plot_dir / "flux_timepix_smeared_{0:.1f}_{1:.1f}_{2:.2f}_{3:.2f}.png".format(args.energy_min,
                                                                                             args.energy_max,
                                                                                             args.bin_width,
                                                                                             args.s_smearing
                                                                                             ),
        title=f"Simulated energy deposition in gas - Timepix smeared (σ={args.s_smearing:.3g})",
        bins=bins,
        y_tot=yh_tot_tp_sm,
        sy_tot=syh_tot_tp_sm,
        y_by_source=yh_source_tp_sm,
        sy_by_source=syh_source_tp_sm,
        sources=sources,
        bin_width_keV=args.bin_width,
    )

    plot_with_band(
        outpath=plot_dir / "flux_timepix_smeared_gagg_{0:.1f}_{1:.1f}_{2:.2f}_{3:.2f}.png".format(args.energy_min,
                                                                                             args.energy_max,
                                                                                             args.bin_width,
                                                                                             args.s_smearing
                                                                                             ),
        title=f"Simulated energy deposition in gas - Timepix smeared + GAGG (σ={args.s_smearing:.3g})",
        bins=bins,
        y_tot=yh_tot_tp_sm_g,
        sy_tot=syh_tot_tp_sm_g,
        y_by_source=yh_source_tp_sm_g,
        sy_by_source=syh_source_tp_sm_g,
        sources=sources,
        bin_width_keV=args.bin_width,
    )

    # --- integrated rates per source (3 cases) ---
    rates_by_case: dict[str, dict[str, tuple[float, float]]] = {}

    # standard
    per_source_std = {s: compute_integrated_rate(yh_source[s], syh_source[s], args.bin_width) for s in sources}
    per_source_std["Total"] = compute_integrated_rate(yh_tot, syh_tot, args.bin_width)
    rates_by_case["standard"] = per_source_std

    # timepix cut
    per_source_tp = {s: compute_integrated_rate(yh_source_tp[s], syh_source_tp[s], args.bin_width) for s in sources}
    per_source_tp["Total"] = compute_integrated_rate(yh_tot_tp, syh_tot_tp, args.bin_width)
    rates_by_case["timepix_cut"] = per_source_tp

    # timepix + diffusion (smearing)
    per_source_tp_sm = {
        s: compute_integrated_rate(yh_source_tp_sm[s], syh_source_tp_sm[s], args.bin_width) for s in sources
    }
    per_source_tp_sm["Total"] = compute_integrated_rate(yh_tot_tp_sm, syh_tot_tp_sm, args.bin_width)
    rates_by_case["timepix_diffusion"] = per_source_tp_sm

    # timepix + diffusion (smearing) + GAGG
    per_source_tp_sm_g = {
        s: compute_integrated_rate(yh_source_tp_sm_g[s], syh_source_tp_sm_g[s], args.bin_width) for s in sources
    }
    per_source_tp_sm_g["Total"] = compute_integrated_rate(yh_tot_tp_sm_g, syh_tot_tp_sm_g, args.bin_width)
    rates_by_case["timepix_diffusion_gagg"] = per_source_tp_sm_g

    rates_csv = plot_dir / "integrated_rates_{0:.1f}_{1:.1f}_{2:.2f}_{3:.2f}.csv".format(
        args.energy_min, args.energy_max, args.bin_width, args.s_smearing
    )
    save_integrated_rates_csv_per_source(
        out_csv=rates_csv,
        directory=directory,
        energy_min_keV=args.energy_min,
        energy_max_keV=args.energy_max,
        bin_width_keV=args.bin_width,
        s_smearing=args.s_smearing,
        rates_by_case=rates_by_case,
    )


    print(f"Wrote plots to: {plot_dir}")
    print(f"Wrote integrated rates to: {rates_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
