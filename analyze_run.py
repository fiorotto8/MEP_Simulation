#!/usr/bin/env python3
import os
import re
import math
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from optparse import OptionParser
from collections import defaultdict
import numpy as np
import uproot

# Optional but strongly recommended for jagged arrays output
try:
    import awkward as ak
except Exception as e:
    ak = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# Config
# -----------------------------
SPHERE_RADIUS_MM = 130.0

BASE_DIR = Path("run_outputs")
SPECTRA_DIR = Path("spectra")
PLOTS_DIR = BASE_DIR / "plots"
#PLOTS_DIR.mkdir(parents=True, exist_ok=True)

OUTFILE = Path("summary_all_branches.root")

ENERGY_RE = re.compile(r"E([\d.]+)MeV")


# -----------------------------
# Logging
# -----------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
log = logging.getLogger("analyze_run")


# -----------------------------
# Data structures
# -----------------------------
@dataclass(frozen=True)
class FileMeta:
    source: str
    particle: str
    energy_MeV: float
    filename: str
    path: Path


@dataclass
class FileData:
    n_beam: int
    edep_keV: np.ndarray  # 1D array
    num_hits: int
    hit_prob: float
    err_hit_prob: float


# -----------------------------
# File discovery / parsing
# -----------------------------
def parse_root_filename(filename: str) -> Optional[Tuple[str, str, float]]:
    """
    Expected pattern: <source>_<particle>_..._E<energy>MeV.root
    Returns (source, particle, energy_MeV) or None if it can't parse.
    """
    if not filename.endswith(".root"):
        return None
    stem = filename[:-5]
    parts = stem.split("_")
    if len(parts) < 2:
        return None

    source = parts[0]
    particle = parts[1]

    m = ENERGY_RE.search(stem)
    energy = float(m.group(1)) if m else float("nan")
    return source, particle, energy


def find_root_files(base_dir: Path = BASE_DIR) -> List[FileMeta]:
    out: List[FileMeta] = []
    for root, _, files in os.walk(base_dir):
        rootp = Path(root)
        for fn in files:
            if not fn.endswith(".root"):
                continue
            parsed = parse_root_filename(fn)
            if not parsed:
                continue
            source, particle, energy = parsed
            out.append(FileMeta(source, particle, energy, fn, rootp / fn))

    out.sort(key=lambda x: (x.source, x.particle, x.energy_MeV, x.filename))
    return out


# -----------------------------
# Spectra utilities
# -----------------------------
def compute_dE_per_entry(out_source: List[str], out_energy: List[float]) -> np.ndarray:
    """
    For each entry (file), compute an effective energy step dE (MeV)
    based on the spacing of simulated energies within each source.
    Midpoint rule for interior points, one-sided for edges.
    """
    N = len(out_source)
    dE = np.zeros(N, dtype=np.float64)

    by_src = defaultdict(list)
    for i, src in enumerate(out_source):
        by_src[src].append(i)

    for src, idxs in by_src.items():
        idxs = sorted(idxs, key=lambda i: out_energy[i])
        E = np.array([out_energy[i] for i in idxs], dtype=np.float64)

        if len(E) < 2:
            dE[idxs[0]] = 0.0
            continue

        dEi = np.empty_like(E)
        dEi[0]  = E[1] - E[0]
        dEi[-1] = E[-1] - E[-2]
        if len(E) > 2:
            dEi[1:-1] = 0.5 * (E[2:] - E[:-2])

        for j, i in enumerate(idxs):
            dE[i] = float(dEi[j])

    return dE


def load_spectra_table(spectra_dir: Path = SPECTRA_DIR) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Returns dict: source_name -> (E_sorted, F_sorted)
    CSV format: E, fluxSpace  (comma-separated)
    Key is filename without extension.
    """
    tables: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for fn in spectra_dir.iterdir():
        if fn.suffix.lower() != ".csv":
            continue
        if fn.name == "manifest.csv":
            continue

        arr = np.loadtxt(fn, delimiter=",")
        if arr.ndim == 1 and arr.size == 2:
            arr = arr.reshape(1, 2)

        E = arr[:, 0].astype(float)
        F = arr[:, 1].astype(float)

        idx = np.argsort(E)
        tables[fn.stem] = (E[idx], F[idx])

    return tables


def interp_flux_loglog(E_grid: np.ndarray, F_grid: np.ndarray, E_query: float) -> float:
    """
    Log-log interpolation. NaN if out of range or invalid.
    Falls back to linear interpolation if logs not possible.
    """
    if not np.isfinite(E_query) or E_query <= 0:
        return float("nan")
    if E_grid.size < 2:
        return float("nan")
    if E_query < E_grid[0] or E_query > E_grid[-1]:
        return float("nan")

    if np.any(E_grid <= 0) or np.any(F_grid <= 0):
        return float(np.interp(E_query, E_grid, F_grid))

    logE = np.log(E_grid)
    logF = np.log(F_grid)
    val = np.interp(np.log(E_query), logE, logF)
    return float(np.exp(val))


def compute_rate(flux_space: float, sphere_radius_mm: float, prob_hit: float, err_prob_hit: float) -> Tuple[float, float]:
    """
    Your original formula:
      area_cm2 = 4*pi*(R_cm)^2
      rate = pi * fluxSpace * area_cm2 * prob_hit
    """
    if not np.isfinite(flux_space):
        return float("nan"), float("nan")
    area_cm2 = 4.0 * math.pi * (sphere_radius_mm / 10.0) ** 2
    rate = math.pi * flux_space * area_cm2 * prob_hit
    err_rate = math.pi * flux_space * area_cm2 * err_prob_hit
    return rate, err_rate


# -----------------------------
# ROOT reading with uproot
# -----------------------------
def read_file_data(path: Path) -> FileData:
    """
    Reads RunInfo/nBeamOnRequested and Events/edepGasTotal_keV.
    Returns hit probability and binomial error.
    """
    with uproot.open(path) as f:
        n_beam = -1
        if "RunInfo" in f and "nBeamOnRequested" in f["RunInfo"]:
            arr = f["RunInfo"]["nBeamOnRequested"].array(library="np")
            if len(arr):
                n_beam = int(arr[0])

        edep = np.empty(0, dtype=np.float64)
        if "Events" in f and "edepGasTotal_keV" in f["Events"]:
            edep = f["Events"]["edepGasTotal_keV"].array(library="np")

        num_hits = int(len(edep))
        hit_prob = (num_hits / n_beam) if n_beam > 0 else 0.0
        err_hit_prob = math.sqrt(hit_prob * (1.0 - hit_prob) / n_beam) if n_beam > 0 else 0.0

        return FileData(n_beam=n_beam, edep_keV=edep, num_hits=num_hits, hit_prob=hit_prob, err_hit_prob=err_hit_prob)


# -----------------------------
# Plotting
# -----------------------------
def plot_all_flux_spectra(spectra_tables: Dict[str, Tuple[np.ndarray, np.ndarray]]) -> None:
    if not spectra_tables:
        return
    plt.figure()
    for src, (E, F) in spectra_tables.items():
        plt.loglog(E, F, marker="o", linestyle="-", label=src)
    plt.xlabel("Primary Particle Energy (MeV)")
    plt.ylabel("Flux Spectra")
    plt.title("Flux Spectra for All Sources")
    plt.legend(fontsize="small")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.savefig(PLOTS_DIR / "flux_spectra_all_sources.png", dpi=150, bbox_inches="tight")
    plt.close()


def plot_hitprob_vs_energy(per_source: Dict[str, Dict[str, List[float]]]) -> None:
    if not per_source:
        return
    plt.figure()
    for src, d in per_source.items():
        plt.errorbar(d["energy"], d["hitProb"], yerr=d["err_hitProb"], fmt="o", label=src)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Primary Particle Energy (MeV)")
    plt.ylabel("Hit Probability")
    plt.title("Hit Probability vs Energy (per source)")
    plt.legend(fontsize="small")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.savefig(PLOTS_DIR / "hit_probability_vs_energy.png", dpi=150, bbox_inches="tight")
    plt.close()


def plot_edep_histograms_weighted(
    sources: List[str],
    out_source: List[str],
    out_edep_lists: List[List[float]],
    out_numhits: List[int],
    w_perHit: np.ndarray,
    plots_dir: Path,
) -> None:

    # find global min/max over all edep values (for common bins)
    global_min, global_max = None, None
    for ed in out_edep_lists:
        if len(ed) == 0:
            continue
        mn, mx = float(np.min(ed)), float(np.max(ed))
        global_min = mn if global_min is None else min(global_min, mn)
        global_max = mx if global_max is None else max(global_max, mx)

    if global_min is None:
        log.warning("No edepGasTotal_keV values found for any source")
        return

    bins = np.linspace(0.99 * global_min, 1.01 * global_max, 101)  # 100 bins
    binw = np.diff(bins)

    total_hist = np.zeros(len(bins) - 1, dtype=np.float64)

    # --- Per-source weighted spectra ---
    for src in sources:
        h_src = np.zeros(len(bins) - 1, dtype=np.float64)

        for i, s in enumerate(out_source):
            if s != src:
                continue
            if out_numhits[i] <= 0:
                continue
            if not np.isfinite(w_perHit[i]) or w_perHit[i] <= 0:
                continue
            ed = np.asarray(out_edep_lists[i], dtype=np.float64)
            if ed.size == 0:
                continue

            # each hit in this entry gets the same weight
            weights = np.full(ed.shape, w_perHit[i], dtype=np.float64)
            hi, _ = np.histogram(ed, bins=bins, weights=weights)
            h_src += hi

        if np.all(h_src == 0):
            log.warning("Weighted spectrum is empty for source %s", src)
            continue

        total_hist += h_src

        # Plot as a density (per keV) in "flux units" (whatever fluxSpace units are)
        h_src_per_keV = h_src / binw

        plt.figure()
        plt.bar(bins[:-1], h_src_per_keV, width=binw, alpha=0.7, align="edge")
        plt.title(f"Flux-weighted Edep spectrum: {src}")
        plt.xlabel("Energy Deposition (keV)")
        plt.ylabel("Weighted counts / keV")
        plt.grid(True, linestyle="--", linewidth=0.5)
        plt.savefig(plots_dir / f"flux_weighted_edep_spectrum_{src}.png", dpi=150, bbox_inches="tight")
        plt.close()

        # Optional: also save a shape-only normalized PDF (area=1)
        area = float(np.sum(h_src_per_keV * binw))
        if area > 0:
            pdf_per_keV = h_src_per_keV / area
            plt.figure()
            plt.bar(bins[:-1], pdf_per_keV, width=binw, alpha=0.7, align="edge")
            plt.title(f"Flux-weighted Edep PDF (area=1): {src}")
            plt.xlabel("Energy Deposition (keV)")
            plt.ylabel("PDF / keV")
            plt.grid(True, linestyle="--", linewidth=0.5)
            plt.savefig(plots_dir / f"flux_weighted_edep_pdf_{src}.png", dpi=150, bbox_inches="tight")
            plt.close()

    # --- Total combined weighted spectrum (all sources) ---
    if np.any(total_hist > 0):
        total_per_keV = total_hist / binw

        plt.figure()
        plt.bar(bins[:-1], total_per_keV, width=binw, alpha=0.7, align="edge")
        plt.title("Total flux-weighted Edep spectrum (all sources)")
        plt.xlabel("Energy Deposition (keV)")
        plt.ylabel("Weighted counts / keV")
        plt.grid(True, linestyle="--", linewidth=0.5)
        plt.savefig(plots_dir / "total_flux_weighted_edep_spectrum.png", dpi=150, bbox_inches="tight")
        plt.close()

        area = float(np.sum(total_per_keV * binw))
        if area > 0:
            total_pdf_per_keV = total_per_keV / area
            plt.figure()
            plt.bar(bins[:-1], total_pdf_per_keV, width=binw, alpha=0.7, align="edge")
            plt.title("Total flux-weighted Edep PDF (area=1, all sources)")
            plt.xlabel("Energy Deposition (keV)")
            plt.ylabel("PDF / keV")
            plt.grid(True, linestyle="--", linewidth=0.5)
            plt.savefig(plots_dir / "total_flux_weighted_edep_pdf.png", dpi=150, bbox_inches="tight")
            plt.close()


# -----------------------------
# Writing output with uproot (no PyROOT)
# -----------------------------
def write_summary_root(outfile: Path, rows: Dict[str, object]) -> None:
    """
    Writes a TTree named FileSummary to `outfile`.
    `rows` is a dict of branchname -> array-like (numpy arrays or awkward arrays).
    """
    with uproot.recreate(outfile) as fout:
        tree = fout.mktree("FileSummary", rows)  # creates a TTree (not an RNTuple)
        tree.extend(rows)          


# -----------------------------
# Main pipeline
# -----------------------------
def main() -> None:

    directory         = None
    spectra_directory = None
    outfile           = None
    plots_dir         = None
    parser = OptionParser(usage='usage: %prog -t <base_directory> -s <spectra_directory -o <output_root_file> -p <plots_directory>')
    parser.add_option('-t', '--directory', dest='directory', default=BASE_DIR, help='Base directory to search for ROOT files (default: %s)' % directory)
    parser.add_option('-s', '--spectra', dest='spectra_directory', default=SPECTRA_DIR, help='Directory to load spectra tables from (default: %s)' % SPECTRA_DIR)
    parser.add_option('-o', '--output', dest='outfile', default=OUTFILE, help='Output ROOT file path (default: %s)' % OUTFILE)
    parser.add_option('-p', '--plots', dest='plots_directory', default=PLOTS_DIR, help='Directory to save plots (default: %s)' % PLOTS_DIR)
    (options, args) = parser.parse_args()

    directory         = options.directory
    spectra_directory = options.spectra_directory
    outfile           = options.outfile
    plots_dir         = options.plots_directory

    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
    
    files = find_root_files(directory)
    if not files:
        log.error("No ROOT files found under %s", directory)
        return

    spectra_tables = load_spectra_table(spectra_directory)

    # Plot all spectra
    plot_all_flux_spectra(spectra_tables)

    # Collect per-source hitProb for plotting
    per_source: Dict[str, Dict[str, List[float]]] = {}

    # Output columns
    out_source: List[str] = []
    out_particle: List[str] = []
    out_filename: List[str] = []
    out_energy: List[float] = []
    out_nbeam: List[int] = []
    out_numhits: List[int] = []
    out_hitprob: List[float] = []
    out_err_hitprob: List[float] = []
    out_flux: List[float] = []
    out_rate: List[float] = []
    out_err_rate: List[float] = []
    out_edep_lists: List[List[float]] = []

    log.info("Processing %d files...", len(files))
    for i, info in enumerate(files, 1):
        try:
            data = read_file_data(info.path)

            # fluxSpace from spectra table matching the "source"
            if info.source in spectra_tables and np.isfinite(info.energy_MeV):
                Egrid, Fgrid = spectra_tables[info.source]
                flux_space = interp_flux_loglog(Egrid, Fgrid, float(info.energy_MeV))
            else:
                flux_space = float("nan")

            rate, err_rate = compute_rate(flux_space, SPHERE_RADIUS_MM, data.hit_prob, data.err_hit_prob)

            # Fill lists
            out_source.append(info.source)
            out_particle.append(info.particle)
            out_filename.append(info.filename)
            out_energy.append(float(info.energy_MeV))
            out_nbeam.append(int(data.n_beam))
            out_numhits.append(int(data.num_hits))
            out_hitprob.append(float(data.hit_prob))
            out_err_hitprob.append(float(data.err_hit_prob))
            out_flux.append(float(flux_space))
            out_rate.append(float(rate))
            out_err_rate.append(float(err_rate))
            out_edep_lists.append([float(x) for x in data.edep_keV])

            # per-source for plot
            per_source.setdefault(info.source, {"energy": [], "hitProb": [], "err_hitProb": []})
            per_source[info.source]["energy"].append(float(info.energy_MeV))
            per_source[info.source]["hitProb"].append(float(data.hit_prob))
            per_source[info.source]["err_hitProb"].append(float(data.err_hit_prob))

            if i % 50 == 0:
                log.info("  %d/%d done", i, len(files))

        except Exception as e:
            log.error("Failed on %s: %s", info.filename, e)

    # Plots using collected values
    plot_hitprob_vs_energy(per_source)

    # Prepare branches for writing
    if ak is None:
        log.error(
            "awkward is not available, so I can't write the jagged edepGasTotal_keV branch safely.\n"
            "Install awkward (usually already with uproot) or tell me and I will write (edep_flat, edep_offsets) instead."
        )
        return

    # --- compute flux weights per entry (one per file) ---
    dE_MeV = compute_dE_per_entry(out_source, out_energy)

    flux_arr  = np.array(out_flux, dtype=np.float64)
    nbeam_arr = np.array(out_nbeam, dtype=np.float64)

    # this is the energy-bin integrated flux weight for that simulated point
    w_flux_dE = flux_arr * dE_MeV

    # per-hit event weight: each hit in that file gets the same weight
    # (so that summing weights of hits reproduces flux*hitProb*dE)
    w_perHit = np.zeros_like(w_flux_dE)
    ok = np.isfinite(w_flux_dE) & (nbeam_arr > 0)
    w_perHit[ok] = w_flux_dE[ok] / nbeam_arr[ok]
    plot_edep_histograms_weighted(
        sorted(per_source.keys()),
        out_source,
        out_edep_lists,
        out_numhits,
        w_perHit,
        plots_dir,
    )


    rows = {
        "source": ak.Array(out_source),
        "particle": ak.Array(out_particle),
        "filename": ak.Array(out_filename),
        "energy_MeV": np.array(out_energy, dtype=np.float64),
        "nBeamOnRequested": np.array(out_nbeam, dtype=np.int64),
        "numHits": np.array(out_numhits, dtype=np.int64),
        "hitProb": np.array(out_hitprob, dtype=np.float64),
        "err_hitProb": np.array(out_err_hitprob, dtype=np.float64),
        "fluxSpace": np.array(out_flux, dtype=np.float64),
        "dE_MeV": dE_MeV,
        "w_flux_dE": w_flux_dE,
        "w_perHit": w_perHit,
        "rate_Hz": np.array(out_rate, dtype=np.float64),
        "err_rate_Hz": np.array(out_err_rate, dtype=np.float64),
        "edepGasTotal_keV": ak.Array(out_edep_lists),
    }

    write_summary_root(outfile, rows)
    log.info("Wrote %s with TTree 'FileSummary' (%d entries)", outfile, len(out_source))


if __name__ == "__main__":
    main()
