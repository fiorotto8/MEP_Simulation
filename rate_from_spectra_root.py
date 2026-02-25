#!/usr/bin/env python3
import os
import re
import math
import csv
import json
import time
import subprocess
import tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Pattern

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- SETTINGS ----------
EXE = "./build/gpd3d"
SPECTRA_DIR = Path("./spectra")
MANIFEST = SPECTRA_DIR / "manifest.csv"
TEMPLATE_MAC = Path("macros/template_bg_sphere.mac")

SPHERE_RADIUS_MM = 130.0
#N_EVENTS = 1000000
THINNING = 1

MAX_WORKERS = None  # None -> os.cpu_count()

# Output
OUTDIR = Path("run_outputs")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_ROOTDIR = OUTDIR / "root"
OUT_ROOTDIR.mkdir(parents=True, exist_ok=True)

OUT_RUNCSV = OUTDIR / "run_points.csv"
OUT_RUNJSON = OUTDIR / "run_points.json"
OUT_SUMCSV = OUTDIR / "run_summary.csv"
OUT_PNG_RATE = OUTDIR / "rates_overview.png"
OUT_PNG_PHIT = OUTDIR / "phit_overview.png"

# ---------- PARSING ----------
RE_P_HIT = re.compile(r"P\(hit\):\s*([0-9eE+\-\.]+)")
RE_NGEN  = re.compile(r"Generated events(?:[^:]*)?:\s*(\d+)")
RE_NHIT  = re.compile(r"Hit events(?:[^:]*)?:\s*(\d+)")


# ---------- MACRO INJECTION HELPERS ----------
RE_SPHERE_RADIUS = re.compile(r"/gpd3d/gen/sphereRadius\s+[0-9eE+\-\.]+\s+mm")
RE_BEAMON       = re.compile(r"/run/beamOn\s+\d+")
RE_WRITE_ROOT = re.compile(r"^\s*/gpd3d/bg/writeROOT\s+\w+.*$", re.IGNORECASE | re.MULTILINE)
RE_ROOT_FILE  = re.compile(r"^\s*/gpd3d/bg/rootFile\s+.*$", re.MULTILINE)

def sanitize_tag(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_\-\.]+", "_", s)
    return s.strip("_")

def energy_tag_mev(E: float) -> str:
    # stable, filesystem-safe tag
    return sanitize_tag(f"E{E:.6g}MeV")

def ensure_line(mac_text: str, key_regex: Pattern, line: str) -> str:
    """Replace existing line matching key_regex, else insert near the top."""
    if key_regex.search(mac_text):
        return key_regex.sub(line, mac_text)
    # insert after verbosity block (or at top)
    lines = mac_text.splitlines(True)
    insert_at = 0
    for i, ln in enumerate(lines):
        if ln.strip().startswith("/tracking/verbose"):
            insert_at = i + 1
            break
    lines.insert(insert_at, line + "\n")
    return "".join(lines)

def inject_gps_block(mac_text: str, particle: str, energy_mev: float) -> str:
    inject = (
        f"/gps/particle {particle}\n"
        f"/gps/ene/type Mono\n"
        f"/gps/ene/mono {energy_mev} MeV\n"
    )
    # Insert right before /run/initialize (as your original)
    return mac_text.replace("/run/initialize", inject + "\n/run/initialize", 1)

def build_root_path(root_dir: Path, spec_path: Path, particle: str, energy_mev: float) -> Path:
    # root/<particle>/<spectrum_stem>/<energy_tag>.root
    part_dir = root_dir / sanitize_tag(particle)
    spec_dir = part_dir / sanitize_tag(spec_path.stem)
    spec_dir.mkdir(parents=True, exist_ok=True)
    fn = f"{sanitize_tag(spec_path.stem)}_{sanitize_tag(particle)}_{energy_tag_mev(energy_mev)}.root"
    return spec_dir / fn

# ---------- SPECTRA ----------
def read_spectrum_csv(path: Path):
    pairs = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = [p.strip() for p in s.split(",")]
        if len(parts) < 2:
            continue
        try:
            e = float(parts[0])
            f = float(parts[1])
        except ValueError:
            continue
        if e <= 0:
            continue
        pairs.append((e, f))
    if len(pairs) < 2:
        raise ValueError(f"Need >=2 valid points in {path}")
    pairs.sort(key=lambda x: x[0])

    E, F = [], []
    for e, f in pairs:
        if E and abs(e - E[-1]) <= 0.0:
            F[-1] += f
        else:
            E.append(e)
            F.append(f)
    if len(E) < 2:
        raise ValueError(f"Need >=2 unique energies in {path}")
    return E, F

def infer_edges_from_centers(E):
    n = len(E)
    edges = [0.0] * (n + 1)
    for i in range(n - 1):
        edges[i + 1] = math.sqrt(E[i] * E[i + 1])
    edges[0] = E[0] * E[0] / edges[1]
    edges[n] = E[n - 1] * E[n - 1] / edges[n - 1]
    for i in range(n):
        if edges[i+1] <= edges[i]:
            raise ValueError(f"Non-increasing inferred edges at i={i}")
    return edges

def load_manifest(manifest_path: Path):
    items = []
    with open(manifest_path, newline="") as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            fn = row["filename"].strip()
            particle = row["particle"].strip()
            n_events = int(row["n_events"])
            if fn:
                items.append((SPECTRA_DIR / fn, particle, n_events))
    if not items:
        raise ValueError(f"No entries found in {manifest_path}")
    return items


# ---------- RUN ONE POINT ----------
def run_one_point(exe: str, template_mac_text: str, spec_path: Path,
                  particle: str, energy_mev: float, sphere_radius_mm: float,
                  n_events: int, root_out: Path):
    mac_text = template_mac_text

    # enforce background ROOT output + file name
    mac_text = ensure_line(mac_text, RE_WRITE_ROOT, f"/gpd3d/bg/writeROOT true")
    mac_text = ensure_line(mac_text, RE_ROOT_FILE,  f"/gpd3d/bg/rootFile {root_out.as_posix()}")

    # override radius + beamOn
    mac_text = RE_SPHERE_RADIUS.sub(f"/gpd3d/gen/sphereRadius {sphere_radius_mm} mm", mac_text)
    mac_text = RE_BEAMON.sub(f"/run/beamOn {n_events}", mac_text)

    # inject GPS block
    mac_text = inject_gps_block(mac_text, particle, energy_mev)

    # (optional) keep one thread per process if your G4 is MT-enabled
    # mac_text = ensure_line(mac_text, re.compile(r"/run/numberOfThreads\s+\d+"), "/run/numberOfThreads 1")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".mac",
                                     prefix=f"gpd3d_{sanitize_tag(particle)}_",
                                     delete=False) as tf:
        tf.write(mac_text)
        tmp_mac = tf.name

    t0 = time.time()
    proc = subprocess.run(
        [exe, tmp_mac],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True
    )
    dt = time.time() - t0
    out = proc.stdout

    #try:
    #    os.remove(tmp_mac)
    #except OSError:
    #    pass

    if proc.returncode != 0:
        raise RuntimeError(f"Run failed for {particle} {energy_mev} MeV\n\n{out[-2000:]}")

    mP = RE_P_HIT.search(out)
    mNgen = RE_NGEN.search(out)
    mNhit = RE_NHIT.search(out)
    if not (mP and mNgen and mNhit):
        raise RuntimeError(f"Could not parse background summary. Output tail:\n\n{out[-2000:]}")

    phit = float(mP.group(1))
    ngen = int(mNgen.group(1))
    nhit = int(mNhit.group(1))

    # If you implement "saveOnlyHitEvents", you may still end up with an empty ROOT file
    # for phit=0. Keep it anyway; it is useful bookkeeping.
    return phit, ngen, nhit, dt, str(root_out)

# ---------- OUTPUT ----------
def save_points_csv(rows, path: Path):
    fieldnames = [
        "source_file","particle","index",
        "E_MeV","flux_differential","dE_MeV","Phi_bin",
        "P_hit","N_gen","N_hit","A_gen_cm2","A_eff_cm2",
        "R_i_per_s","runtime_s",
        "root_file"
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def save_summary_csv(rows, path: Path):
    fieldnames = [
        "source_file","particle","n_points_used",
        "rate_per_s","sum_Phi_bin","weighted_Phit"
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def make_overview_plots(points_rows):
    groups = {}
    for r in points_rows:
        key = (r["source_file"], r["particle"])
        groups.setdefault(key, []).append(r)

    fig1 = plt.figure()
    ax1 = plt.gca()
    any1 = False
    for (src, part), rows in groups.items():
        rows = sorted(rows, key=lambda x: x["E_MeV"])
        E = [x["E_MeV"] for x in rows if x["R_i_per_s"] > 0 and x["E_MeV"] > 0]
        R = [x["R_i_per_s"] for x in rows if x["R_i_per_s"] > 0 and x["E_MeV"] > 0]
        if len(E) >= 2:
            ax1.loglog(E, R, label=f"{Path(src).stem} [{part}]")
            any1 = True
    ax1.set_xlabel("Primary particle Energy (MeV)")
    ax1.set_ylabel("Rate contribution per point Ri (1/s)")
    ax1.set_title("Rate contributions vs energy (log-log)")
    ax1.grid(True, which="both")
    if any1:
        ax1.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(str(OUT_PNG_RATE), dpi=200)
    plt.close(fig1)

    fig2 = plt.figure()
    ax2 = plt.gca()
    any2 = False
    for (src, part), rows in groups.items():
        rows = sorted(rows, key=lambda x: x["E_MeV"])
        E = [x["E_MeV"] for x in rows if x["E_MeV"] > 0]
        P = [x["P_hit"] for x in rows if x["E_MeV"] > 0]
        if len(E) >= 2:
            ax2.semilogx(E, P, label=f"{Path(src).stem} [{part}]")
            any2 = True
    ax2.set_xlabel("Primary particle Energy (MeV)")
    ax2.set_ylabel("P(hit)  (EdepGas>0)")
    ax2.set_title("Hit probability vs energy")
    ax2.grid(True, which="both")
    if any2:
        ax2.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(str(OUT_PNG_PHIT), dpi=200)
    plt.close(fig2)

# ---------- WORKER ----------
def _worker_run_manifest_item(args):
    (exe, template_mac_text, spec_path_str, particle, sphere_radius_mm, n_events, thinning, out_rootdir_str) = args
    spec_path = Path(spec_path_str)
    out_rootdir = Path(out_rootdir_str)

    R_cm = sphere_radius_mm / 10.0
    A_gen = 4.0 * math.pi * (R_cm ** 2)

    E, F = read_spectrum_csv(spec_path)
    edges = infer_edges_from_centers(E)

    rate_file = 0.0
    sum_phi_bin = 0.0
    sum_phi_bin_phit = 0.0
    n_used = 0
    points = []

    for i in range(0, len(E), thinning):
        Ei = E[i]
        Fi = F[i]
        dE = edges[i+1] - edges[i]
        Phi_bin = Fi * dE

        root_out = build_root_path(out_rootdir, spec_path, particle, Ei)

        phit, ngen, nhit, dt, root_file = run_one_point(
            exe, template_mac_text, spec_path, particle, Ei,
            sphere_radius_mm, n_events, root_out
        )

        A_eff = A_gen * phit
        R_i = math.pi * A_eff * Phi_bin

        rate_file += R_i
        sum_phi_bin += Phi_bin
        sum_phi_bin_phit += Phi_bin * phit
        n_used += 1

        points.append({
            "N_events_cfg": n_events,
            "source_file": spec_path.name,
            "particle": particle,
            "index": i,
            "E_MeV": Ei,
            "flux_differential": Fi,
            "dE_MeV": dE,
            "Phi_bin": Phi_bin,
            "P_hit": phit,
            "N_gen": ngen,
            "N_hit": nhit,
            "A_gen_cm2": A_gen,
            "A_eff_cm2": A_eff,
            "R_i_per_s": R_i,
            "runtime_s": dt,
            "root_file": root_file
        })

    weighted_phit = (sum_phi_bin_phit / sum_phi_bin) if sum_phi_bin > 0 else 0.0
    summary = {
        "n_events_cfg": n_events,
        "source_file": spec_path.name,
        "particle": particle,
        "n_points_used": n_used,
        "rate_per_s": rate_file,
        "sum_Phi_bin": sum_phi_bin,
        "weighted_Phit": weighted_phit
    }
    return summary, points

def main():
    if not TEMPLATE_MAC.exists():
        raise FileNotFoundError(f"Missing {TEMPLATE_MAC}")
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Missing {MANIFEST}")
    if not Path(EXE).exists():
        raise FileNotFoundError(f"Missing executable {EXE}")

    manifest_items = load_manifest(MANIFEST)
    template_mac_text = TEMPLATE_MAC.read_text()

    for spec_path, _, _ in manifest_items:
        if not spec_path.exists():
            raise FileNotFoundError(f"Missing spectrum file {spec_path}")


    worker_args = []
    for spec_path, particle, n_events in manifest_items:
        worker_args.append((
            EXE,
            template_mac_text,
            str(spec_path),
            particle,
            SPHERE_RADIUS_MM,
            n_events,
            THINNING,
            str(OUT_ROOTDIR)
        ))


    results = []
    all_points = []
    total_rate = 0.0

    print(f"Running {len(worker_args)} manifest entries in parallel "
          f"(max_workers={MAX_WORKERS or os.cpu_count()})")
    print(f"ROOT outputs under: {OUT_ROOTDIR}")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        fut_map = {pool.submit(_worker_run_manifest_item, a): a for a in worker_args}
        for fut in as_completed(fut_map):
            a = fut_map[fut]
            spec_path_str = a[2]
            particle = a[3]
            try:
                summary, points = fut.result()
            except Exception as e:
                raise RuntimeError(f"Worker failed for {Path(spec_path_str).name} [{particle}]: {e}") from e

            print(f"\n=== DONE: {summary['source_file']}   particle: {summary['particle']} ===")
            print(f"--> Rate = {summary['rate_per_s']:.6g} 1/s "
                  f"(points={summary['n_points_used']}, weighted P(hit)={summary['weighted_Phit']:.3g})")

            results.append(summary)
            all_points.extend(points)
            total_rate += float(summary["rate_per_s"])

    results.sort(key=lambda r: (r["source_file"], r["particle"]))
    all_points.sort(key=lambda r: (r["source_file"], r["particle"], r["E_MeV"]))

    save_points_csv(all_points, OUT_RUNCSV)
    with open(OUT_RUNJSON, "w") as f:
        json.dump(all_points, f, indent=2)

    save_summary_csv(results, OUT_SUMCSV)
    make_overview_plots(all_points)

    print("\n================ SUMMARY ================")
    for r in results:
        print(f"{r['source_file']:25s}  {r['particle']:7s}  {r['rate_per_s']:.6g} 1/s "
              f"(points={r['n_points_used']}, weighted P(hit)={r['weighted_Phit']:.3g})")
    print("----------------------------------------")
    print(f"TOTAL HIT RATE: {total_rate:.6g} 1/s")
    print("========================================\n")

    print(f"Saved point-by-point: {OUT_RUNCSV}")
    print(f"Saved point-by-point: {OUT_RUNJSON}")
    print(f"Saved summary:        {OUT_SUMCSV}")
    print(f"Saved plots:          {OUT_PNG_RATE}")
    print(f"Saved plots:          {OUT_PNG_PHIT}")
    print(f"Saved ROOT files:     {OUT_ROOTDIR}")

if __name__ == "__main__":
    main()
