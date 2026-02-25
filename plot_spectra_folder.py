#!/usr/bin/env python3
import csv
from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")  # <-- add this before importing pyplot
import matplotlib.pyplot as plt


SPECTRA_DIR = Path("./spectra")
MANIFEST = SPECTRA_DIR / "manifest.csv"
OUT_PNG = "spectra_overview.png"

def load_manifest(manifest_path: Path):
    items = []
    with open(manifest_path, newline="") as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            fn = row.get("filename", "").strip()
            particle = row.get("particle", "").strip()
            if fn:
                items.append((SPECTRA_DIR / fn, particle))
    if not items:
        raise ValueError(f"No entries found in {manifest_path}")
    return items

def read_spectrum_raw(path: Path):
    """Return raw (E, F, raw_line_no, raw_line_text) for diagnostics."""
    rows = []
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines, start=1):
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
        rows.append((e, f, i, line))
    return rows

def diagnose(rows):
    """Find ordering/duplicate/invalid issues."""
    issues = {
        "nonpositive_E": [],
        "nonpositive_F": [],
        "nonincreasing_pairs": [],   # (idx, E_i, E_next)
        "duplicates": [],            # (idx, E_i)
    }
    # raw order checks
    for (e, f, ln, txt) in rows:
        if e <= 0:
            issues["nonpositive_E"].append((ln, e, txt))
        if f <= 0:
            issues["nonpositive_F"].append((ln, f, txt))

    # order + duplicates (in file order)
    for k in range(len(rows) - 1):
        e1, f1, ln1, t1 = rows[k]
        e2, f2, ln2, t2 = rows[k + 1]
        if e2 < e1:
            issues["nonincreasing_pairs"].append((k, (ln1, e1, t1), (ln2, e2, t2)))
        if e2 == e1:
            issues["duplicates"].append((k, (ln1, e1, t1), (ln2, e2, t2)))

    return issues

def main():
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Missing {MANIFEST}")
    items = load_manifest(MANIFEST)

    fig = plt.figure()
    ax = plt.gca()

    any_plotted = False

    print("\n=== SPECTRA DIAGNOSTICS ===")
    for spec_path, particle in items:
        if not spec_path.exists():
            print(f"[MISSING] {spec_path} (particle={particle})")
            continue

        rows = read_spectrum_raw(spec_path)
        if len(rows) < 2:
            print(f"[TOO FEW POINTS] {spec_path.name} (particle={particle})")
            continue

        issues = diagnose(rows)

        # For plotting: sort by energy, keep only positive E and positive F (needed for log-log)
        rows_sorted = sorted(rows, key=lambda x: x[0])
        E = [e for (e, f, ln, txt) in rows_sorted if e > 0 and f > 0]
        F = [f for (e, f, ln, txt) in rows_sorted if e > 0 and f > 0]

        # Print diagnostics summary
        tag = f"{spec_path.name} ({particle})"
        has_issue = any(len(v) > 0 for v in issues.values())
        if has_issue:
            print(f"\n[ISSUES] {tag}")
            if issues["nonincreasing_pairs"]:
                k, a, b = issues["nonincreasing_pairs"][0]
                print(f"  - Non-increasing energy at file order pair around lines {a[0]} -> {b[0]}: {a[1]} -> {b[1]}")
                print(f"    line {a[0]}: {a[2]}")
                print(f"    line {b[0]}: {b[2]}")
            if issues["duplicates"]:
                k, a, b = issues["duplicates"][0]
                print(f"  - Duplicate energy at lines {a[0]} and {b[0]}: {a[1]}")
                print(f"    line {a[0]}: {a[2]}")
                print(f"    line {b[0]}: {b[2]}")
            if issues["nonpositive_E"]:
                ln, e, txt = issues["nonpositive_E"][0]
                print(f"  - Non-positive E at line {ln}: E={e}")
                print(f"    {txt}")
            if issues["nonpositive_F"]:
                ln, f, txt = issues["nonpositive_F"][0]
                print(f"  - Non-positive flux at line {ln}: F={f}")
                print(f"    {txt}")
        else:
            print(f"[OK] {tag}")

        # Plot (if we have data after filtering)
        if len(E) >= 2:
            ax.loglog(E, F, label=f"{spec_path.stem} [{particle}]")
            any_plotted = True
        else:
            print(f"  -> Not plotted (insufficient positive points for log-log).")

    print("\n==========================\n")

    if not any_plotted:
        print("No spectra plotted. Check manifest paths and file contents.")
        return

    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel(r"Flux (MeV$^{-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$)")
    ax.set_title("Background spectra overview (log-log)")
    ax.grid(True, which="both")
    ax.legend(fontsize=7)

    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=200)
    print(f"Saved: {OUT_PNG}")
    # no plt.show() in headless mode


if __name__ == "__main__":
    main()
