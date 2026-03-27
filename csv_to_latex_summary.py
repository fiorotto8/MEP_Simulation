#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


CASE_LABELS_DEFAULT = {
    "standard": "Standard",
    "timepix_cut": "Timepix cut",
    "timepix_diffusion": "Timepix+diffusion",
    "timepix_diffusion_gagg": "Timepix+diffusion+GAGG veto",
    "timepix_diffusion_gaggIN": "Timepix+diffusion+GAGG IN",
}

# Hard-coded detector area in cm^2
# Change here to whatever active area you want to use.
AREA_CM2 = 1.408 * 1.408


def humanize_case_name(case: str) -> str:
    return case.replace("_", " ")


def fmt_pm(
    val: float,
    err: float,
    digits: int = 2,
    sci_large: float = 1e4,
) -> str:
    """
    Format as LaTeX: value ± error.

    Rule:
    - fixed notation if the number can be represented with at most `digits`
      decimals;
    - scientific notation if it would require more than `digits` decimals;
    - scientific notation also for very large numbers.
    """
    if val is None or err is None:
        return r"--"
    if not (math.isfinite(val) and math.isfinite(err)):
        return r"--"

    err = abs(err)

    if val == 0 and err == 0:
        z = f"{0:.{digits}f}"
        return rf"${z} \pm {z}$"

    scale = max(abs(val), abs(err))
    use_sci = False

    if scale != 0:
        exp = int(math.floor(math.log10(scale)))
        if exp < -digits or scale >= sci_large:
            use_sci = True

    if use_sci:
        exp = 0 if scale == 0 else int(math.floor(math.log10(scale)))
        v = val / (10 ** exp)
        e = err / (10 ** exp)
        return rf"$({v:.{digits}f} \pm {e:.{digits}f})\times 10^{{{exp}}}$"

    return rf"${val:.{digits}f} \pm {err:.{digits}f}$"


def discover_case_order(df: pd.DataFrame) -> list[str]:
    preferred = [
        "standard",
        "timepix_cut",
        "timepix_diffusion",
        "timepix_diffusion_gagg",
        "timepix_diffusion_gaggIN",
    ]
    found = list(pd.unique(df["case"]))
    ordered = [c for c in preferred if c in found]
    ordered += sorted([c for c in found if c not in ordered])
    return ordered


def build_header_cols(case_order: list[str], case_labels: dict[str, str], unit_label: str) -> list[str]:
    return ["Source"] + [f"Rate {case_labels.get(c, humanize_case_name(c))} ({unit_label})" for c in case_order]


def make_latex_table(
    df: pd.DataFrame,
    case_order: list[str],
    case_labels: dict[str, str],
    digits: int = 2,
    include_total_first: bool = True,
    caption: str | None = None,
    label: str | None = None,
    unit_label: str = r"Hz",
) -> str:
    df = df[["case", "source", "integrated_rate_1_over_s", "sigma_1_over_s"]].copy()

    pivot_val = df.pivot(index="source", columns="case", values="integrated_rate_1_over_s")
    pivot_err = df.pivot(index="source", columns="case", values="sigma_1_over_s")

    for c in case_order:
        if c not in pivot_val.columns:
            pivot_val[c] = float("nan")
            pivot_err[c] = float("nan")

    pivot_val = pivot_val[case_order]
    pivot_err = pivot_err[case_order]

    sources = list(pivot_val.index)
    if include_total_first and "Total" in sources:
        sources_sorted = ["Total"] + sorted([s for s in sources if s != "Total"])
    else:
        sources_sorted = sorted(sources)

    rows = []
    for src in sources_sorted:
        cells = [src]
        for c in case_order:
            cells.append(fmt_pm(pivot_val.loc[src, c], pivot_err.loc[src, c], digits=digits))
        rows.append(cells)

    header_cols = build_header_cols(case_order, case_labels, unit_label)
    col_spec = "l" + "c" * len(case_order)

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\resizebox{\textwidth}{!}{%")
    lines.append(rf"\begin{{tabular}}{{{col_spec}}}")
    lines.append(r"\toprule")
    lines.append(" & ".join(header_cols) + r" \\")
    lines.append(r"\midrule")

    for r in rows:
        lines.append(" & ".join(r) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}%")
    lines.append(r"}")

    if caption:
        lines.append(rf"\caption{{{caption}}}")
    if label:
        lines.append(rf"\label{{{label}}}")

    lines.append(r"\end{table}")
    return "\n".join(lines)


def normalize_rate_density(df: pd.DataFrame, area_cm2: float) -> pd.DataFrame:
    """
    Convert integrated_rate_1_over_s and sigma_1_over_s
    from Hz to Hz/cm^2/keV using:
        rate_density = rate / (area_cm2 * (Emax - Emin))
    """
    out = df.copy()

    if "energy_min_keV" not in out.columns or "energy_max_keV" not in out.columns:
        raise KeyError("CSV must contain energy_min_keV and energy_max_keV columns")

    delta_e = out["energy_max_keV"].astype(float) - out["energy_min_keV"].astype(float)
    if (delta_e <= 0).any():
        bad = out.loc[delta_e <= 0, ["case", "source", "energy_min_keV", "energy_max_keV"]]
        raise ValueError(f"Found non-positive energy range in rows:\n{bad}")

    norm = area_cm2 * delta_e

    out["integrated_rate_1_over_s"] = out["integrated_rate_1_over_s"].astype(float) / norm
    out["sigma_1_over_s"] = out["sigma_1_over_s"].astype(float) / norm
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Print LaTeX summary tables from integrated_rates CSV.")
    p.add_argument("csv", type=Path, help="Path to integrated_rates CSV.")
    p.add_argument("--digits", type=int, default=4, help="Digits after the decimal point in fixed notation or mantissa.")
    p.add_argument("--out", type=Path, default=None, help="Write the Hz table to this file instead of stdout.")
    p.add_argument("--out-density", type=Path, default=None, help="Write the Hz/cm^2/keV table to this file.")
    p.add_argument("--caption", type=str, default="Integrated rates per source.", help="Caption for the Hz table.")
    p.add_argument(
        "--caption-density",
        type=str,
        default="Integrated rates per source normalized to detector area and energy interval.",
        help="Caption for the Hz/cm^2/keV table.",
    )
    p.add_argument("--label", type=str, default="tab:integrated_rates", help="LaTeX label for the Hz table.")
    p.add_argument(
        "--label-density",
        type=str,
        default="tab:integrated_rate_density",
        help="LaTeX label for the Hz/cm^2/keV table.",
    )
    p.add_argument("--no-total-first", action="store_true", help="Do not force 'Total' row first.")
    p.add_argument(
        "--cases",
        type=str,
        default=None,
        help="Comma-separated case order. If omitted, all distinct cases in the CSV are used.",
    )
    p.add_argument(
        "--area-cm2",
        type=float,
        default=AREA_CM2,
        help=f"Detector area in cm^2 for normalized table (default: {AREA_CM2:.6g}).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    df = pd.read_csv(args.csv)

    if args.cases:
        case_order = [c.strip() for c in args.cases.split(",") if c.strip()]
    else:
        case_order = discover_case_order(df)

    case_labels = dict(CASE_LABELS_DEFAULT)
    for c in case_order:
        if c not in case_labels:
            case_labels[c] = humanize_case_name(c)

    latex_hz = make_latex_table(
        df=df,
        case_order=case_order,
        case_labels=case_labels,
        digits=args.digits,
        include_total_first=not args.no_total_first,
        caption=args.caption if args.caption else None,
        label=args.label if args.label else None,
        unit_label=r"Hz",
    )

    df_density = normalize_rate_density(df, area_cm2=args.area_cm2)

    latex_density = make_latex_table(
        df=df_density,
        case_order=case_order,
        case_labels=case_labels,
        digits=args.digits,
        include_total_first=not args.no_total_first,
        caption=args.caption_density if args.caption_density else None,
        label=args.label_density if args.label_density else None,
        unit_label=r"Hz/cm$^2$/keV",
    )

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(latex_hz)
    else:
        print(latex_hz)
        print()
        print(latex_density)

    if args.out_density:
        args.out_density.parent.mkdir(parents=True, exist_ok=True)
        args.out_density.write_text(latex_density)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())