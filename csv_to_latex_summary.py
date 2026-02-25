#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


CASE_ORDER_DEFAULT = ["standard", "timepix_cut", "timepix_diffusion"]
CASE_LABELS_DEFAULT = {
    "standard": "Standard",
    "timepix_cut": "Timepix cut",
    "timepix_diffusion": "Timepix+diffusion",
}


def _round_to_sig(x: float, sig: int) -> float:
    """Round x to `sig` significant digits."""
    if x == 0 or not math.isfinite(x):
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))

def _decimal_places(x: float) -> int:
    """
    Number of decimal places needed to represent x (already rounded) without scientific notation,
    based on its order of magnitude.
    """
    if x == 0 or not math.isfinite(x):
        return 0
    # e.g. x=0.012 -> log10=-1.92 -> floor=-2 -> need 2 decimals? Actually 0.01 is 2 dp
    # decimals = max(0, -floor(log10(|x|)))
    return max(0, -int(math.floor(math.log10(abs(x)))))

# def fmt_pm(val: float, err: float, err_sig: int = 2) -> str:
#     """
#     Format as LaTeX: value ± error
#     where error has `err_sig` significant digits (default 2),
#     and value is rounded to the same decimal place as the rounded error.
#     """
#     if val is None or err is None:
#         return r"--"
#     if not (math.isfinite(val) and math.isfinite(err)):
#         return r"--"
#     if err < 0:
#         err = abs(err)
#     if err == 0:
#         return rf"${val:g} \pm 0$"

#     # 1) round uncertainty to 2 significant digits
#     err_r = _round_to_sig(err, err_sig)

#     # 2) round value to the same decimal place as err_r
#     dp = _decimal_places(err_r)
#     val_r = round(val, dp)

#     # 3) choose fixed-point vs scientific notation
#     # Use scientific if numbers are too small/large (optional heuristic)
#     scale = max(abs(val_r), abs(err_r))
#     if scale != 0 and (scale < 1e-3 or scale >= 1e4):
#         # scientific notation, same exponent for both
#         exp = int(math.floor(math.log10(scale)))
#         v = val_r / (10 ** exp)
#         e = err_r / (10 ** exp)

#         # keep dp in mantissa consistent with rounding from err
#         # (dp refers to absolute decimals; convert to mantissa decimals)
#         # mantissa decimals ~ max(0, dp - exp)
#         mant_dp = max(0, dp - exp)
#         v_str = f"{v:.{mant_dp}f}"
#         e_str = f"{e:.{mant_dp}f}"
#         return rf"${v_str}\pm {e_str}\times 10^{{{exp}}}$"

#     # fixed-point
#     v_str = f"{val_r:.{dp}f}"
#     e_str = f"{err_r:.{dp}f}"
#     return rf"${v_str} \pm {e_str}$"

def fmt_pm(val: float, err: float, err_sig: int = 2) -> str:
    """
    LaTeX wrapper around format_avg_std_sci:
    returns: $(a \pm s)\times 10^{e}$
    with uncertainty at `err_sig` significant digits.
    """
    try:
        s = format_avg_std_sci(val, err, unc_sig_digits=err_sig)
    except Exception:
        return r"--"

    # Make it LaTeX-friendly
    s = s.replace("×", r"\times")
    return f"${s}$"

def format_avg_std_sci(avg: float, std: float, *, unc_sig_digits: int = 2) -> str:
    """
    Always scientific notation:
      (a ± s) × 10^e
    where:
      - e is chosen from avg (so the mantissa is in [1, 10) for nonzero avg)
      - uncertainty has `unc_sig_digits` significant digits
      - avg is rounded to the same decimal place as the uncertainty (in mantissa space)

    Examples:
      (1234.5, 67.89)  -> "(1.235 ± 0.068) × 10^3"
      (0.001234, 5.6e-5)-> "(1.234 ± 0.056) × 10^-3"
      (-2.345, 0.01234) -> "(-2.345 ± 0.012) × 10^0"
    """
    if not (math.isfinite(avg) and math.isfinite(std)):
        raise ValueError("avg and std must be finite numbers")
    if std < 0:
        raise ValueError("std must be non-negative")
    if std == 0:
        # Still force scientific notation based on avg (or 0 if avg==0)
        e = 0 if avg == 0 else math.floor(math.log10(abs(avg)))
        m = avg / (10 ** e) if avg != 0 else 0.0
        return f"({m} ± 0) × 10^{e}"

    # Choose exponent from avg if possible; otherwise from std (avg==0)
    if avg != 0:
        e = math.floor(math.log10(abs(avg)))
    else:
        e = math.floor(math.log10(std))

    scale = 10 ** e
    m_avg = avg / scale
    m_std = std / scale  # uncertainty in mantissa units

    # Determine decimals so mantissa-uncertainty has `unc_sig_digits` significant digits
    exp_u = math.floor(math.log10(m_std))  # can be negative
    decimals = max(0, unc_sig_digits - 1 - exp_u)

    m_std_r = round(m_std, decimals)

    # Handle rounding bump (e.g. 9.95 -> 10.0) for m_std
    if m_std_r != 0:
        new_exp_u = math.floor(math.log10(abs(m_std_r)))
        if new_exp_u != exp_u:
            exp_u = new_exp_u
            decimals = max(0, unc_sig_digits - 1 - exp_u)
            m_std_r = round(m_std, decimals)

    m_avg_r = round(m_avg, decimals)

    # Fixed-point formatting for consistent decimals
    if decimals > 0:
        fmt = f"{{:.{decimals}f}}"
        return f"({fmt.format(m_avg_r)} \\pm {fmt.format(m_std_r)}) × 10^{e}"
    else:
        return f"({int(round(m_avg_r))} \\pm {int(round(m_std_r))}) × 10^{e}"


def make_latex_table(
    df: pd.DataFrame,
    case_order: list[str],
    case_labels: dict[str, str],
    sig: int = 3,
    include_total_first: bool = True,
    caption: str | None = None,
    label: str | None = None,
) -> str:
    # Keep only needed columns
    df = df[["case", "source", "integrated_rate_1_over_s", "sigma_1_over_s"]].copy()

    # Pivot into columns per case
    pivot_val = df.pivot(index="source", columns="case", values="integrated_rate_1_over_s")
    pivot_err = df.pivot(index="source", columns="case", values="sigma_1_over_s")

    # Ensure column order and presence
    for c in case_order:
        if c not in pivot_val.columns:
            pivot_val[c] = float("nan")
            pivot_err[c] = float("nan")

    pivot_val = pivot_val[case_order]
    pivot_err = pivot_err[case_order]

    # Build formatted table rows
    sources = list(pivot_val.index)

    # Optional ordering: Total first, then alphabetical
    if include_total_first and "Total" in sources:
        sources_sorted = ["Total"] + sorted([s for s in sources if s != "Total"])
    else:
        sources_sorted = sorted(sources)

    rows = []
    for src in sources_sorted:
        cells = [src]
        for c in case_order:
            cells.append(fmt_pm(pivot_val.loc[src, c], pivot_err.loc[src, c], err_sig=sig))
        rows.append(cells)

    # LaTeX table string (booktabs style)
    header_cols = ["Source"] + [case_labels.get(c, c) for c in case_order]
    col_spec = "l" + "c" * len(case_order)

    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(rf"\begin{{tabular}}{{{col_spec}}}")
    lines.append(r"\toprule")
    lines.append(" & ".join(header_cols) + r" \\")
    lines.append(r"\midrule")

    for r in rows:
        lines.append(" & ".join(r) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")

    if caption:
        lines.append(rf"\caption{{{caption}}}")
    if label:
        lines.append(rf"\label{{{label}}}")

    lines.append(r"\end{table}")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Print a LaTeX summary table from integrated_rates CSV.")
    p.add_argument("csv", type=Path, help="Path to integrated_rates CSV.")
    p.add_argument("--sig", type=int, default=2, help="Significant digits for formatting (default: 3).")
    p.add_argument("--out", type=Path, default=None, help="Write LaTeX to this file instead of stdout.")
    p.add_argument("--caption", type=str, default="Integrated rates per source.", help="LaTeX caption.")
    p.add_argument("--label", type=str, default="tab:integrated_rates", help="LaTeX label.")
    p.add_argument("--no-total-first", action="store_true", help="Do not force 'Total' row first.")
    p.add_argument(
        "--cases",
        type=str,
        default=",".join(CASE_ORDER_DEFAULT),
        help="Comma-separated case order (default: standard,timepix_cut,timepix_diffusion)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    df = pd.read_csv(args.csv)

    case_order = [c.strip() for c in args.cases.split(",") if c.strip()]
    latex = make_latex_table(
        df=df,
        case_order=case_order,
        case_labels=CASE_LABELS_DEFAULT,
        sig=args.sig,
        include_total_first=not args.no_total_first,
        caption=args.caption if args.caption else None,
        label=args.label if args.label else None,
    )

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(latex)
    else:
        print(latex)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
