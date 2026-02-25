#!/usr/bin/env bash
set -euo pipefail

SCRIPT="analyze_run_spectra.py"
DIRECTORY="run_outputs_continous_spectra_x10/"
PLOT_BASE="${DIRECTORY%/}/plots/"   # e.g. run_outputs_continous_spectra_x10/plots

# ---- parameter grids (edit these) ----
ENERGY_MINS=(3 4 5 6)
ENERGY_MAXS=(30 32 34 36)
BIN_WIDTHS=(0.25 0.5 1 2)
S_SMEARINGS=(0.10 0.15 0.20)
# --------------------------------------

for emin in "${ENERGY_MINS[@]}"; do
  for emax in "${ENERGY_MAXS[@]}"; do
    for bw in "${BIN_WIDTHS[@]}"; do
      for smear in "${S_SMEARINGS[@]}"; do

        echo "Running: emin=$emin emax=$emax bw=$bw smear=$smear"

        python3 "$SCRIPT" \
          --energy-min "$emin" \
          --energy-max "$emax" \
          --bin-width "$bw" \
          --s-smearing "$smear" \
          --directory "$DIRECTORY" \
          --plot-dir "$PLOT_BASE"

      done
    done
  done
done

echo "Done. Results under: $PLOT_BASE"
