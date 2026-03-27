# GPD3D MedEnergy Background Simulation

Geant4-based simulation of the GPD3D medium-energy detector background in LEO.
In background mode, primary particles are generated on a sphere around the detector using the spectra listed in `spectra/manifest.csv`. Event-level and hit-level information are written to ROOT files for later analysis.

## Environment

This project is intended to run inside the provided Docker image.

Build the image:

```bash
docker build -t u1804-root62004-g4_1051 .
```

Run the container:

```bash
xhost +local:

docker run --rm -it \
  -e DISPLAY \
  -e LIBGL_ALWAYS_SOFTWARE=1 \
  -e MESA_LOADER_DRIVER_OVERRIDE=llvmpipe \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -v $HOME/.Xauthority:/root/.Xauthority:ro \
  -v "$(pwd):/work" -w /work \
  u1804-root62004-g4_1051
```

## Build

Inside the container:

```bash
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

Interactive visualization:

```bash
./gpd3d macros/vis.mac
```

## Typical Workflow

1. Generate a simulation run from the spectra manifest:

```bash
python3 launch_run.py --outdir my_run
```

This reads:
- `spectra/manifest.csv`
- `macros/template_bg_sphere_spectrum.mac`

and writes ROOT outputs plus run summaries under the configured output directory. The default output directory is `run_outputs/`.

2. Analyze one run:

```bash
python3 analyze_run.py --directory run_outputs/
```

This produces plots and integrated rate summaries for the selected run.

3. Convert the integrated-rate CSV to LaTeX tables:

```bash
python3 csv_to_latex_summary.py run_outputs/plots/integrated_rates_*.csv
```

## Main Files

- `launch_run.py`: run the Geant4 simulation over all entries in the manifest
- `analyze_run.py`: analyze one simulation output directory
- `csv_to_latex_summary.py`: generate LaTeX summary tables
- `spectra/manifest.csv`: source list, particle type, event count, angular range
- `macros/`: Geant4 macro templates and visualization macros

## Notes

- The repository still contains some legacy scripts under `old_scripts/`.
- Background mode is enabled with `/gpd3d/bg/enable true`.
- The main custom generator options are sphere generation, spectrum sampling from CSV, and angular selection through `posThetaMin` and `posThetaMax`.
