# GPD3D medium-energy background simulation

This repository contains a Geant4 model of the GPD3D medium-energy gas pixel
detector and the scripts used to estimate its particle-induced background.
The main workflow:

1. reads differential particle spectra from CSV files;
2. generates primaries on a sphere surrounding the detector;
3. transports them through the detector, GAGG shield and aluminium enclosure;
4. writes gas energy-deposition events and hits to ROOT files;
5. converts the simulated hit probabilities into physical rates and plots;
6. optionally compares three tungsten geometries using identical inputs and
   matched Geant4 random seeds.

The current production path is the background workflow driven by
`launch_run.py`. Some older polarimetry/FITS code is still present, but it is
not used by this pipeline.

## Simulation model

### Detector and shielding

The geometry is assembled in `src/GPD3D_DetectorConstruction.cc`. Its main
elements are:

- a 1 m vacuum world;
- a layered detector stack with aluminium support, FR4 PCB, Kapton/Vespel,
  copper anode, PEEK drift body, titanium frame and beryllium window;
- a 14.08 mm × 14.08 mm Timepix/foil footprint with SU-8 pillars;
- a 30 mm deep drift region filled with Ar/DME 78/22 at 3 atm and 293.15 K;
- a uniform electric field of 1.85 kV/cm;
- a 10 mm thick cylindrical GAGG wall and a solid GAGG cap on the -Z side;
- a 15 cm aluminium enclosure with 1.5 mm walls and a circular +Z entrance
  hole of 7.5 mm radius;
- an optional 1 mm tungsten plug over the +Z hole;
- an optional 100 mm long tungsten tube, with 7.5 mm inner radius and 0.2 mm
  wall thickness, extending from the plug toward +Z.

The tungsten components are selected at runtime:

| Configuration | `GPD3D_ENABLE_W_PLUG` | `GPD3D_ENABLE_W_TUBE` |
| --- | ---: | ---: |
| plug and tube | `1` | `1` |
| plug only | `1` | `0` |
| neither | `0` | `0` |

Both variables default to enabled when they are unset. A tube requested without
a plug is automatically disabled because the tube is defined as extending from
the plug.

### Physics

The custom physics list combines low-energy polarized electromagnetic models
with decay, hadron elastic, FTFP_BERT, stopping and ion physics. A step limiter
is applied to all particles. The production-cut energy range is 10 eV to
200 GeV; the configured range cuts are 50 µm for photons and 0.2 nm for
electrons.

The repository is tied to the APIs available in Geant4 10.5.1, including the
polarized Livermore photoelectric model used by `GPD3D_Physics.cc`. The supplied
Docker image is therefore the recommended build environment; a substantially
newer host Geant4 installation may not compile this source unchanged.

### Primary generation

The default background generator uses a sphere of radius 130 mm centred on the
detector. Generation points are uniform over the selected polar-angle patch,
and incoming directions follow the cosine law appropriate to an isotropic
intensity crossing a surface.

For every event, the primary energy is sampled from the source CSV using
piecewise log-log interpolation. The spectrum integral and the sampled
simulation are later combined to obtain a physical rate.

## Repository structure

```text
.
├── main.cc                         Geant4 executable entry point
├── CMakeLists.txt                  C++ build configuration
├── Dockerfile                     ROOT 6.20/04 + Geant4 10.5.1 environment
├── include/                        C++ class declarations
├── src/                            geometry, physics, actions and writers
├── macros/                         Geant4 macros and background template
├── spectra/                        spectra, manifest and overview plot
├── launch_run.py                   run every manifest entry
├── analyze_run.py                  ROOT analysis, rates and plots
├── csv_to_latex_summary.py         LaTeX tables from integrated-rate CSV
├── plot_spectra_folder.py          spectrum validation and overview plot
├── scripts/                        three-geometry Docker/Screen launcher
├── debug/                          notebooks and development diagnostics
└── old_scripts/                    legacy scripts, not the production path
```

The central C++ components are:

- `GPD3D_DetectorConstruction`: materials, detector, shields, sensitive gas
  region and electric field;
- `GPD3D_PrimaryGeneratorAction`: sphere/beam generation and spectrum
  sampling;
- `GPD3D_PhysicsList` and `GPD3D_Physics`: particle processes and cuts;
- `GPD3D_GasCellSD` and `GPD3D_GasCellHit`: gas energy-deposition scoring;
- `GPD3D_EventAction` and `GPD3D_SteppingAction`: event summaries,
  primary-containment flags and GAGG tagging;
- `GPD3D_RootWriter`: `RunInfo`, `Events` and `Hits` ROOT trees.

## Requirements

Simulation and compilation:

- Docker with access to the Docker daemon;
- the image built from this repository;
- GNU Screen on the host for the three-geometry launcher.

Analysis:

- Python 3.9 or newer;
- `numpy`, `awkward`, `uproot`, `matplotlib`, `scipy`, `tqdm` and `pandas`.

For example, in a host virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install numpy awkward uproot matplotlib scipy tqdm pandas
```

Ubuntu 18.04 in the supplied image provides the compatible C++ toolchain, ROOT
and Geant4. Its system Python is suitable for `launch_run.py`, but the current
analysis scripts should be run from the newer host environment described
above.

## Build the Docker image

From the repository root:

```bash
docker build -t u1804-root62004-g4_1051 .
```

Docker daemon access is independent of image ownership. If the current user
cannot access `/var/run/docker.sock`, an administrator can grant access once:

```bash
sudo usermod -aG docker "$USER"
```

Then log out and back in. For a temporary refreshed shell, use:

```bash
newgrp docker
```

Rebuilding the image as the current user does not fix Docker socket
permissions.

## Compile the executable

Start a container with the repository mounted at `/work`:

```bash
docker run --rm -it \
  -v "$(pwd):/work" \
  -w /work \
  u1804-root62004-g4_1051
```

Inside the container:

```bash
mkdir -p build
cd build
cmake ..
cmake --build . -- -j"$(nproc)"
```

This creates `build/gpd3d` and copies the Geant4 macros into
`build/macros/`.

### Interactive visualization

On a Linux host with X11, allow the local container and start it with the X11
socket:

```bash
xhost +local:

docker run --rm -it \
  -e DISPLAY \
  -e LIBGL_ALWAYS_SOFTWARE=1 \
  -e MESA_LOADER_DRIVER_OVERRIDE=llvmpipe \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -v "$HOME/.Xauthority:/root/.Xauthority:ro" \
  -v "$(pwd):/work" \
  -w /work \
  u1804-root62004-g4_1051
```

After building, run from the repository root:

```bash
./build/gpd3d macros/vis.mac
```

Running `./build/gpd3d` without arguments also opens the default
`macros/vis.mac`. Supplying one or more macro paths runs them in order; add
`-i` or `--interactive` to keep the UI open afterwards.

## Input spectra and manifest

### Spectrum files

Each source spectrum in `spectra/` is a headerless, two-column CSV:

```text
energy_MeV,differential_flux
```

The flux used by the current workflow has units
MeV^-1 cm^-2 s^-1 sr^-1. Energies and fluxes must be positive, and each file
must contain at least two points. Comment lines beginning with `#` are ignored
by the C++ loader.

The spectrum filename stem is the source identifier. Keep it free of
underscores: `analyze_run.py` recovers the source from the first
underscore-delimited field of the ROOT filename.

Validate and plot all manifest spectra before a long run:

```bash
python3 plot_spectra_folder.py
```

This prints ordering, duplicate and non-positive-value diagnostics and writes
`spectra_overview.png` in the current directory.

### Manifest

`spectra/manifest.csv` controls which simulations are launched:

```csv
filename,particle,n_events,minTheta,maxTheta
CXB.csv,gamma,1000000000,0,114
photonAlbedo.csv,gamma,100000000,114,180
```

The columns are:

| Column | Meaning |
| --- | --- |
| `filename` | CSV filename relative to the spectra directory |
| `particle` | Geant4 particle name, for example `gamma`, `neutron`, `proton`, `alpha`, `e-` or `e+` |
| `n_events` | number of generated primaries for this source |
| `minTheta` | minimum source-position polar angle in degrees |
| `maxTheta` | maximum source-position polar angle in degrees |

`launch_run.py` starts one independent Geant4 process per manifest row, up to
`--max-workers` processes at a time. Choose the event counts with care: the
committed manifest contains production-scale values.

## Run a simulation

Run this command inside the Docker environment after compiling:

```bash
python3 launch_run.py \
  --outdir run_outputs/example \
  --exe ./build/gpd3d \
  --max-workers 4 \
  --base-seed 20260728
```

Important options:

| Option | Default | Purpose |
| --- | --- | --- |
| `--outdir` | `run_outputs` | run output directory |
| `--exe` | `./build/gpd3d` | Geant4 executable |
| `--spectra-dir` | `./spectra` | source spectrum directory |
| `--manifest` | `<spectra-dir>/manifest.csv` | manifest to execute |
| `--template-mac` | `macros/template_bg_sphere_spectrum.mac` | macro patched for each source |
| `--max-workers` | host CPU count | maximum simultaneous Geant4 processes |
| `--base-seed` | unset | deterministic seed from which a pair is derived per manifest row |

Without `--base-seed`, `main.cc` seeds Geant4 from the current time. Reusing
the same positive base seed and the same manifest order gives matched random
streams across geometry configurations.

For a manual geometry selection:

```bash
GPD3D_ENABLE_W_PLUG=1 \
GPD3D_ENABLE_W_TUBE=0 \
python3 launch_run.py \
  --outdir run_outputs/plug_only \
  --exe ./build/gpd3d \
  --max-workers 4 \
  --base-seed 20260728
```

The launcher patches the macro template with the particle, absolute spectrum
path, angular interval, ROOT path, event count, progress interval and optional
random seeds.

### Run output

A completed run has this layout:

```text
run_outputs/example/
├── run_configuration.json
├── run_points.csv
├── run_points.json
├── run_summary.csv
└── root/
    └── <particle>/
        └── <source>/
            └── <source>_<particle>.root
```

`run_configuration.json` records the executable, inputs, worker count, base
seed and tungsten switches. `run_points.*` contains source hit probabilities,
effective areas and rates. `run_summary.csv` contains the per-source summary.

For a source polar-angle interval `[theta_min, theta_max]`, the launcher uses:

```text
A_patch = 2 pi R^2 [cos(theta_min) - cos(theta_max)]
A_eff   = A_patch P(hit)
rate    = pi A_eff integral(Phi(E) dE)
```

Here `R = 13 cm`, `P(hit)` is the fraction of generated events with positive
energy deposition in the sensitive gas, and `Phi` is the differential source
spectrum. The reported launch-time rate is therefore an all-energy gas-hit
rate before the energy window, smearing or GAGG veto used by the later
analysis.

### ROOT schema

Each source ROOT file contains:

- `RunInfo`: run identifier, requested event count, output controls and
  mirrored run configuration;
- `Events`: primary metadata, total/ionizing/non-ionizing gas deposition,
  positive-hit count, gas-containment flags and `GAGGHit`;
- `Hits`: positive-energy-deposition gas hits with track, parent, position,
  deposited/kinetic energy, time, direction, vertex and process identifiers.

Only events with at least one positive gas hit are stored in `Events`, and only
positive-energy-deposition hits are stored in `Hits`. The normalization
denominator is `RunInfo.nBeamOnRequested`, not the number of rows in
`Events`.

`GAGGHit` is a bit-like integer:

| Value | Meaning |
| ---: | --- |
| `0` | no GAGG tag |
| `1` | GAGG cylindrical wall |
| `2` | GAGG -Z cap |
| `3` | both wall and cap |

For launch metadata, prefer `run_configuration.json`: several `RunInfo`
configuration fields are available for manual macro use and are not all
explicitly populated by the generated macro.

## Analyze a run

Run the analysis from the Python 3.9+ environment after all ROOT files have
finished:

```bash
source .venv/bin/activate

python3 analyze_run.py \
  --directory run_outputs/example \
  --spectra-dir spectra \
  --energy-min 6 \
  --energy-max 30 \
  --bin-width 1 \
  --s-smearing 0.15
```

Defaults are 0–100 keV, 1 keV bins and 15% relative Gaussian energy smearing.
ROOT files are discovered recursively beneath `--directory`.

The analysis produces five selections:

| Case in CSV | Selection |
| --- | --- |
| `standard` | deposited energy in the sensitive gas |
| `timepix_cut` | at least one hit inside the 14.08 mm × 14.08 mm Timepix footprint |
| `timepix_diffusion` | Timepix selection with Gaussian energy smearing |
| `timepix_diffusion_gagg` | smeared Timepix selection with `GAGGHit == 0`, the veto sample |
| `timepix_diffusion_gaggIN` | smeared Timepix selection with `GAGGHit != 0`, the tagged sample |

The historical word `diffusion` in these case names currently means Gaussian
energy-resolution smearing:

```text
E_measured ~ Normal(E_deposited, s_smearing * E_deposited)
```

It is not a spatial charge-diffusion model. The analysis currently does not set
its own random seed, so smeared histograms can vary slightly between repeated
analysis runs even when the Geant4 outputs are identical.

For each source and energy bin, the physical differential rate is calculated
from the number of selected events, the requested primary count, the integrated
source spectrum and the generation patch. The bin probability uses
`(k + 1) / (N + 2)`, with the corresponding finite-sample variance. Source
rates are summed and uncertainties are added in quadrature.

Results are written below `<run>/plots/`:

```text
flux_gas_volume_*.png
flux_timepix_*.png
flux_timepix_smeared_*.png
flux_timepix_smeared_gagg_*.png
flux_timepix_smeared_gaggIN_*.png
integrated_rates_*.csv
```

The plot ordinate is a detector rate density in s^-1 keV^-1. The CSV contains
the integrated rate and its uncertainty in s^-1 for every case and source,
plus a `Total` row.

### LaTeX summary tables

Convert the integrated-rate CSV into Hz and normalized
Hz cm^-2 keV^-1 tables:

```bash
python3 csv_to_latex_summary.py \
  run_outputs/example/plots/integrated_rates_6.0_30.0_1.00_0.15.csv \
  --out run_outputs/example/plots/rates_hz.tex \
  --out-density run_outputs/example/plots/rates_density.tex
```

The normalized table divides by the detector area, 14.08 mm × 14.08 mm =
1.982464 cm² by default, and by the selected energy interval. Override the
area with `--area-cm2`. The generated LaTeX uses `booktabs` and `graphicx`.

## Compare the three tungsten geometries

The repository provides a host launcher that builds once in Docker, freezes
one copy of the input spectra and starts all three configurations in parallel:

```bash
./scripts/launch_geometry_comparison.sh
```

The two optional positional arguments are the output root and the worker count
per configuration:

```bash
./scripts/launch_geometry_comparison.sh my_comparison 4
```

The default output directory contains a timestamp. The default worker count is
one third of the available logical CPUs for each configuration. The launcher
refuses to reuse a non-empty output directory.

The following environment variables can override its defaults:

| Variable | Default | Meaning |
| --- | --- | --- |
| `GPD3D_BASE_SEED` | `20260728` | common deterministic base seed |
| `GPD3D_BUILD_JOBS` | host CPU count | parallel compiler jobs |
| `GPD3D_DOCKER_IMAGE` | `u1804-root62004-g4_1051` | Docker image tag |
| `GPD3D_DOCKER_BUILD_DIR` | `build_docker_geometry` | shared build directory name |

For example:

```bash
GPD3D_BASE_SEED=123456 \
GPD3D_BUILD_JOBS=8 \
./scripts/launch_geometry_comparison.sh comparison_01 4
```

The launcher:

1. verifies host `screen` and Docker access;
2. compiles one shared executable inside the Docker image;
3. copies `spectra/` to `<output>/input_spectra/`;
4. writes `manifest.sha256`;
5. starts one detached Docker container per geometry as the host UID/GID;
6. starts one detached Screen log monitor per geometry;
7. writes the session/container/path mapping to `screen_sessions.txt`.

The output layout is:

```text
my_comparison/
├── input_spectra/
├── manifest.sha256
├── screen_sessions.txt
├── plug_and_tube/
├── plug_only/
└── no_plug_no_tube/
```

Each configuration directory has the normal run layout, plus `screen.log` and
an `exit_status` file when it finishes.

Useful monitoring commands:

```bash
screen -ls
screen -r <session-name>
tail -f my_comparison/plug_and_tube/screen.log
docker ps --filter name=gpd3d-sim
```

Detach from Screen without stopping the run with `Ctrl-A`, then `D`.
If group membership was added after the current login, either open a new login
session or run the launcher through a refreshed group shell:

```bash
sg docker -c './scripts/launch_geometry_comparison.sh my_comparison 4'
```

To compare like with like, analyze each configuration with the frozen spectra,
not a subsequently edited working copy:

```bash
python3 analyze_run.py \
  --directory my_comparison/plug_and_tube \
  --spectra-dir my_comparison/input_spectra \
  --energy-min 6 --energy-max 30 --bin-width 1 --s-smearing 0.15

python3 analyze_run.py \
  --directory my_comparison/plug_only \
  --spectra-dir my_comparison/input_spectra \
  --energy-min 6 --energy-max 30 --bin-width 1 --s-smearing 0.15

python3 analyze_run.py \
  --directory my_comparison/no_plug_no_tube \
  --spectra-dir my_comparison/input_spectra \
  --energy-min 6 --energy-max 30 --bin-width 1 --s-smearing 0.15
```

## Current implementation notes

- Gas hits are scored only in `GAS_SENSITIVE_LV`, a box spanning the Timepix
  footprint through the drift depth, not in the full cylindrical
  `GAS_CHAMBER_LV`. Consequently, `standard` already represents the sensitive
  prism, and `timepix_cut` is normally equivalent or very close to it. The
  historical plot title “entire gas volume” should be interpreted with this
  implementation detail in mind.
- `GAGGHit` is currently formed in `GPD3D_SteppingAction` from energy deposited
  by the primary track (`trackID == 1`) in the wall or cap. Energy deposited
  there only by secondary tracks does not set the tag.
- `RunInfo.nBeamOnRequested` is the authoritative generated-event denominator
  for analysis. Do not normalize by the sparse `Events` tree length.
- The simulation base seed makes Geant4 comparisons reproducible, but the
  Gaussian analysis smearing currently remains stochastic.
- `launch_run.py` still prints paths for two overview PNG files, but its
  overview plotting function is presently a no-op; the authoritative launcher
  products are the CSV/JSON summaries and ROOT files.
- Files in `old_scripts/`, `debug/` and local convenience shell scripts are
  retained for reference and may assume old paths or data formats. Use
  `launch_run.py`, `analyze_run.py` and `scripts/launch_geometry_comparison.sh`
  for the documented workflow.

## Troubleshooting

### Docker permission denied

Check:

```bash
docker info
id
```

If the user appears in the `docker` group but `docker info` still fails, the
current shell has stale group membership. Log in again, use `newgrp docker`, or
use the `sg docker -c ...` command shown above.

### CMake fails in a host build directory

Use the provided Docker image. The comparison launcher deliberately builds in
`build_docker_geometry/` inside Docker, avoiding an incompatible or
root-owned host `build/` tree.

### Analysis cannot match a ROOT file to a spectrum

Check that:

- `--spectra-dir` contains both the relevant CSV files and `manifest.csv`;
- the spectrum stem matches the source part of the ROOT filename;
- source names contain no underscores;
- the manifest contains angular limits for every simulated source.

### Screen session ends or reports failure

Inspect the saved log and status:

```bash
tail -n 100 my_comparison/plug_only/screen.log
cat my_comparison/plug_only/exit_status
```

An `exit_status` of `0` means success. A missing status while the container is
still listed by `docker ps` means the simulation is still running.
