# GEANT4 for GPD3D MedEnergy

This code has been developed to simulate the backgorund in LEO orbit for GPD3D medium energy 6-30keV
If background mode is true, particles from energy and type defined in the spectra/manifest.csv file are generated inwards a sphere inside the detector.
One particle shor is an Event, general info about events are saved along with detailed information on every hit in the Event

NOTE: this code is derived Daweoon code, it contasins also the old digitazion part used for polarimetry mode: /gpd3d/bg/enable false some day in the future should becleand up

## Installation

You can run this code inside in a Docker container
- Install Docker https://docs.docker.com/get-started/get-docker/
- docker build -t u1804-root62004-g4_1051 .
- xhost +local:

  docker run --rm -it \
    -e DISPLAY \
    -e LIBGL_ALWAYS_SOFTWARE=1 \
    -e MESA_LOADER_DRIVER_OVERRIDE=llvmpipe \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v $HOME/.Xauthority:/root/.Xauthority:ro \
    -v "$(pwd):/work" -w /work \
    u1804-root62004-g4_1051

Build the Geant4
- mkdir build && cd build
- cmake ..
- make
- run with ./gpd3d macros/vis.mac

## Geometry
- Timepix 37cm3 gas volume with no shielding
- added passive 15 cm side cube of 1.5mm thick Al with only the circular window open
- Added anticoincidence 1cm GAGG cylinder surroungin the detector 
TBA:
- Make GAGG senitive ad store info

## Added commands

- /gpd3d/bg/enable true [background mode or polarimetry mode]
- /gpd3d/bg/writeFITS false [wrtie FITS]
- /gpd3d/gen/useSphere [true|false] [If true, primaries are generated on a sphere around the detector with cosine-law inward directions (4π isotropic environment). If false, you use the Gaussian beam mode.]
- /gpd3d/gen/sphereRadius <value> mm [Radius of the generation sphere (must fit inside the world volume).]
- /gpd3d/gen/sphereCenter x y z mm [Center of the generation sphere.]
- /gpd3d/gen/beamSigma <value> mm [Gaussian σ for the beam spot in x and y.]
- /gpd3d/gen/beamCenter x y z mm [Beam center; x and y are used, z is ignored (use beamZ).]
- /gpd3d/gen/beamZ <value> mm [Fixed z position of the beam source.]
- /gpd3d/bg/rootFile test.root [name of output file root]
- /gpd3d/gen/useSpectrum [true|false] [If true, primaries are generated with energies distributed as input spectrum]
- /gpd3d/gen/spectrumFile [/ABS/PATH/TO/spectrum.csv] [Path to the input spectrum for energy distribution]

## Analysis

- rate_from_spectra_onlyRate.py run mulitple instances of the geant4 code using macros/template_bg_sphere.mac and the spectra/manifest.csv file to get the rate of events in the detector for each particle type and energy range. The output is a csv file with the rates.

- rate_from_spectra_root.py run mulitple instances of the geant4 code using macros/template_bg_sphere.mac and the spectra/manifest.csv file but save root files with detail inforation, if used well can generate an entire "run"

- rate_from_spectra_root_fromcsv.py run mulitple instances of the geant4 code using macros/template_bg_sphere_spectra.mac and the spectra/manifest.csv file but save root files with detail inforation, if used well can generate an entire "run"

- analyze_run.py is analyzing one "run" by gahtering togheder all the info, computingexpected rate, hit probability and for now also histograms of deposited energy.

- analyze_run_spectra.py is analyzing one "run" by gahtering togheder all the info, computing expected rate, hit probability and for now also histograms of deposited energy. Suggested usage:

```
python3 analyze_run_spectra.py --energy-min 6 --energy-max 30 --bin-width 1 --s-smearing 0.15 --directory run_outputs_continous_spectra_x10/ --plot-dir run_outputs_continous_spectra_x10/plots/
```

-run_analyzer_spectra.sh is a script looping over analyze_run_spectra.py arguments and running it
