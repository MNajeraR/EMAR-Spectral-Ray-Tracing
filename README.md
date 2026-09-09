# EMAR Spectral Ray-Tracing Analysis

## Overview

This repository contains a Python-based ray-tracing and spectral-analysis
framework developed for the Espectrógrafo Mexicano de Alta Resolución (EMAR) 
optical system. The code interfaces with Zemax OpticStudio through PyZDDE and provides
tools for sampling the normalized pupil, tracing rays through selected
surfaces of the optical system, modifying wavelength values in the Zemax
Multi-Configuration Editor (MCE), and analyzing the resulting spectral
format.

The current Zemax model contains 13 representative configurations
associated with echelle diffraction orders spanning from m = 60 to
m = 144. Each configuration defines 11 reference wavelengths through
the Multi-Configuration Editor. The analysis implemented in this repository was 
developed progressively, starting from individual pupil and wavelength tests and 
extending to the simultaneous analysis of all original configurations. 
The main workflow currently includes:

- normalized-pupil sampling using random or hexapolar distributions;
- ray tracing of individual and multiple Zemax configurations;
- dense wavelength sampling within the spectral interval of each
  configuration;
- visualization of ray footprints and spectral centroid traces at
  selected optical surfaces;
- detector-style visualization of the spectral format;
- analysis of the relation between echelle diffraction order and
  wavelength; and
- reconstruction of the wavelength sampling for the complete sequence
  of integer echelle orders from m = 60 to m = 144.

The original Zemax optical model is intentionally not distributed with
this repository. Users with authorized access to the EMAR model can
place their local copy in the `Zemax/` directory and use the analysis
scripts without modifying the repository structure.

## Repository Structure

The repository is organized to separate the main analysis scripts,
reusable functions, validation experiments, generated results, and the
local Zemax optical model.

```text
EMAR/
├── README.md
├── requirements.txt
├── .gitignore
│
├── scripts/
│   ├── EMAR_All_confs_surf26.py
│   ├── EMAR_All_confs_surf54.py
│   └── EMAR_missing_orders.py
│
├── experiments/
│   ├── 01_pupil_sampling.py
│   ├── 02_single_config_original_waves.py
│   ├── 03_single_config_100_waves.py
│   └── 04_two_configs_100_waves.py
│
├── utils/
│   ├── __init__.py
│   ├── emar_utils.py
│   └── emar_plots.py
│
├── Results/
│   └── echelle_orders_60_144.csv
│
└── Zemax/
    └── .gitkeep
```
### `scripts/`

Contains the main analysis scripts used to study the EMAR optical
system at selected stages of the spectrograph.

- `EMAR_All_confs_surf26.py` traces the 13 original Zemax
  configurations to surface 26 using dense wavelength sampling.
  Surface 26 corresponds to the second focal plane of the system. Before
  reaching this surface, the beam is collimated by the first off-axis
  parabola (OAP1), dispersed by the echelle diffraction grating, and
  refocused by the second off-axis parabola (OAP2). This surface is
  therefore used to inspect the spectral behavior produced by the
  echelle dispersion before the subsequent cross-dispersion stage.

- `EMAR_All_confs_surf54.py` propagates the same configurations to
  surface 54, which corresponds to the final image plane. After the
  second focal plane, the beam is recollimated by a third off-axis
  parabola (OAP3), cross-dispersed by two prisms, and finally focused
  by a camera composed of one triplet and three doublets. At this
  surface, the combined echelle and cross-dispersion produces the
  two-dimensional spectral format analyzed by the script.

- `EMAR_missing_orders.py` reads the original MCE configurations,
  analyzes the `m * lambda` relation, and reconstructs the wavelength
  sampling for every integer echelle order from `m = 60` to `m = 144`.


### `experiments/`

Contains the progressive validation experiments used during development
of the ray-tracing methodology. These scripts document the transition 
from simple pupil-sampling tests to the complete multi-configuration 
analysis:

1. Visualization and validation of pupil-sampling methods.
2. Tracing the original wavelengths of a single Zemax configuration.
3. Dense spectral sampling of a single configuration.
4. Simultaneous dense spectral sampling of two configurations.

They are retained to document and reproduce the validation process
leading to the main analysis scripts.


### `utils/`

Contains reusable functions shared by the analysis and experimental
scripts.

- `emar_utils.py` contains project paths, pupil-sampling routines,
  PyZDDE ray-tracing functions, multi-configuration tracing tools, and
  spectral-centroid calculations.

- `emar_plots.py` contains the visualization routines for pupil
  sampling, footprints, spectral traces, and detector-style
  representations.

Keeping these functions separate from the analysis scripts avoids code
duplication and provides a common implementation for all ray-tracing
experiments.


### `Results/`

Contains compact numerical products generated by the analysis.

The current `echelle_orders_60_144.csv` file contains the reconstructed
wavelength table for the complete sequence of echelle diffraction
orders from `m = 60` to `m = 144`.

Generated PNG and PDF figures are excluded from version control because
they can be reproduced directly from the analysis scripts.


### `Zemax/`

Provides the expected local location for the EMAR Zemax optical model.
The `.gitkeep` file preserves this directory in Git while its optical
model contents remain excluded through `.gitignore`.

## Zemax Model Interaction

The analysis scripts interact with the Zemax multi-configuration model
through PyZDDE. The workflow is organized as a nested process in which
each configuration is selected, spectrally sampled, and ray traced
before proceeding to the next one:

1. Select a Zemax configuration.
2. Read its Multi-Configuration Editor (MCE) data.
3. Generate the wavelength sampling for that configuration.
4. Sequentially modify and trace each sampled wavelength.
5. Restore the original MCE values.
6. Continue with the next configuration.

### Configuration selection

The configurations to be analyzed are defined in Python, for example:

    configs = range(1, 14)

The script iterates over this sequence and activates each configuration
with `zSetConfig()`. After changing the active configuration,
`zGetUpdate()` updates the optical system according to the parameters
associated with that configuration.

### Reading the MCE

Once a configuration is selected, `zGetMulticon()` reads the required
values directly from its MCE entries. Depending on the analysis, these
values can include the echelle diffraction order and the reference
wavelengths associated with the selected configuration.

For dense spectral sampling, the first and last reference wavelengths
define the spectral interval:

    WAVE 1  -> minimum wavelength
    WAVE 11 -> maximum wavelength

A set of 100 equally spaced wavelengths is generated between these two
limits for the active configuration.

### Temporary wavelength modification

Each sampled wavelength is temporarily assigned to `WAVE 1` using
`zSetMulticon()`. The model is then refreshed with `zGetUpdate()` before
the ray trace is performed.

The sequence for a sampled wavelength is therefore:

    zSetConfig()
         |
         v
    zGetMulticon()
         |
         v
    zSetMulticon()
         |
         v
    zGetUpdate()
         |
         v
    zGetTrace()

### Ray tracing

For each sampled wavelength, a set of points from the normalized pupil
is traced to the selected optical surface using `zGetTrace()`. The pupil
sampling is configurable and can be generated using either a random or
a hexapolar distribution. The number of sampled pupil points is also
defined by the user.

In the current analysis scripts, 100 pupil points are used for each of
the 100 sampled wavelengths. The valid X-Y ray coordinates are stored
in Python for subsequent analysis and visualization. Conceptually, 
the complete multi-configuration analysis follows:

    Configuration 1
        ├── wavelength 1   -> pupil sample -> ray tracing
        ├── wavelength 2   -> pupil sample -> ray tracing
        ├── ...
        └── wavelength 100 -> pupil sample -> ray tracing

    Configuration 2
        ├── wavelength 1   -> pupil sample -> ray tracing
        ├── wavelength 2   -> pupil sample -> ray tracing
        ├── ...
        └── wavelength 100 -> pupil sample -> ray tracing

        ...

    Configuration 13
        ├── wavelength 1   -> pupil sample -> ray tracing
        ├── wavelength 2   -> pupil sample -> ray tracing
        ├── ...
        └── wavelength 100 -> pupil sample -> ray tracing

After the sampled wavelengths of a configuration have been traced, the
original `WAVE 1` value is restored in the MCE and the model is updated
again with `zGetUpdate()` before proceeding to the next configuration.
This restoration prevents temporary wavelength modifications from
affecting subsequent configurations or analyses.