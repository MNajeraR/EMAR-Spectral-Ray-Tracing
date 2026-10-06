# EMAR Spectral Ray-Tracing Analysis
[![DOI](https://zenodo.org/badge/1363089025.svg)](https://doi.org/10.5281/zenodo.23198727)
## Overview

This repository contains a Python-based ray-tracing and spectral-analysis
framework developed for the Espectrógrafo Mexicano de Alta Resolución (EMAR)
optical system. The code interfaces with Zemax OpticStudio through PyZDDE and
provides tools for configurable field and normalized-pupil sampling,
multi-configuration ray tracing, wavelength and diffraction-order control
through the Zemax Multi-Configuration Editor (MCE), spectral-dispersion
analysis, and characterization of extended input sources. The Zemax model
contains 13 representative configurations associated with echelle diffraction
orders spanning from `m = 60` to `m = 144`, with each configuration defining
11 reference wavelengths through the MCE. The framework was developed
progressively, beginning with individual pupil and wavelength tests and
extending to the reconstruction and analysis of the complete integer
echelle-order sequence.

The current workflow includes:

- configurable field and normalized-pupil sampling using random or hexapolar
  distributions;
- ray tracing of individual and multiple Zemax configurations;
- reconstruction of the complete integer echelle-order sequence from `m = 60`
  to `m = 144`;
- dense wavelength sampling across individual and complete echelle orders;
- visualization of ray footprints and spectral centroid traces at selected
  optical surfaces;
- detector-style visualization of the complete spectral format;
- calculation of linear spectral dispersion and resolving power;
- tracing and characterization of finite slit and circular-fiber input sources;
- comparison of dense and reduced fiber-field sampling strategies; and
- analysis of the free separation between adjacent echelle-order fiber
  envelopes.

A circular `100 µm` fiber is used to establish a single-fiber reference across
all 85 integer diffraction orders, providing the basis for subsequent studies
of multiple fiber traces and polarimetric configurations. The original Zemax
optical model is intentionally not distributed with this repository; users
with authorized access to the EMAR model can place their local copy in the
`Zemax/` directory and use the analysis scripts without modifying the
repository structure.

## Contents

- [Repository Structure](#repository-structure)
  - [`scripts/`](#scripts)
  - [`experiments/`](#experiments)
  - [`utils/`](#utils)
  - [`Results/`](#results)
  - [`Zemax/`](#zemax)
- [Running the Analysis](#running-the-analysis)
- [Zemax Model Interaction](#zemax-model-interaction)
  - [Configuration selection](#configuration-selection)
  - [Reading the MCE](#reading-the-mce)
  - [Temporary wavelength modification](#temporary-wavelength-modification)
  - [Ray tracing](#ray-tracing)
- [Current Analysis Status](#current-analysis-status)

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
│   ├── EMAR_missing_orders.py
│   └── EMAR_full_orders.py
│
├── experiments/
│   ├── 01_pupil_sampling.py
│   ├── 02_single_config_original_waves.py
│   ├── 03_single_config_100_waves.py
│   ├── 04_two_configs_100_waves.py
│   ├── 05_intermediate_order_61.py
│   ├── 06_missing_orders_61_66.py
│   ├── 07_all_orders_11_wavelengths.py
│   ├── 08_single_order_dispersion.py
│   ├── 09_single_order_spectral_resolution.py
│   ├── 10_theoretical_slit.py
│   ├── 11_fiber_resolution.py
│   ├── 12_full_spectral_format.py
│   ├── 13_five_field_fiber_spectrum.py
│   ├── 14_compare_fiber_sampling.py
│   └── 15_order_separation.py
│
├── utils/
│   ├── __init__.py
│   ├── emar_utils.py
│   └── emar_plots.py
│
├── Results/
│   ├── echelle_orders_60_144.csv
│   ├── order_102_dispersion.csv
│   ├── order_102_rays.csv
│   ├── order_102_theoretical_slit_resolution.csv
│   ├── order_102_fiber_resolution.csv
│   ├── full_orders_fiber_resolution.csv
│   ├── full_orders_fiber_resolution_5fields.csv
│   ├── fiber_sampling_comparison.csv
│   ├── fiber_sampling_statistics_by_order.csv
│   └── order_separation_91fields.csv
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
  analyzes the `m * lambda` relation, and reconstructs the spectral
  definition of every integer echelle order from `m = 60` to `m = 144`.

- `EMAR_full_orders.py` performs dense ray tracing of the complete integer
  echelle-order sequence from `m = 60` to `m = 144`. The script reads the
  order definitions from `Results/echelle_orders_60_144.csv`, which identifies
  each order as original or reconstructed and provides its 11 reference
  wavelengths. For each order, a dense wavelength grid is generated between
  the first and last reference wavelengths and propagated to the final image
  plane. Reconstructed orders are traced directly through the Zemax optical
  model rather than obtained by interpolating detector coordinates. The
  resulting ray data can be visualized as complete ray footprints, spectral
  centroid traces, or a simplified detector representation.

### `experiments/`

Contains the progressive validation experiments used during development of the
EMAR ray-tracing and spectral-analysis framework. The experiments document the
transition from basic pupil-sampling and single-configuration tests to
diffraction-order reconstruction, spectral-dispersion and resolution analysis,
finite-source modeling, complete-format fiber tracing, and adjacent-order
separation analysis. Each experiment isolates a specific step of the
methodology before extending it to the complete echelle format.

1. `01_pupil_sampling.py` visualizes and validates the normalized-pupil
   sampling methods used by the ray-tracing routines. Random and deterministic
   pupil distributions can be inspected independently of the Zemax model,
   providing a basic validation of the coordinates subsequently passed to
   `zGetTrace()`.

2. `02_single_config_original_waves.py` traces the original reference
   wavelengths of a single Zemax configuration. This experiment verifies the
   basic PyZDDE connection, configuration selection, wavelength handling, and
   extraction of valid ray coordinates at the selected optical surface.

3. `03_single_config_100_waves.py` introduces dense wavelength sampling
   within a single original configuration. Intermediate wavelengths are
   generated across the spectral interval defined by the original reference
   wavelengths, allowing the evolution of the spectral trace to be examined
   continuously rather than only at the discrete wavelengths stored in the
   Zemax model.

4. `04_two_configs_100_waves.py` extends the dense wavelength procedure to
   two original Zemax configurations. The experiment verifies that the same
   tracing sequence can be applied while changing configurations and that
   temporary wavelength modifications can be performed without altering the
   original MCE definition.

5. `05_intermediate_order_61.py` reconstructs and ray traces the first
   intermediate diffraction order, `m = 61`, between the original `m = 60`
   and `m = 67` configurations. The experiment validates temporary
   modification of both the echelle diffraction order and wavelength through
   the MCE and demonstrates that an integer order not explicitly stored as an
   original Zemax configuration can be reconstructed.

6. `06_missing_orders_61_66.py` extends the reconstruction to all
   intermediate integer orders between `m = 60` and `m = 67`. Their detector
   traces continuously populate the spectral region between the two original
   configurations, validating the order-reconstruction procedure over a
   complete configuration interval.

7. `07_all_orders_11_wavelengths.py` applies the validated reconstruction
   procedure to the complete integer sequence from `m = 60` to `m = 144`.
   All 85 diffraction orders are traced using their 11 reference wavelengths
   at the final image plane. This provides the complete reconstructed echelle
   format and constitutes the final validation step before introducing dense
   wavelength sampling and quantitative spectral analysis.

8. `08_single_order_dispersion.py` introduces quantitative spectral-dispersion
   analysis using the representative order `m = 102` and a point source at
   the spectrograph input. The order is sampled with 100 wavelengths while a
   common normalized-pupil distribution is traced at each wavelength. Detector
   ray coordinates are used to calculate the wavelength-dependent spectral
   centroid and its trajectory across the image plane. The centroid displacement 
   is then used to derive the linear dispersion as a function of wavelength.

9. `09_single_order_spectral_resolution.py` extends the point-source analysis
   of `m = 102` from dispersion to spectral resolution. At each wavelength,
   the monochromatic ray footprint produced by the point source is projected
   onto the local direction of spectral dispersion. The projected footprint
   width is used to estimate the wavelength interval associated with the
   monochromatic image and the corresponding resolving power. This experiment
   establishes the resolution-analysis procedure before replacing the point
   source with a finite entrance source.

10. `10_theoretical_slit.py` replaces the point source with a finite
    rectangular slit at the spectrograph input. A theoretical `50 × 15 µm`
    slit is mapped onto the corresponding Zemax field coordinates and sampled
    across its extent. Rays from each slit position are traced through the
    pupil and onto the detector. The resulting footprints are projected along
    the local spectral direction to determine the projected slit width,
    spectral purity, and resolving power across order `m = 102`.

11. `11_fiber_resolution.py` replaces the rectangular slit with the circular
    fiber geometry relevant to the instrument. A `100 µm` input fiber is
    represented by a deterministic five-ring hexapolar distribution containing
    91 source positions. Each fiber position is traced using an independent
    five-ring hexapolar pupil sample of 91 rays, producing 8,281
    rays per wavelength. The experiment first validates the reconstructed
    fiber geometry at the input surface and then traces 100 wavelengths across
    `m = 102`. The detector distributions are used to measure the
    wavelength-dependent X and Y dimensions of the projected fiber image and
    to calculate its corresponding spectral resolving power.

12. `12_full_spectral_format.py` extends the `100 µm` fiber analysis from
    `m = 102` to all integer orders from `m = 60` to `m = 144`. Each order is
    sampled at 100 wavelengths, retaining the 91 fiber positions and 91 pupil
    rays established in the preceding experiment. The resulting detector
    coordinates are used to calculate the fiber centroid, X-Y envelope, linear
    spectral dispersion, projected spectral width, and resolving power across
    all 85 diffraction orders. These results also provide the fiber envelopes
    required for the subsequent analysis of spacing between neighboring orders.

13. `13_five_field_fiber_spectrum.py` evaluates a reduced source-sampling
    strategy by replacing the 91 fiber positions with five representative
    fields: the center and the positive and negative extrema along X and Y.
    The 91-ray hexapolar pupil and 100 wavelengths per order are retained.
    The complete `m = 60–144` sequence is then traced to provide a
    computationally faster alternative to the 91-field calculation.

14. `14_compare_fiber_sampling.py` compares the five-field and 91-field
    calculations without additional Zemax ray tracing. Their results are
    evaluated wavelength by wavelength across all 85 orders to quantify
    differences in projected fiber width and resolving power. Per-order
    statistics are also generated to determine how well the reduced sampling
    preserves the spectral behavior of the 91-field calculation. The comparison
    establishes which sampling strategy is appropriate for subsequent
    fiber-envelope and order-separation analyses.

15. `15_order_separation.py` uses the 91-field fiber envelopes to quantify
    the free space between adjacent echelle orders. For every neighboring pair
    from `m = 60–61` through `m = 143–144`, the two envelopes are evaluated
    over their common detector-X interval. Their Y boundaries are interpolated
    at common X coordinates, and the edge-to-edge separation is calculated
    between the neighboring envelope boundaries. Minimum, maximum, mean, and
    median separations are determined for each order pair together with the
    detector-X locations of the extrema. Representative separation profiles
    are also evaluated across the detector, providing the reference for
    subsequent multi-fiber and polarimetric configurations.


### `utils/`

Contains reusable functions shared by the analysis and experimental
scripts.

- `emar_utils.py` contains project paths, pupil-sampling routines,
  PyZDDE ray-tracing functions, multi-configuration tracing tools, and
  reusable functions for tracing original and reconstructed echelle
  diffraction orders while preserving and restoring the original MCE
  state.

- `emar_plots.py` contains reusable visualization routines for pupil
  sampling, multi-configuration footprints, and the complete echelle
  spectral format. Echelle-order results can be displayed as full ray
  footprints, centroid spectral traces, or detector-style
  representations. Scientific echelle views use a continuous
  diffraction-order color scale, with lower orders represented toward
  red and higher orders toward blue.

### `Results/`

Contains numerical products and reproducible analysis outputs generated by the
scripts and validation experiments. The master spectral definition,
`echelle_orders_60_144.csv`, contains the complete integer echelle-order
sequence from `m = 60` to `m = 144`, identifying each order as `ORIGINAL` or
`NEW` and providing the 11 reference wavelengths used for subsequent spectral
sampling.

Single-order analysis products include:

- `order_102_dispersion.csv`, containing the wavelength-dependent detector
  centroid and spectral-dispersion measurements for `m = 102`;
- `order_102_rays.csv`, containing the corresponding individual detector ray
  coordinates;
- `order_102_theoretical_slit_resolution.csv`, containing the
  spectral-resolution analysis for the finite theoretical slit; and
- `order_102_fiber_resolution.csv`, containing the projected `100 µm` fiber
  dimensions and resolving-power analysis for `m = 102`.

Complete-format fiber products include:

- `full_orders_fiber_resolution.csv`, the 91-field fiber analysis for
  all 85 integer diffraction orders;
- `full_orders_fiber_resolution_5fields.csv`, the corresponding reduced
  five-field calculation;
- `fiber_sampling_comparison.csv`, containing the wavelength-by-wavelength
  comparison between both fiber-sampling strategies;
- `fiber_sampling_statistics_by_order.csv`, containing the comparison
  statistics summarized by diffraction order; and
- `order_separation_91fields.csv`, containing the free-separation
  statistics for adjacent echelle-order fiber envelopes.

The analysis scripts also generate figures for local inspection of the
spectral format, dispersion, projected fiber dimensions, resolving power,
sampling comparisons, and adjacent-order separation. These generated image
files are excluded from version control through `.gitignore` and are not
distributed with the repository.


### `Zemax/`

Provides the expected local location for the EMAR Zemax optical model.
The `.gitkeep` file preserves this directory in Git while its optical
model contents remain excluded through `.gitignore`.

## Running the Analysis

Before running the EMAR analysis scripts, Zemax OpticStudio must be open and
the PyZDDE connection must be properly configured. Installation instructions and connection examples are available in the
[official PyZDDE repository](https://github.com/xzos/PyZDDE) and can be used
to verify the setup before executing the EMAR workflow. A local
copy of the EMAR Zemax model must also be available in the `Zemax/` directory.

Once the setup has been verified, `experiments/02_single_config_original_waves.py`
provides a minimal starting point for the EMAR workflow by selecting a single
configuration, reading its reference wavelengths, and tracing rays to the
selected optical surface. The subsequent experiments progressively extend this
procedure to dense wavelength sampling, diffraction-order reconstruction,
spectral analysis, and finite-source tracing.

## Zemax Model Interaction

The analysis scripts interact with the Zemax multi-configuration model through
PyZDDE. Configurations and MCE parameters are read and temporarily modified as
required by each analysis, while the original model state is restored after
ray tracing.

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
wavelengths associated with the selected configuration. For dense spectral 
sampling, the first and last reference wavelengths define the spectral interval:

    WAVE 1  -> minimum wavelength
    WAVE 11 -> maximum wavelength

For dense spectral sampling, 100 equally spaced wavelengths are generated
between these two limits.

### Temporary wavelength modification

Each sampled wavelength is temporarily assigned to `WAVE 1` using
`zSetMulticon()`. The model is then refreshed with `zGetUpdate()` before
the ray trace is performed. The sequence for a sampled wavelength is therefore:

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

For each sampled wavelength, rays are traced to the selected optical surface
using `zGetTrace()`. The sampling strategy is configurable independently for
the source field and the normalized pupil, with random and hexapolar
distributions available depending on the analysis. For point-source analyses,
a single field position is combined with random sampling of the normalized
pupil. Finite-slit analyses use random sampling across both the source field
and the pupil. For the circular-fiber analysis, both are sampled using
hexapolar distributions.

In the full-format fiber calculation, the `100 µm` circular fiber is sampled
with 91 field positions arranged in a five-ring hexapolar distribution, and
each field position is traced through 91 normalized-pupil positions generated
with an independent five-ring hexapolar distribution. Each diffraction order
is sampled at 100 wavelengths. A reduced five-field representation of the
fiber is also implemented and benchmarked against the 91-field calculation.
In this case, only the source-field sampling is reduced to the fiber center and
the four extrema along X and Y, while the 91-point hexapolar pupil sampling is
retained.

For the initial multi-configuration spectral reconstruction, the tracing
sequence can be represented conceptually as:

    Configuration 1
        ├──> wavelength 1   -> pupil sample -> ray tracing
        ├──> wavelength 2   -> pupil sample -> ray tracing
        ├── ...
        └──> wavelength 100 -> pupil sample -> ray tracing

    Configuration 2
        ├──> wavelength 1   -> pupil sample -> ray tracing
        ├──> wavelength 2   -> pupil sample -> ray tracing
        ├── ...
        └──> wavelength 100 -> pupil sample -> ray tracing

        ...

    Configuration 13
        ├──> wavelength 1   -> pupil sample -> ray tracing
        ├──> wavelength 2   -> pupil sample -> ray tracing
        ├── ...
        └──> wavelength 100 -> pupil sample -> ray tracing

For the extended-source fiber analysis, the tracing sequence includes the
additional fiber-position sampling:

    Diffraction order
        ├──> wavelength 1
        │       ├──> fiber position 1 -> pupil sample -> ray tracing
        │       ├──> fiber position 2 -> pupil sample -> ray tracing
        │       ├── ...
        │       └──> fiber position N -> pupil sample -> ray tracing
        │
        ├──> wavelength 2
        │       └──> fiber positions -> pupil samples -> ray tracing
        │
        ├── ...
        │
        └──> wavelength 100
                └──> fiber positions -> pupil samples -> ray tracing

where `N = 91` for the full hexapolar fiber sampling and `N = 5` for the
reduced fiber-field sampling experiment.

After the sampled wavelengths of a configuration have been traced, the
original `WAVE 1` value is restored in the MCE and the model is updated again
with `zGetUpdate()` before proceeding to the next configuration. This
restoration prevents temporary wavelength modifications from affecting
subsequent configurations or analyses.

## Current Analysis Status

The current workflow provides a complete single-fiber analysis of the EMAR
echelle format from `m = 60` to `m = 144`. A circular `100 µm` input fiber is
traced across all 85 integer diffraction orders using the 91-field sampling,
providing wavelength-dependent projected fiber dimensions, spectral dispersion,
resolving power, and adjacent-order separation. A reduced five-field sampling
was benchmarked against this solution, reproducing the wavelength-dependent
behavior of the projected fiber and resolving power at substantially lower
computational cost. The 91-field results are retained for fiber-envelope and
order-separation analyses.

Across the complete spectral format, adjacent single-fiber envelopes remain
separated, with a minimum edge-to-edge distance of approximately `141 µm`
between orders `m = 60` and `m = 61`. These results provide the reference for
the next development stage: tracing multiple fiber images within each echelle
order and evaluating their compatibility with a dual-beam polarimetric
configuration.
