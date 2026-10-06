# EMAR Spectral Ray-Tracing Analysis
## Overview

This repository contains a Python-based ray-tracing and spectral-analysis
framework developed for the Espectrógrafo Mexicano de Alta Resolución (EMAR)
optical system. The code interfaces with Zemax OpticStudio through PyZDDE and
provides tools for normalized-pupil sampling, multi-configuration ray tracing,
wavelength and diffraction-order control through the Zemax Multi-Configuration
Editor (MCE), spectral-dispersion analysis, and geometrical characterization of
extended input sources. The Zemax model contains 13 representative configurations associated with
echelle diffraction orders spanning from `m = 60` to `m = 144`. Each
configuration defines 11 reference wavelengths through the MCE. The analysis
implemented in this repository was developed progressively, beginning with
individual pupil and wavelength tests and extending to the reconstruction and
analysis of the complete integer echelle-order sequence.

The current workflow includes:

- normalized-pupil sampling using random or hexapolar distributions;
- ray tracing of individual and multiple Zemax configurations;
- reconstruction of intermediate integer echelle orders from `m = 60` to
  `m = 144`;
- dense wavelength sampling across individual and complete echelle orders;
- visualization of ray footprints and spectral centroid traces at selected
  optical surfaces;
- detector-style visualization of the complete spectral format;
- calculation of linear spectral dispersion;
- spectral-resolution analysis;
- tracing of finite slit and circular-fiber input sources;
- wavelength-dependent characterization of the projected fiber image;
- calculation of spectral resolving power across the complete echelle format;
- comparison of dense and reduced fiber-field sampling strategies; and
- analysis of the free separation between adjacent echelle-order
  fiber envelopes.

The current extended-source analysis uses a circular `100 µm` fiber and
provides a single-fiber reference across all 85 integer diffraction
orders. This framework is intended to support subsequent studies of multiple
fiber traces and polarimetric configurations. The original Zemax optical model is intentionally not distributed with this
repository. Users with authorized access to the EMAR model can place their
local copy in the `Zemax/` directory and use the analysis scripts without
modifying the repository structure.

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
  analyzes the `m * lambda` relation, and reconstructs the wavelength
  sampling for every integer echelle order from `m = 60` to `m = 144`.

- `EMAR_full_orders.py` performs dense ray tracing of the complete
  integer echelle-order sequence from `m = 60` to `m = 144`. The
  script reads the previously reconstructed order definitions directly
  from `Results/echelle_orders_60_144.csv`, which identifies each order
  as original or reconstructed and provides its 11 reference
  wavelengths. For each order, a dense wavelength grid is generated
  between its first and last reference wavelengths and propagated to
  the final image plane. The resulting ray data can be visualized as
  complete ray footprints, spectral centroid traces, or a simplified
  detector representation.

### `experiments/`

Contains the progressive validation experiments used during development
of the ray-tracing and order reconstruction methodology. These
scripts document the transition from simple pupil-sampling tests to the
validated reconstruction of the complete spectral format:

1. `01_pupil_sampling.py` visualizes and validates the normalized-pupil
   sampling methods.

2. `02_single_config_original_waves.py` traces the original reference
   wavelengths of a single Zemax configuration.

3. `03_single_config_100_waves.py` introduces dense wavelength sampling
   within a single original configuration.

4. `04_two_configs_100_waves.py` extends the dense sampling procedure
   to two Zemax configurations.

5. `05_intermediate_order_61.py` reconstructs and ray traces the first
   intermediate diffraction order, `m = 61`, between the original
   `m = 60` and `m = 67` configurations. This experiment validates
   temporary modification of both the echelle diffraction order and
   wavelength through the MCE.

6. `06_missing_orders_61_66.py` extends the reconstruction to all
   intermediate orders between `m = 60` and `m = 67`, demonstrating
   that the reconstructed spectral traces continuously fill the gap
   between the two original Zemax orders.

7. `07_all_orders_11_wavelengths.py` applies the validated reconstruction
   procedure to the complete integer order sequence from `m = 60` to
   `m = 144`. All 85 orders are traced using their 11 reference
   wavelengths at the final image plane, providing the final validation
   step before dense wavelength sampling.

8. `08_single_order_dispersion.py` introduces the spectral-dispersion
   analysis for a representative echelle order. The script densely samples
   order `m = 102`, traces a common normalized-pupil distribution, calculates
   the wavelength-dependent spectral centroid at the detector, and derives
   the corresponding linear and reciprocal spectral dispersion.

9. `09_single_order_spectral_resolution.py` extends the single-order
   analysis to geometrical spectral resolution. Monochromatic ray footprints
   are projected along the local spectral-dispersion direction to characterize
   their extent and evaluate the wavelength-dependent resolving behavior of
   the optical system.

10. `10_theoretical_slit.py` introduces a finite theoretical slit at the
    spectrograph input. A `50 × 15 µm` rectangular slit is mapped to Zemax
    field coordinates and sampled across its physical extent. The resulting
    detector footprints are used to evaluate the projected slit width,
    spectral purity, and geometrical resolving power.

11. `11_fiber_resolution.py` replaces the rectangular slit with a circular
    `100 µm` fiber. The fiber surface is sampled using a deterministic
    hexapolar distribution, while each fiber position is independently
    sampled across the normalized pupil. The experiment characterizes the
    wavelength-dependent projected fiber dimensions and spectral resolving
    power for the representative order `m = 102`.

12. `12_full_spectral_format.py` extends the circular-fiber analysis to the
    complete integer echelle sequence from `m = 60` to `m = 144`. The robust
    calculation uses 91 hexapolar fiber positions and 91 pupil rays at each
    wavelength, with 100 wavelengths sampled per order. The resulting data
    provide the projected fiber envelope, spectral dispersion, and resolving
    power throughout the complete spectral format.

13. `13_five_field_fiber_spectrum.py` implements a reduced fiber-field
    sampling strategy using five representative source positions: the fiber
    center and the four extrema along the X and Y directions. The same
    91-ray hexapolar pupil sampling and complete `m = 60–144` wavelength
    coverage are retained, providing a computationally faster approximation
    to the dense fiber-envelope calculation.

14. `14_compare_fiber_sampling.py` directly compares the five-field
    approximation with the robust 91-field calculation. The experiment
    quantifies wavelength-dependent and per-order differences in projected
    fiber width and spectral resolving power, providing a benchmark for the
    reduced sampling strategy.

15. `15_order_separation.py` uses the robust 91-field fiber envelopes to
    measure the free geometrical separation between adjacent echelle orders.
    Neighboring order boundaries are interpolated at common detector-X
    coordinates, allowing the minimum, maximum, mean, and median free
    separations to be calculated for every adjacent order pair. Representative
    separation profiles are also evaluated across the detector to identify
    the spatial location of limiting cases.


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
scripts and validation experiments. The master spectral definition is stored in
`echelle_orders_60_144.csv`. This table defines the complete integer
echelle-order sequence from `m = 60` to `m = 144`, identifying each order as
`ORIGINAL` or `NEW` and providing the 11 reference wavelengths used for
subsequent dense spectral sampling.

Single-order analysis products include:

- `order_102_dispersion.csv`, containing the wavelength-dependent detector
  centroid and spectral-dispersion measurements for `m = 102`;
- `order_102_rays.csv`, containing the corresponding individual detector ray
  coordinates;
- `order_102_theoretical_slit_resolution.csv`, containing the geometrical
  spectral-resolution analysis for the finite theoretical slit; and
- `order_102_fiber_resolution.csv`, containing the projected `100 µm` fiber
  dimensions and resolving-power analysis for `m = 102`.

Complete-format fiber products include:

- `full_orders_fiber_resolution.csv`, the robust 91-field fiber analysis for
  all 85 integer diffraction orders;
- `full_orders_fiber_resolution_5fields.csv`, the corresponding reduced
  five-field calculation;
- `fiber_sampling_comparison.csv`, containing the wavelength-by-wavelength
  comparison between both fiber-sampling strategies;
- `fiber_sampling_statistics_by_order.csv`, containing the comparison
  statistics summarized by diffraction order; and
- `order_separation_91fields.csv`, containing the geometrical free-separation
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
wavelengths associated with the selected configuration. For dense spectral 
sampling, the first and last reference wavelengths define the spectral interval:

    WAVE 1  -> minimum wavelength
    WAVE 11 -> maximum wavelength

A set of 100 equally spaced wavelengths is generated between these two
limits for the active configuration.

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

For each sampled wavelength, a set of points from the normalized pupil
is traced to the selected optical surface using `zGetTrace()`. The pupil
sampling is configurable and can be generated using either a random or
a hexapolar distribution. The number of sampled pupil points is also
defined by the user. The pupil-sampling strategy depends on the analysis stage. Initial spectral
reconstruction experiments use random normalized-pupil samples, typically
with 100 pupil points for each of the 100 sampled wavelengths. Extended-source
analyses use hexapolar sampling to provide repeatable coverage
of both the source geometry and the pupil boundary.

In the robust full-format fiber calculation, the `100 µm` circular fiber is
represented by 91 hexapolar field positions and the pupil by 91 hexapolar
rays. Each diffraction order is sampled at 100 wavelengths. A reduced
five-field source representation is also implemented and benchmarked against
the dense 91-field calculation. For the initial multi-configuration spectral reconstruction, 
the tracing sequence can be represented conceptually as:

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

For the extended-source fiber analysis, an additional sampling level is
introduced because rays must be traced from multiple positions across the
physical fiber:

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

where `N = 91` for the robust hexapolar fiber representation and `N = 5`
for the reduced sampling experiment.

After the sampled wavelengths of a configuration have been traced, the
original `WAVE 1` value is restored in the MCE and the model is updated
again with `zGetUpdate()` before proceeding to the next configuration.
This restoration prevents temporary wavelength modifications from
affecting subsequent configurations or analyses.

### Complete order tracing

The workflow extends the original configuration-based
procedure to every integer diffraction order from `m = 60` to
`m = 144`. The spectral definition of these orders is not recalculated during the
dense ray-tracing stage. Instead, `EMAR_full_orders.py` reads the
validated `echelle_orders_60_144.csv` table produced by
`EMAR_missing_orders.py`.

For each order:

1. Read the order number, reconstruction status, and 11 reference
   wavelengths from the CSV table.
2. Use the first and last reference wavelengths to define the spectral
   interval of the order.
3. Generate a dense set of equally spaced wavelengths across this
   interval.
4. Select the appropriate original Zemax configuration as the optical
   template.
5. Temporarily assign the required diffraction order and wavelength
   through the MCE.
6. Trace the common normalized-pupil sample to the selected optical
   surface.
7. Restore the original diffraction-order and wavelength values before
   continuing.

The procedure therefore reconstructs the spectral parameters required
by Zemax rather than interpolating detector coordinates. The X-Y
positions of the reconstructed orders are obtained directly from ray
tracing through the optical model. For the current dense analysis, 100 
wavelengths are sampled for each of the 85 diffraction orders, 
with 100 normalized pupil rays traced at each wavelength.

## Current Analysis Status

The current workflow provides a complete single-fiber analysis of
the EMAR echelle format from `m = 60` to `m = 144`. The robust 91-field
calculation traces a circular `100 µm` input fiber across all 85 integer
diffraction orders and provides wavelength-dependent projected fiber
dimensions, spectral dispersion, resolving power, and adjacent-order free
separation. A reduced five-field source representation was benchmarked against
the 91-field reference, reproducing the wavelength-dependent behavior of the
projected fiber and resolving power at substantially lower computational cost.
The dense 91-field solution is retained as the reference for geometrical
envelope and order-separation analyses.

Across the complete spectral format, the projected single-fiber envelopes of
adjacent echelle orders remain geometrically separated, with the limiting
edge-to-edge separation of approximately `141 µm` occurring between orders
`m = 60` and `m = 61`. This single-fiber analysis establishes the geometrical
reference for the next development stage: tracing multiple fiber images within
each echelle order and evaluating their compatibility with a dual-beam
polarimetric configuration