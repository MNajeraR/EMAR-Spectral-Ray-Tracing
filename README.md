# EMAR Spectral Ray-Tracing Analysis
## Overview

This repository contains a Python-based ray-tracing and spectral-analysis
framework developed for the Espectrógrafo Mexicano de Alta Resolución (EMAR)
optical system. The code interfaces with Zemax OpticStudio through PyZDDE and
provides tools for normalized-pupil sampling, multi-configuration ray tracing,
wavelength and diffraction-order control through the Zemax Multi-Configuration
Editor (MCE), spectral-dispersion analysis, and geometrical characterization of
extended input sources. The interaction between Python and the Zemax
multi-configuration model, including configuration selection, temporary MCE
modification, and ray tracing, is described in detail in the
[Zemax Model Interaction](#zemax-model-interaction) section. The Zemax model
contains 13 representative configurations associated with echelle diffraction
orders spanning from `m = 60` to `m = 144`, with each configuration defining
11 reference wavelengths through the MCE. The framework was developed
progressively, beginning with individual pupil and wavelength tests and
extending to the reconstruction and analysis of the complete integer
echelle-order sequence.

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

A circular `100 µm` fiber is used to establish a single-fiber reference across
all 85 integer diffraction orders, providing the basis for subsequent studies
of multiple fiber traces and polarimetric configurations. The original Zemax
optical model is intentionally not distributed with this repository; users
with authorized access to the EMAR model can place their local copy in the
`Zemax/` directory and use the analysis scripts without modifying the
repository structure.

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
   centroid and its trajectory across the image plane. The centroid
   displacement is then used to derive the local linear and reciprocal
   spectral dispersion as functions of wavelength.

9. `09_single_order_spectral_resolution.py` extends the point-source analysis
   of `m = 102` from dispersion to geometrical spectral resolution. At each
   wavelength, the monochromatic ray footprint produced by the point source is
   projected onto the local direction of spectral dispersion. The projected
   footprint width provides a geometrical estimate of the wavelength interval
   associated with the monochromatic image and therefore of the
   wavelength-dependent resolving power. This experiment establishes the
   resolution-analysis procedure before replacing the point source with a
   finite physical entrance source.

10. `10_theoretical_slit.py` introduces a finite rectangular source at the
    spectrograph input to replace the point-field approximation used in the
    preceding experiments. A theoretical `50 × 15 µm` slit is mapped onto the
    corresponding Zemax field coordinates and sampled over its physical
    extent. Rays from each slit position are traced through the pupil and onto
    the detector. The resulting finite-source footprints are projected along
    the local spectral direction to determine the projected slit width,
    spectral purity, and geometrical resolving power across order `m = 102`.

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

12. `12_full_spectral_format.py` extends the validated `100 µm` fiber
    calculation from the representative order to the complete reconstructed
    echelle format. All integer orders from `m = 60` to `m = 144` are sampled
    at 100 wavelengths. At every wavelength, 91 fiber positions are combined
    with 91 pupil rays, retaining the dense source and pupil
    sampling established in the single-order experiment. The resulting
    detector coordinates are used to calculate the fiber centroid, X-Y
    envelope, linear spectral dispersion, projected spectral width, and
    resolving power throughout all 85 diffraction orders. The same results
    also provide the projected fiber envelopes required for subsequent
    analysis of the physical spacing between neighboring spectral orders.

13. `13_five_field_fiber_spectrum.py` evaluates whether the computational
    cost of the complete-format fiber calculation can be reduced by replacing
    the 91 source positions with five representative fiber fields. The reduced
    geometry samples the fiber center together with the positive and negative
    extrema along X and Y, while retaining the 91-ray hexapolar pupil and
    100 wavelengths per order. The complete sequence from `m = 60` to
    `m = 144` is traced using this representation, providing a direct
    alternative to the dense 91-field calculation without changing the
    wavelength or pupil sampling.

14. `14_compare_fiber_sampling.py` benchmarks the five-field representation
    against the robust 91-field reference without performing additional Zemax
    ray tracing. The wavelength-by-wavelength results of both calculations are
    compared across all 85 orders to quantify absolute and relative
    differences in projected fiber width and resolving power. Per-order
    statistics are also generated to determine whether the reduced source
    representation preserves the spectral behavior observed with dense
    sampling. This comparison establishes the 91-field calculation as the
    reference for analyses that depend on the extrema of the projected fiber
    envelope while providing a faster reduced representation for exploratory
    calculations.

15. `15_order_separation.py` uses the robust 91-field results to quantify the
    free geometrical space between projected fiber envelopes in adjacent
    echelle orders. For every neighboring pair from `m = 60–61` through
    `m = 143–144`, the two order envelopes are evaluated over their common
    detector-X interval. Their Y boundaries are interpolated at common X
    coordinates, and the edge-to-edge free separation is calculated between
    the upper boundary of the lower trace and the lower boundary of the upper
    trace. Minimum, maximum, mean, and median separations are determined for
    every adjacent-order pair together with the detector-X locations of the
    extrema. Representative separation profiles are also evaluated across the
    detector, providing the geometrical reference required for subsequent
    multi-fiber and polarimetric configurations.


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

For each sampled wavelength, a set of points from the normalized pupil is
traced to the selected optical surface using `zGetTrace()`. The pupil is
sampled using a random distribution within the normalized circular aperture,
with the number of rays defined according to the analysis. Initial spectral
reconstruction experiments typically use 100 pupil rays for each of the
100 sampled wavelengths.

For the extended-source analysis, an additional sampling level is introduced
to represent the physical extent of the input source. In the robust
full-format fiber calculation, the `100 µm` circular fiber is represented by
91 field positions arranged in a five-ring hexapolar distribution. A
normalized-pupil sample is independently traced from each fiber position at
each wavelength. Each diffraction order is sampled at 100 wavelengths. A
reduced five-field representation of the fiber is also implemented and
benchmarked against the dense 91-field calculation.

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

where `N = 91` for the robust hexapolar representation of the circular fiber
and `N = 5` for the reduced fiber-field sampling experiment.

After the sampled wavelengths of a configuration have been traced, the
original `WAVE 1` value is restored in the MCE and the model is updated again
with `zGetUpdate()` before proceeding to the next configuration. This
restoration prevents temporary wavelength modifications from affecting
subsequent configurations or analyses.

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
adjacent echelle orders remain separated, with the limiting
edge-to-edge separation of approximately `141 µm` occurring between orders
`m = 60` and `m = 61`. This single-fiber analysis establishes the geometrical
reference for the next development stage: tracing multiple fiber images within
each echelle order and evaluating their compatibility with a dual-beam
polarimetric configuration