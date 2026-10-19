# EMAR Spectral Ray-Tracing Analysis

[![DOI](https://zenodo.org/badge/1363089025.svg)](https://doi.org/10.5281/zenodo.23198727)

## Overview

This repository contains a Python-based ray-tracing and spectral-analysis
framework developed for the Espectrógrafo Mexicano de Alta Resolución (EMAR),
a fiber-fed, high-resolution echelle spectrograph. The framework interfaces
with Zemax OpticStudio through PyZDDE, providing computational tools for
investigating the optical performance and spectral characteristics of the
instrument. The Zemax optical model contains 13 representative configurations associated
with echelle diffraction orders spanning from `m = 60` to `m = 144`, each
defining 11 reference wavelengths through the Multi-Configuration Editor
(MCE). The framework was developed progressively, beginning with individual
ray-tracing experiments and extending to the reconstruction of the complete
sequence of 85 integer diffraction orders. A circular `100 µm` fiber is used
as the reference input source for characterizing the resulting spectral
format on the detector.

The current workflow includes:

- Configurable field and normalized-pupil sampling using random or hexapolar
  distributions.
- Multi-configuration ray tracing and dense wavelength sampling across
  individual and complete echelle orders.
- Visualization of ray footprints, spectral centroid traces, and the
  complete detector spectral format.
- Calculation of linear dispersion and resolving power.
- Characterization of extended input sources, including finite slits and
  circular fibers.
- Comparison of fiber-field sampling strategies and evaluation of
  inter-order separation.
- Modeling of spatially separated polarimetric images and their mapping
  onto the detector.
- Geometric clearance and overlap analysis between adjacent diffraction
  orders.
- Reconstruction and visualization of the complete polarimetric spectral
  format across all 85 integer echelle orders.

The framework has subsequently been extended to investigate the geometric
requirements of a potential spectropolarimetric configuration. In this
approach, the two polarization components are represented as spatially
separated images of the same input fiber, following the action of a
polarization-splitting element. Their relative displacement at the
spectrograph entrance is treated as a geometrical parameter, allowing
the available detector space and potential overlap between adjacent
spectral orders to be evaluated without explicitly modeling the
polarization-splitting optics. The original Zemax optical model is intentionally not distributed with
this repository. Users with authorized access to the EMAR model should
place the file `WP - Con prismas diseñados - camara.zmx` in the `Zemax/`
directory to reproduce the analyses without modifying the repository
structure.

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
reusable functions, experiments, results, and the
local Zemax optical model.

```text
EMAR/
├── README.md
├── requirements.txt
├── CITATION.cff
├── LICENSE
├── .gitignore
│
├── scripts/
│   ├── EMAR_All_conf_surf_26.py
│   ├── EMAR_All_confs_surf54.py
│   ├── EMAR_missing_orders.py
│   ├── EMAR_full_orders.py
│   └── EMAR_polarimetric_all_orders.py
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
│   ├── 12_full_fiber_spectrum.py
│   ├── 13_five_field_fiber_spectrum.py
│   ├── 14_compare_fiber_sampling.py
│   ├── 15_order_separation.py
│   ├── 16_polarimetric_image_geometry.py
│   ├── 17_polarimetric_field_mapping.py
│   ├── 18_polarimetric_separation_scan.py
│   └── 19_polarimetric_adjacent_orders.py
│
├── utils/
│   ├── __init__.py
│   ├── emar_utils.py
│   └── emar_plots.py
│
├── Results/
│   ├── echelle_orders_60_144.csv
│   ├── ...
│   ├── polarimetric_envelopes_m60_61_d*.csv
│   ├── polarimetric_clearance_m60_61_d*.csv
│   ├── polarimetric_first_contact_m60_61_*
│   └── polarimetric_all_orders_d350um/
│       ├── envelopes/
│       │   ├── m60_d350um_100w_9fields.csv
│       │   ├── ...
│       │   └── m144_d350um_100w_9fields.csv
│       ├── 01_full_spectral_format.png
│       ├── 02_orders_60_61_zoom.png
│       ├── 03_detector_true_aspect.png
│       ├── adjacent_order_clearance.csv
│       └── polarimetric_envelopes_m60_144.csv
│
└── Zemax/
    └── .gitkeep
```

### `scripts/`

Contains the main analysis scripts used to study the EMAR optical
system at selected stages of the spectrograph.

- `EMAR_All_confs_surf_26.py` traces the 13 original Zemax
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

- `EMAR_polarimetric_all_orders.py` extends the complete echelle-order
  analysis to a potential spectropolarimetric configuration. Two spatially
  separated images of a single circular fiber are traced across all 85
  diffraction orders (`m = 60–144`). The script reconstructs their detector
  envelopes, evaluates geometric clearance and overlap between adjacent
  orders, and generates three complementary visualizations: the complete
  polarimetric spectral format, a detailed view of orders 60 and 61, and
  a detector representation preserving the physical aspect ratio.
  Ray-tracing results are cached to avoid unnecessary repeated calculations.

### `experiments/`

Contains the individual experiments used to develop and validate the
computational procedures implemented in the main analysis scripts.
Rather than executing the complete spectral analysis, each experiment
focuses on a specific numerical or geometrical problem, allowing the
underlying assumptions, sampling strategies, and ray-tracing results
to be examined independently. The sequence also documents how the
methods were refined before being extended to more complex
configurations.

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

12. `12_full_fiber_spectrum.py` extends the `100 µm` fiber analysis from
    `m = 102` to all integer orders from `m = 60` to `m = 144`. Each order is
    sampled at 100 wavelengths, retaining the 91 fiber positions and 91 pupil
    rays established in the preceding experiment. The resulting detector
    coordinates are used to calculate the fiber centroid, X-Y envelope, linear
    dispersion, projected spectral width, and resolving power across
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

16. `16_polarimetric_image_geometry.py` introduces the geometrical
    representation of two spatially separated polarimetric images at
    the spectrograph entrance. The experiment establishes the
    relationship between their center-to-center displacement and
    edge-to-edge separation, considering the finite diameter of the
    input fiber. Each circular image is represented by its center
    and eight boundary points, providing a consistent geometrical
    description for the subsequent ray-tracing experiments.

17. `17_polarimetric_field_mapping.py` investigates the transformation
    between physical coordinates at the spectrograph entrance and
    the normalized field coordinates required by Zemax. The
    geometrical representation established in the preceding
    experiment is translated into field positions and propagated
    through the optical system. The resulting detector footprints
    are examined to determine how the relative position of the
    input images is transformed by the spectrograph.

18. `18_polarimetric_separation_scan.py` examines how the separation
    between two polarimetric images affects their projected footprints
    on the detector for diffraction order `m = 60` at a single
    reference wavelength. Center-to-center separations from 110 to
    115 µm are evaluated in increments of 0.5 µm, using nine field
    points per image and 91 hexapolar pupil rays per field.
    Convex-hull representations of the resulting ray distributions
    are used to calculate the overlap area and minimum distance
    between the projected images. This experiment establishes how
    their geometrical overlap varies with entrance separation before
    extending the analysis to neighboring diffraction orders.

19. `19_polarimetric_adjacent_orders.py` extends the geometrical
    analysis to the neighboring diffraction orders `m = 60` and
    `m = 61`, identified in Experiment 15 as the pair with the
    minimum inter-order separation across the complete echelle
    format. This pair provides a reference case for investigating
    the available detector space when two polarimetric images are
    introduced. Both images are traced across 100 wavelengths per
    order, using nine field points and 91 pupil rays per field.
    The detector envelopes are reconstructed from the sampled
    ray distributions and interpolated over a common detector-X
    interval. The four possible image combinations between the
    two orders are evaluated to determine their minimum geometric
    clearance and identify regions of envelope overlap. The
    resulting procedure provides the basis for extending the
    clearance analysis to all adjacent diffraction-order pairs.
    
### `utils/`

Provides the shared computational and visualization functions used
throughout the repository, separating reusable operations from the
individual analysis workflows.

- `emar_utils.py` contains project paths, normalized-pupil sampling
  routines, PyZDDE ray-tracing functions, and multi-configuration
  tracing utilities. It also provides functions for tracing original
  and reconstructed echelle diffraction orders while preserving and
  restoring the original MCE state.

- `emar_plots.py` contains visualization routines for pupil sampling,
  multi-configuration ray footprints, and the complete echelle spectral
  format. Results can be displayed as individual ray distributions,
  spectral centroid traces, or detector-style representations.
  Echelle-order visualizations use a continuous color scale, with
  lower diffraction orders represented toward red and higher orders
  toward blue.

### `Results/`

Contains numerical datasets and analysis products generated by the
main scripts and validation experiments. The master spectral definition, 
`echelle_orders_60_144.csv`, identifies
each diffraction order as `ORIGINAL` or `NEW` and provides the 11
reference wavelengths used to construct the corresponding spectral
sampling grids.

Single-order analysis products include:

- `order_102_dispersion.csv`, containing the wavelength-dependent
  detector centroid and linear spectral dispersion for `m = 102`;
- `order_102_rays.csv`, containing the corresponding individual
  detector ray coordinates;
- `order_102_theoretical_slit_resolution.csv`, containing the
  spectral-resolution analysis for the finite theoretical slit; and
- `order_102_fiber_resolution.csv`, containing the projected fiber
  dimensions and resolving-power analysis for `m = 102`.

Complete-format fiber products include:

- `full_orders_fiber_resolution.csv`, containing the 91-field fiber
  analysis across all 85 diffraction orders;
- `full_orders_fiber_resolution_5fields.csv`, containing the
  corresponding reduced five-field calculation;
- `fiber_sampling_comparison.csv`, containing the wavelength-dependent
  comparison between the two sampling strategies;
- `fiber_sampling_statistics_by_order.csv`, containing the comparison
  statistics summarized by diffraction order; and
- `order_separation_91fields.csv`, containing the free-separation
  statistics between adjacent echelle-order fiber envelopes.

Polarimetric analysis products include:

- `polarimetric_envelopes_m60_61_d*_100w_9fields.csv`, containing
  the detector envelopes of the two polarimetric images for orders
  `m = 60` and `m = 61` at different input-image separations.
  Each calculation uses 100 wavelengths per order and nine field
  points per image.

- `polarimetric_clearance_m60_61_d*.csv`, containing the geometric
  clearance measurements between adjacent-order polarimetric
  envelopes for the corresponding input-image separations.

- `polarimetric_first_contact_m60_61_d500um.csv`, containing
  results associated with the investigation of the first
  geometrical contact between the polarimetric envelopes
  of orders `m = 60` and `m = 61`.

- `polarimetric_all_orders_d350um/`, containing the complete
  polarimetric spectral-format analysis for a center-to-center
  input-image separation of 350 µm.

Within `polarimetric_all_orders_d350um/`, the numerical results
are organized as follows:

- `envelopes/` contains 85 individual CSV files, one for each
  diffraction order from `m = 60` to `m = 144`. Each file follows
  the naming convention `m{order}_d350um_100w_9fields.csv` and
  stores the corresponding detector-envelope data. These files
  allow previously computed order envelopes to be reused without
  repeating the Zemax ray tracing.

- `polarimetric_envelopes_m60_144.csv` consolidates the
  reconstructed polarimetric envelopes across the complete
  diffraction-order sequence.

- `adjacent_order_clearance.csv` contains the geometric clearance
  measurements between neighboring polarimetric order envelopes
  across the complete spectral format.

The directory also contains three complementary visualizations:

- `01_full_spectral_format.png`, showing the complete polarimetric
  spectral format across all 85 diffraction orders;
- `02_orders_60_61_zoom.png`, providing a detailed view of the
  neighboring orders identified as having the smallest
  inter-order separation; and
- `03_detector_true_aspect.png`, displaying the complete spectral
  format while preserving the physical detector aspect ratio.

The analysis scripts also generate diagnostic figures and intermediate
ray-tracing products. Selected figures and numerical datasets are
included in the repository, while other intermediate outputs may
remain local depending on the analysis configuration.

### `Zemax/`

Provides the expected local location for the EMAR Zemax optical model.
The `.gitkeep` file preserves this directory in Git while its optical
model contents remain excluded through `.gitignore`.

## Running the Analysis

Before running the EMAR analysis scripts, Zemax OpticStudio must be
open and the PyZDDE connection must be properly configured.
Installation instructions and connection examples are available in
the [official PyZDDE repository](https://github.com/xzos/PyZDDE)
and can be used to verify the setup before executing the EMAR
workflow. A local copy of the EMAR Zemax model must also be
available in the `Zemax/` directory. Once the setup has been verified,
`experiments/02_single_config_original_waves.py` provides a minimal
starting point by selecting a single configuration, reading its
reference wavelengths, and tracing rays to the selected optical
surface. The numbered experiments can then be executed individually to
examine the development and validation of specific computational
procedures. The scripts in `scripts/` provide the corresponding
complete-format analyses, including the reconstruction of all
integer echelle orders and the evaluation of a potential
spectropolarimetric configuration.

## Zemax Model Interaction

The computational workflow relies on reading and temporarily modifying
the Multi-Configuration Editor (MCE) parameters to perform the required
ray-tracing analyses. This approach allows the original spectral
definitions to be extended while preserving the optical model through
the restoration of its initial configuration and parameters.

### Configuration selection

The original Zemax configurations can be selected directly in Python:

```python
configs = range(1, 14)
```

The script iterates over the selected configurations using
`zSetConfig()`, followed by `zGetUpdate()` to refresh the optical
system. For reconstructed diffraction orders, the corresponding
original configuration is selected as a reference, and the required
echelle order and wavelength parameters are temporarily modified
through the MCE. This allows intermediate integer orders to be traced
using the optical model rather than interpolating detector coordinates.

### Reading the MCE

Once a configuration is selected, `zGetMulticon()` reads the
corresponding MCE entries, including the echelle diffraction order
and its reference wavelengths. The first and last reference
wavelengths define the spectral interval used to generate a grid
of 100 equally spaced wavelengths for dense spectral sampling:

```text
WAVE 1  -> minimum wavelength
WAVE 11 -> maximum wavelength
```

### Temporary wavelength modification

Each sampled wavelength is temporarily assigned to `WAVE 1` using
`zSetMulticon()`, and the optical model is refreshed with
`zGetUpdate()` before tracing the corresponding rays. When
reconstructing intermediate diffraction orders, the echelle order
parameter is also modified through the MCE. The basic sequence is:

```text
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
```

After completing the ray tracing, the original MCE values are
restored and the optical model is updated, preventing temporary
modifications from affecting subsequent calculations.

### Ray tracing

Ray tracing is performed using `zGetTrace()`, with the field
coordinates, normalized-pupil coordinates, wavelength, and target
surface specified according to the analysis. The sampling strategy
depends on the input-source geometry:

- Point-source analyses use a single field position and random
  normalized-pupil sampling.
- Finite-slit analyses sample the source extent and pupil using
  random distributions.
- Circular-fiber analyses use deterministic hexapolar sampling
  of the source field and normalized pupil.

For the full-format fiber calculation, the circular fiber is
represented by 91 field positions arranged in a five-ring
hexapolar distribution, each traced through 91 normalized-pupil
positions, producing 8,281 rays per wavelength. A reduced sampling
strategy uses five representative fiber positions corresponding
to the center and four extrema along X and Y, while retaining
the 91-point pupil distribution. This allows the effect of
reducing the source-field sampling to be evaluated independently.

The general tracing sequence for a finite input source is:

```text
Diffraction order
    |
    ├──> wavelength 1
    │       ├──> field position 1 -> pupil sample -> ray tracing
    │       ├──> field position 2 -> pupil sample -> ray tracing
    │       ├── ...
    │       └──> field position N -> pupil sample -> ray tracing
    │
    ├──> wavelength 2
    │       └──> field positions -> pupil samples -> ray tracing
    │
    ├── ...
    │
    └──> wavelength 100
            └──> field positions -> pupil samples -> ray tracing
```

Here, `N = 91` for the full fiber sampling and `N = 5` for
the reduced sampling strategy.

For the spectropolarimetric analysis, the source-field sampling
is extended to two spatially displaced circular images, each
represented by nine field positions: its center and eight equally
spaced points along its circumference. The physical input
coordinates are converted into normalized Zemax field coordinates,
and 91 hexapolar pupil rays are traced for each field position.
The resulting tracing sequence becomes:

```text
Diffraction order
    |
    ├──> wavelength 1
    │       ├──> polarimetric image 1
    │       │       └──> 9 field positions -> 91 pupil rays each
    │       │
    │       └──> polarimetric image 2
    │               └──> 9 field positions -> 91 pupil rays each
    │
    ├──> wavelength 2
    │       └──> both images -> field sampling -> pupil sampling
    │
    ├── ...
    │
    └──> wavelength 100
            └──> both images -> field sampling -> pupil sampling
```

This procedure produces 1,638 traced rays per wavelength,
corresponding to two images, nine field positions per image,
and 91 pupil rays per field. The detector coordinates are
subsequently grouped by diffraction order and polarimetric
image to reconstruct their spectral envelopes.

## Current Analysis Status

The current workflow provides a complete single-fiber analysis of the EMAR
echelle format from `m = 60` to `m = 144`. A circular `100 µm` input fiber is
traced across all 85 integer diffraction orders using 91-field hexapolar
sampling, providing wavelength-dependent projected fiber dimensions, lineal
dispersion, resolving power, and adjacent-order separation. The calculated
resolving power ranges from approximately `28,403` to `54,863` across the
complete spectral format. A reduced five-field sampling strategy was
benchmarked against the full calculation, reproducing the wavelength-dependent
behavior of the projected fiber and resolving power at substantially lower
computational cost. The 91-field results are retained for fiber-envelope
and order-separation analyses. Across the complete spectral format, adjacent
single-fiber envelopes remain separated, with a minimum edge-to-edge distance
of approximately `141 µm` between orders `m = 60` and `m = 61`.

Building on these results, the framework has been extended to investigate
the geometric feasibility of a potential dual-beam spectropolarimetric
configuration. Two spatially separated images of the same `100 µm` input
fiber are propagated through the optical model, with their center-to-center
separation treated as an adjustable parameter. The analysis initially
focuses on diffraction orders `m = 60` and `m = 61`, identified as the
adjacent pair with the smallest inter-order separation. Their detector
envelopes are reconstructed using 100 wavelengths per order, nine field
positions per image, and 91 pupil rays per field. Different input-image
separations are evaluated to characterize the geometric clearance and
potential overlap between the polarimetric envelopes of these orders.

The procedure is subsequently extended to all 85 integer diffraction
orders using `EMAR_polarimetric_all_orders.py`. For a center-to-center
input-image separation of `350 µm`, the complete polarimetric spectral
format is reconstructed, and the geometric clearance between adjacent
orders is evaluated. Across 336 comparisons, corresponding to the four
possible polarimetric-image combinations for each pair of neighboring
orders, no envelope overlaps are detected. The minimum geometric
clearance is approximately `44.45 µm`, occurring between the upper image
of order `m = 60` and the lower image of order `m = 61`. These results
indicate that the investigated configuration maintains positive
inter-order clearance within the adopted geometrical model, providing
a basis for further evaluation of the detector-space requirements
without explicitly modeling the polarization-splitting optics.
