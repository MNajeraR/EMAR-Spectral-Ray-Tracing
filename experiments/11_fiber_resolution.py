# -*- coding: utf-8 -*-

"""
11_fiber_resolution.py
======================

Geometrical spectral resolution of EMAR using a
100 um circular fiber.

The circular fiber and the optical pupil are sampled using
deterministic 5-ring hexapolar distributions.

For every wavelength of diffraction order m = 102:

    1. Sample the 100 um circular fiber with 91 field points.
    2. Trace 91 pupil rays for every fiber point.
    3. Measure the complete geometrical fiber image at S54.
    4. Calculate:

           Delta X = Xmax - Xmin
           Delta Y = Ymax - Ymin

    5. Since the spectral dispersion is predominantly along X:

           omega' = Delta X

    6. Read the previously calculated reciprocal dispersion
       from order_102_dispersion.csv.

    7. Calculate:

           delta_lambda = omega' * d(lambda)/dX

           R = lambda / delta_lambda

The detector pixel size is NOT required for this calculation.

Fiber sampling:
    5 hexapolar rings -> 91 points

Pupil sampling:
    5 hexapolar rings -> 91 rays

Total rays per wavelength:
    91 x 91 = 8281

Author
------
Morgan Rhaí Nájera Roa
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyzdde.zdde as pyz


##############################################################
# Project paths
##############################################################

PROJECT_DIR = Path(__file__).resolve().parents[1]

if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(
        0,
        str(PROJECT_DIR)
    )


from utils import emar_utils


##############################################################
# Progress bar
##############################################################

def print_progress(
    current,
    total,
    bar_length=40
):
    """
    Display a simple terminal progress bar.
    """

    fraction = current / total

    filled = int(
        bar_length * fraction
    )

    bar = (
        "=" * filled
        + "-" * (bar_length - filled)
    )

    percent = (
        100.0 * fraction
    )

    print(
        f"\r[{bar}] "
        f"{percent:6.2f}% "
        f"({current}/{total})",
        end="",
        flush=True
    )

    if current == total:
        print()


##############################################################
# Experiment parameters
##############################################################

target_order = 102

surf_detector = 54


##############################################################
# Fiber geometry
##############################################################

fiber_diameter_um = 100.0


##############################################################
# Hexapolar sampling
##############################################################

n_fiber_rings = 5

n_pupil_rings = 5


##############################################################
# Calibrated field limits
##############################################################

# Previous calibration at surface 6:
#
# Hx = +/- 0.286480
# Hy = +/- 0.285400
#
# These limits represent a circular fiber with a physical
# radius of approximately 50 um at the input plane.

hx_max = 0.286480

hy_max = 0.285400


##############################################################
# MCE information
##############################################################

row_order = 1

row_wave1 = 2


##############################################################
# Configuration associated with order 102
##############################################################

base_config = (
    (target_order - 60) // 7 + 1
)


##############################################################
# Files
##############################################################

zemax_file = (
    PROJECT_DIR
    / "Zemax"
    / "WP - Con prismas diseñados - camara - theoretical slit.zmx"
)


dispersion_file = (
    PROJECT_DIR
    / "Results"
    / "order_102_dispersion.csv"
)


output_file = (
    PROJECT_DIR
    / "Results"
    / "order_102_fiber_resolution.csv"
)


##############################################################
# Read previously calculated dispersion
##############################################################

dispersion_data = pd.read_csv(
    dispersion_file
)


##############################################################
# Check required columns
##############################################################

required_columns = [
    "wavelength_nm",
    "reciprocal_dispersion_nm_per_mm"
]


for column in required_columns:

    if column not in dispersion_data.columns:

        raise KeyError(
            f"Column '{column}' not found in "
            f"{dispersion_file.name}"
        )


##############################################################
# Wavelengths
##############################################################

wavelengths_nm = (
    dispersion_data[
        "wavelength_nm"
    ].to_numpy(
        dtype=float
    )
)


wavelengths_um = (
    wavelengths_nm
    / 1000.0
)


##############################################################
# Reciprocal linear dispersion
##############################################################

# Units:
#
#       nm / mm

reciprocal_dispersion = (
    dispersion_data[
        "reciprocal_dispersion_nm_per_mm"
    ].to_numpy(
        dtype=float
    )
)


##############################################################
# Generate normalized circular-fiber sampling
##############################################################

# sample_hexapolar_pupil() generates points over a normalized
# unit disk.
#
# Here we reuse that geometrical distribution for the fiber.

fiber_x_norm, fiber_y_norm = (
    emar_utils.sample_hexapolar_pupil(
        n_rings=n_fiber_rings,
        include_center=True
    )
)


##############################################################
# Convert normalized fiber coordinates to field coordinates
##############################################################

hx_fiber = (
    hx_max
    * fiber_x_norm
)

hy_fiber = (
    hy_max
    * fiber_y_norm
)


##############################################################
# Optical pupil sampling
##############################################################

px, py = (
    emar_utils.sample_hexapolar_pupil(
        n_rings=n_pupil_rings,
        include_center=True
    )
)


##############################################################
# Sampling information
##############################################################

n_fiber_points = len(
    hx_fiber
)

n_pupil_rays = len(
    px
)

n_rays_per_wavelength = (
    n_fiber_points
    * n_pupil_rays
)

n_wavelengths = len(
    wavelengths_um
)

n_total_rays = (
    n_rays_per_wavelength
    * n_wavelengths
)


##############################################################
# Experiment information
##############################################################

print()

print(
    "EMAR 100 um fiber spectral resolution"
)

print(
    "-------------------------------------"
)

print(
    f"Order:                    m = {target_order}"
)

print(
    f"Base config:              {base_config}"
)

print(
    f"Detector surface:         {surf_detector}"
)

print()

print(
    f"Fiber diameter:           {fiber_diameter_um:.1f} um"
)

print(
    f"Hx limit:                 +/- {hx_max:.6f}"
)

print(
    f"Hy limit:                 +/- {hy_max:.6f}"
)

print()

print(
    f"Fiber rings:              {n_fiber_rings}"
)

print(
    f"Fiber points:             {n_fiber_points}"
)

print(
    f"Pupil rings:              {n_pupil_rings}"
)

print(
    f"Pupil rays:               {n_pupil_rays}"
)

print()

print(
    f"Wavelengths:              {n_wavelengths}"
)

print(
    f"Rays per wavelength:      {n_rays_per_wavelength}"
)

print(
    f"Maximum total rays:       {n_total_rays}"
)

print()

print(
    f"Wavelength range:         "
    f"{wavelengths_nm.min():.6f} - "
    f"{wavelengths_nm.max():.6f} nm"
)

print()


##############################################################
# Result arrays
##############################################################

x_centroid = np.full(
    n_wavelengths,
    np.nan
)

y_centroid = np.full(
    n_wavelengths,
    np.nan
)


x_minimum = np.full(
    n_wavelengths,
    np.nan
)

x_maximum = np.full(
    n_wavelengths,
    np.nan
)


y_minimum = np.full(
    n_wavelengths,
    np.nan
)

y_maximum = np.full(
    n_wavelengths,
    np.nan
)


delta_x_um = np.full(
    n_wavelengths,
    np.nan
)

delta_y_um = np.full(
    n_wavelengths,
    np.nan
)


valid_ray_count = np.zeros(
    n_wavelengths,
    dtype=int
)


##############################################################
# Open Zemax
##############################################################

ln = pyz.createLink()

ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Trace fiber for every wavelength
##############################################################

try:

    print(
        "Tracing fiber images..."
    )

    print()


    for wave_index, wavelength_um in enumerate(
        wavelengths_um
    ):

        ######################################################
        # Complete detector footprint for this wavelength
        ######################################################

        wavelength_x = []

        wavelength_y = []


        ######################################################
        # Loop over fiber points
        ######################################################

        for hx, hy in zip(
            hx_fiber,
            hy_fiber
        ):

            result = emar_utils.trace_all_order(
                ln=ln,
                wavelengths=np.array([
                    wavelength_um
                ]),
                px=px,
                py=py,
                surf=surf_detector,
                target_order=target_order,
                base_config=base_config,
                row_order=row_order,
                row_wave1=row_wave1,
                hx=float(hx),
                hy=float(hy)
            )


            ##################################################
            # Detector footprint
            ##################################################

            spot = next(
                iter(
                    result["spots"].values()
                )
            )


            x = np.asarray(
                spot["x"],
                dtype=float
            )

            y = np.asarray(
                spot["y"],
                dtype=float
            )


            ##################################################
            # Store valid rays
            ##################################################

            if len(x) > 0:

                wavelength_x.extend(
                    x
                )

                wavelength_y.extend(
                    y
                )


        ######################################################
        # Convert to arrays
        ######################################################

        wavelength_x = np.asarray(
            wavelength_x,
            dtype=float
        )

        wavelength_y = np.asarray(
            wavelength_y,
            dtype=float
        )


        ######################################################
        # Number of valid rays
        ######################################################

        valid_ray_count[
            wave_index
        ] = len(
            wavelength_x
        )


        ######################################################
        # Measure geometrical fiber image
        ######################################################

        if len(wavelength_x) > 0:

            ##############################################
            # Centroid
            ##############################################

            x_centroid[
                wave_index
            ] = np.mean(
                wavelength_x
            )

            y_centroid[
                wave_index
            ] = np.mean(
                wavelength_y
            )


            ##############################################
            # Limits
            ##############################################

            x_minimum[
                wave_index
            ] = np.min(
                wavelength_x
            )

            x_maximum[
                wave_index
            ] = np.max(
                wavelength_x
            )


            y_minimum[
                wave_index
            ] = np.min(
                wavelength_y
            )

            y_maximum[
                wave_index
            ] = np.max(
                wavelength_y
            )


            ##############################################
            # Full image amplitudes
            ##############################################

            delta_x_um[
                wave_index
            ] = (
                x_maximum[
                    wave_index
                ]
                -
                x_minimum[
                    wave_index
                ]
            ) * 1000.0


            delta_y_um[
                wave_index
            ] = (
                y_maximum[
                    wave_index
                ]
                -
                y_minimum[
                    wave_index
                ]
            ) * 1000.0


        ######################################################
        # Progress bar
        ######################################################

        print_progress(
            current=wave_index + 1,
            total=n_wavelengths
        )


finally:

    ##########################################################
    # Close Zemax
    ##########################################################

    pyz.closeLink()


##############################################################
# Fiber image width used for spectral resolution
##############################################################

# For order 102 the spectral trace is predominantly oriented
# along detector X.
#
# Therefore, following the method being reproduced:
#
#       omega' = Delta X

omega_prime_um = (
    delta_x_um.copy()
)


omega_prime_mm = (
    omega_prime_um
    / 1000.0
)


##############################################################
# Spectral purity
##############################################################

# omega':
#
#       mm
#
# reciprocal dispersion:
#
#       nm / mm
#
# therefore:
#
#       delta_lambda = nm

delta_lambda_nm = (
    omega_prime_mm
    * reciprocal_dispersion
)


delta_lambda_pm = (
    delta_lambda_nm
    * 1000.0
)


##############################################################
# Spectral resolving power
##############################################################

resolving_power = (
    wavelengths_nm
    / delta_lambda_nm
)


##############################################################
# Build output table
##############################################################

results = pd.DataFrame(
    {
        "wavelength_nm":
            wavelengths_nm,

        "x_centroid_mm":
            x_centroid,

        "y_centroid_mm":
            y_centroid,

        "x_min_mm":
            x_minimum,

        "x_max_mm":
            x_maximum,

        "y_min_mm":
            y_minimum,

        "y_max_mm":
            y_maximum,

        "delta_x_um":
            delta_x_um,

        "delta_y_um":
            delta_y_um,

        "omega_prime_um":
            omega_prime_um,

        "reciprocal_dispersion_nm_per_mm":
            reciprocal_dispersion,

        "delta_lambda_nm":
            delta_lambda_nm,

        "delta_lambda_pm":
            delta_lambda_pm,

        "resolving_power":
            resolving_power,

        "valid_ray_count":
            valid_ray_count
    }
)


##############################################################
# Save results
##############################################################

output_file.parent.mkdir(
    parents=True,
    exist_ok=True
)


results.to_csv(
    output_file,
    index=False
)


##############################################################
# Summary: fiber image
##############################################################

print()

print(
    "Fiber image"
)

print(
    "-----------"
)


print(
    f"Delta X range:       "
    f"{np.nanmin(delta_x_um):.3f} - "
    f"{np.nanmax(delta_x_um):.3f} um"
)

print(
    f"Mean Delta X:        "
    f"{np.nanmean(delta_x_um):.3f} um"
)


print()


print(
    f"Delta Y range:       "
    f"{np.nanmin(delta_y_um):.3f} - "
    f"{np.nanmax(delta_y_um):.3f} um"
)

print(
    f"Mean Delta Y:        "
    f"{np.nanmean(delta_y_um):.3f} um"
)


##############################################################
# Summary: spectral resolution
##############################################################

print()

print(
    "Spectral resolution"
)

print(
    "-------------------"
)


print(
    f"omega' range:        "
    f"{np.nanmin(omega_prime_um):.3f} - "
    f"{np.nanmax(omega_prime_um):.3f} um"
)

print(
    f"Mean omega':         "
    f"{np.nanmean(omega_prime_um):.3f} um"
)


print()


print(
    f"Reciprocal dispersion range: "
    f"{np.nanmin(reciprocal_dispersion):.6f} - "
    f"{np.nanmax(reciprocal_dispersion):.6f} nm/mm"
)

print(
    f"Mean reciprocal dispersion:  "
    f"{np.nanmean(reciprocal_dispersion):.6f} nm/mm"
)


print()


print(
    f"Delta lambda range:  "
    f"{np.nanmin(delta_lambda_pm):.3f} - "
    f"{np.nanmax(delta_lambda_pm):.3f} pm"
)

print(
    f"Mean Delta lambda:   "
    f"{np.nanmean(delta_lambda_pm):.3f} pm"
)


print()


print(
    f"R range:             "
    f"{np.nanmin(resolving_power):.0f} - "
    f"{np.nanmax(resolving_power):.0f}"
)

print(
    f"Mean R:              "
    f"{np.nanmean(resolving_power):.0f}"
)


print()

print(
    "Results saved to:"
)

print(
    output_file
)


##############################################################
# Plot omega'
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    wavelengths_nm,
    omega_prime_um,
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mathrm{nm})$",
    fontsize=18
)

ax.set_ylabel(
    r"$\omega'\;(\mu\mathrm{m})$",
    fontsize=18
)


ax.minorticks_on()


ax.tick_params(
    which="major",
    direction="in",
    top=True,
    right=True,
    length=6,
    width=2.2,
    labelsize=16
)


ax.tick_params(
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.2
)


fig.tight_layout()

plt.show()


##############################################################
# Plot spectral purity
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    wavelengths_nm,
    delta_lambda_pm,
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mathrm{nm})$",
    fontsize=18
)

ax.set_ylabel(
    r"$\Delta\lambda\;(\mathrm{pm})$",
    fontsize=18
)


ax.minorticks_on()


ax.tick_params(
    which="major",
    direction="in",
    top=True,
    right=True,
    length=6,
    width=2.2,
    labelsize=16
)


ax.tick_params(
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.2
)


fig.tight_layout()

plt.show()


##############################################################
# Plot resolving power
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    wavelengths_nm,
    resolving_power,
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mathrm{nm})$",
    fontsize=18
)

ax.set_ylabel(
    r"$R=\lambda/\Delta\lambda$",
    fontsize=18
)


ax.minorticks_on()


ax.tick_params(
    which="major",
    direction="in",
    top=True,
    right=True,
    length=6,
    width=2.2,
    labelsize=16
)


ax.tick_params(
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.2
)


fig.tight_layout()

plt.show()


##############################################################
# End
##############################################################