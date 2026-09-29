# -*- coding: utf-8 -*-

"""
10_theoretical_slit.py
======================

Spectral resolving power of EMAR using a theoretical finite slit.

Order:
    m = 102

Theoretical slit:
    50 x 15 um

The slit is represented through normalized field coordinates:

    Hx_max = 0.14324
    Hy_max = 0.04281

These limits were previously calibrated and validated at surface 6.

For every wavelength:

1. Sample the theoretical slit.
2. Trace the same pupil distribution from every slit point.
3. Propagate all rays to the detector.
4. Determine the local spectral-dispersion direction.
5. Project the complete finite-slit image onto that direction.
6. Calculate the projected slit-image width omega'.
7. Calculate spectral purity:

       delta_lambda = omega' * (d lambda / ds)

8. Calculate resolving power:

       R = lambda / delta_lambda

The wavelength grid and reciprocal linear dispersion are taken from
the results of experiment 08.

Author
------
Morgan Rhaí Nájera Roa
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pyzdde.zdde as pyz


##############################################################
# Project path
##############################################################

PROJECT_DIR = Path(
    __file__
).resolve().parent.parent


if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(
        0,
        str(PROJECT_DIR)
    )


##############################################################
# EMAR utilities
##############################################################

from utils import emar_utils


##############################################################
# Experiment parameters
##############################################################

target_order = 102

surf_detector = 54

n_rays = 100
seed = 123


##############################################################
# Theoretical slit
##############################################################

slit_width_x = 0.050   # mm = 50 um
slit_width_y = 0.015   # mm = 15 um

hx_max = 0.14324
hy_max = 0.04281


##############################################################
# Slit sampling
##############################################################

n_slit_x = 11
n_slit_y = 5


hx_values = np.linspace(
    -hx_max,
    +hx_max,
    n_slit_x
)


hy_values = np.linspace(
    -hy_max,
    +hy_max,
    n_slit_y
)


field_points = [
    (hx, hy)
    for hy in hy_values
    for hx in hx_values
]


##############################################################
# MCE rows
##############################################################

row_order = 1
row_wave1 = 2


##############################################################
# Read dispersion results from experiment 08
##############################################################

dispersion_file = (
    emar_utils.RESULTS_DIR
    / "order_102_dispersion.csv"
)


dispersion_data = np.genfromtxt(
    dispersion_file,
    delimiter=",",
    names=True,
    dtype=float
)


##############################################################
# Extract spectral quantities
##############################################################

wavelengths_nm = np.asarray(
    dispersion_data["wavelength_nm"],
    dtype=float
)


x_centroid = np.asarray(
    dispersion_data["x_centroid_mm"],
    dtype=float
)


y_centroid = np.asarray(
    dispersion_data["y_centroid_mm"],
    dtype=float
)


reciprocal_dispersion = np.asarray(
    dispersion_data["reciprocal_dispersion_nm_per_mm"],
    dtype=float
)


##############################################################
# Convert wavelengths to micrometers for Zemax
##############################################################

wavelengths_um = (
    wavelengths_nm
    / 1000.0
)


##############################################################
# Local dispersion direction
#
# Important:
# dx/dlambda and dy/dlambda use the same wavelength units,
# therefore their normalization removes the wavelength unit.
##############################################################

dx_dlambda = np.gradient(
    x_centroid,
    wavelengths_nm
)


dy_dlambda = np.gradient(
    y_centroid,
    wavelengths_nm
)


tangent_norm = np.sqrt(
    dx_dlambda**2
    + dy_dlambda**2
)


tx = (
    dx_dlambda
    / tangent_norm
)


ty = (
    dy_dlambda
    / tangent_norm
)


##############################################################
# Select Zemax template configuration
##############################################################

base_config = (
    (target_order - 60) // 7
    + 1
)


##############################################################
# Random pupil sampling
#
# The same normalized pupil is used for every slit point
# and every wavelength.
##############################################################

px, py = emar_utils.sample_random_pupil(
    n_rays=n_rays,
    seed=seed
)


##############################################################
# Zemax model
##############################################################

zemax_file = (
    emar_utils.ZEMAX_DIR
    / "WP - Con prismas diseñados - camara - theoretical slit.zmx"
)


##############################################################
# Output arrays
##############################################################

omega_slit_um = np.full(
    len(wavelengths_nm),
    np.nan
)


delta_lambda_nm = np.full(
    len(wavelengths_nm),
    np.nan
)


resolving_power = np.full(
    len(wavelengths_nm),
    np.nan
)


valid_ray_count = np.zeros(
    len(wavelengths_nm),
    dtype=int
)


##############################################################
# Print experiment setup
##############################################################

print()

print(
    "EMAR theoretical-slit spectral-resolution experiment"
)

print(
    "----------------------------------------------------"
)

print(
    f"Order:              m = {target_order}"
)

print(
    f"Base config:        {base_config}"
)

print(
    f"Detector surface:   {surf_detector}"
)

print(
    f"Slit size:          "
    f"{slit_width_x * 1000.0:.1f} x "
    f"{slit_width_y * 1000.0:.1f} um"
)

print(
    f"Hx limits:          +/- {hx_max:.6f}"
)

print(
    f"Hy limits:          +/- {hy_max:.6f}"
)

print(
    f"Slit sampling:      "
    f"{n_slit_x} x {n_slit_y}"
)

print(
    f"Slit points:        {len(field_points)}"
)

print(
    f"Pupil rays/point:   {n_rays}"
)

print(
    f"Wavelengths:        {len(wavelengths_nm)}"
)

print(
    f"Maximum rays:       "
    f"{len(wavelengths_nm) * len(field_points) * n_rays}"
)

print()


##############################################################
# Open PyZDDE connection
##############################################################

ln = pyz.createLink()


if ln is None:
    raise RuntimeError(
        "Could not establish a PyZDDE connection "
        "to Zemax OpticStudio."
    )


##############################################################
# Load Zemax model
##############################################################

ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Trace finite slit for every wavelength
##############################################################

try:

    for i, wavelength_um in enumerate(
        wavelengths_um
    ):

        print(
            f"\rWavelength "
            f"{i + 1:03d}/"
            f"{len(wavelengths_um):03d}  "
            f"({wavelengths_nm[i]:.3f} nm)",
            end="",
            flush=True
        )


        ######################################################
        # Complete finite-slit image at this wavelength
        ######################################################

        wavelength_x = []
        wavelength_y = []


        for hx, hy in field_points:

            result = emar_utils.trace_all_order(
                ln=ln,
                base_config=base_config,
                target_order=target_order,
                wavelengths=np.array(
                    [wavelength_um]
                ),
                surf=surf_detector,
                px=px,
                py=py,
                row_order=row_order,
                row_wave1=row_wave1,
                hx=float(hx),
                hy=float(hy)
            )


            spot = result[
                "spots"
            ][float(wavelength_um)]


            x = np.asarray(
                spot["x"],
                dtype=float
            )


            y = np.asarray(
                spot["y"],
                dtype=float
            )


            if len(x) == 0:
                continue


            wavelength_x.extend(
                x
            )


            wavelength_y.extend(
                y
            )


        ######################################################
        # Convert complete slit image to arrays
        ######################################################

        wavelength_x = np.asarray(
            wavelength_x,
            dtype=float
        )


        wavelength_y = np.asarray(
            wavelength_y,
            dtype=float
        )


        valid_ray_count[i] = len(
            wavelength_x
        )


        if len(wavelength_x) == 0:
            continue


        ######################################################
        # Finite-slit centroid
        ######################################################

        xc = np.mean(
            wavelength_x
        )


        yc = np.mean(
            wavelength_y
        )


        ######################################################
        # Project every detector ray onto the local
        # dispersion direction
        ######################################################

        u = (
            (wavelength_x - xc) * tx[i]
            +
            (wavelength_y - yc) * ty[i]
        )


        ######################################################
        # Projected finite-slit image width
        #
        # u is in mm.
        ######################################################

        omega_mm = (
            np.max(u)
            - np.min(u)
        )


        omega_slit_um[i] = (
            omega_mm
            * 1000.0
        )


        ######################################################
        # Spectral purity
        #
        # reciprocal_dispersion:
        #
        #     nm / mm
        #
        # omega_mm:
        #
        #     mm
        #
        # therefore:
        #
        #     delta_lambda -> nm
        ######################################################

        delta_lambda_nm[i] = (
            omega_mm
            * reciprocal_dispersion[i]
        )


        ######################################################
        # Resolving power
        ######################################################

        if delta_lambda_nm[i] > 0.0:

            resolving_power[i] = (
                wavelengths_nm[i]
                / delta_lambda_nm[i]
            )


finally:

    ln.close()


print()
print()


##############################################################
# Valid results
##############################################################

valid = (
    np.isfinite(omega_slit_um)
    &
    np.isfinite(delta_lambda_nm)
    &
    np.isfinite(resolving_power)
)


if not np.any(valid):
    raise RuntimeError(
        "No valid resolving-power results were obtained."
    )


##############################################################
# Print summary
##############################################################

print(
    "Finite-slit results"
)

print(
    "----------------------------------------------------"
)

print(
    f"Omega' range:       "
    f"{np.nanmin(omega_slit_um):.3f} - "
    f"{np.nanmax(omega_slit_um):.3f} um"
)

print(
    f"Mean omega':        "
    f"{np.nanmean(omega_slit_um):.3f} um"
)

print(
    f"Delta lambda range: "
    f"{np.nanmin(delta_lambda_nm) * 1000.0:.3f} - "
    f"{np.nanmax(delta_lambda_nm) * 1000.0:.3f} pm"
)

print(
    f"Mean delta lambda:  "
    f"{np.nanmean(delta_lambda_nm) * 1000.0:.3f} pm"
)

print(
    f"R range:            "
    f"{np.nanmin(resolving_power):.0f} - "
    f"{np.nanmax(resolving_power):.0f}"
)

print(
    f"Mean R:             "
    f"{np.nanmean(resolving_power):.0f}"
)

print()


##############################################################
# Save results
##############################################################

output_file = (
    emar_utils.RESULTS_DIR
    / "order_102_theoretical_slit_resolution.csv"
)


output_data = np.column_stack(
    (
        wavelengths_nm,
        omega_slit_um,
        reciprocal_dispersion,
        delta_lambda_nm,
        resolving_power,
        valid_ray_count
    )
)


header = (
    "wavelength_nm,"
    "omega_slit_um,"
    "reciprocal_dispersion_nm_per_mm,"
    "delta_lambda_nm,"
    "resolving_power,"
    "valid_ray_count"
)


np.savetxt(
    output_file,
    output_data,
    delimiter=",",
    header=header,
    comments="",
    fmt="%.10e"
)


print(
    f"Results saved to:"
)

print(
    output_file
)

print()


##############################################################
# Resolving-power plot
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    wavelengths_nm[valid],
    resolving_power[valid] / 1.0e4,
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    "Wavelength (nm)",
    fontsize=18
)


ax.set_ylabel(
    r"$R$ ($10^5$)",
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


plt.tight_layout()

plt.show()