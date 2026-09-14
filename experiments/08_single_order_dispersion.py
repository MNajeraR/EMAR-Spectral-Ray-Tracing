# -*- coding: utf-8 -*-

"""
08_single_order_dispersion.py
=============================

Single-order spectral-dispersion experiment for EMAR.

This experiment traces one echelle diffraction order at the final
image plane and evaluates how the spectral centroid changes with
wavelength.

The complete order definition is read from:

    Results/echelle_orders_60_144.csv

For the selected order:

1. Read the 11 reference wavelengths from the master order table.
2. Generate a dense wavelength sampling between wave_1 and wave_11.
3. Trace the same normalized pupil sample at every wavelength.
4. Compute the spectral centroid X(lambda), Y(lambda).
5. Construct the cumulative coordinate s along the spectral trace.
6. Estimate the local linear dispersion ds/dlambda.
7. Compute the reciprocal linear dispersion dlambda/ds.

This experiment is intended as the first step toward the quantitative
spectral characterization of EMAR.

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

surf = 54

target_order = 102

n_rays = 100

n_wavelengths = 100

seed = 123


##############################################################
# MCE rows
##############################################################

row_order = 1

row_wave1 = 2


##############################################################
# Load complete echelle-order table
##############################################################

orders_file = (
    emar_utils.RESULTS_DIR
    / "echelle_orders_60_144.csv"
)


order_table = np.genfromtxt(
    orders_file,
    delimiter=",",
    names=True,
    dtype=None,
    encoding="utf-8"
)


##############################################################
# Select target order
##############################################################

selected_rows = order_table[
    order_table["order"] == target_order
]


if len(selected_rows) == 0:

    raise ValueError(
        f"Order m={target_order} "
        "was not found in the order table."
    )


row = selected_rows[0]


status = row["status"]


reference_wavelengths = np.asarray(
    [
        row[f"wave_{j}"]
        for j in range(1, 12)
    ],
    dtype=float
)


##############################################################
# Dense wavelength sampling
##############################################################

wavelengths = np.linspace(
    reference_wavelengths[0],
    reference_wavelengths[-1],
    n_wavelengths
)


##############################################################
# Select Zemax template configuration
##############################################################

base_config = (
    (target_order - 60) // 7
    + 1
)


##############################################################
# Common normalized pupil sample
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
    / "WP - Con prismas diseñados - camara.zmx"
)


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
# Trace selected echelle order
##############################################################

try:

    order_data = emar_utils.trace_all_order(
        ln=ln,
        base_config=base_config,
        target_order=target_order,
        wavelengths=wavelengths,
        surf=surf,
        px=px,
        py=py,
        row_order=row_order,
        row_wave1=row_wave1
    )


finally:

    ln.close()


##############################################################
# Compute spectral centroids
##############################################################

x_centroid = []

y_centroid = []

valid_wavelengths = []


for wavelength in wavelengths:

    wavelength = float(
        wavelength
    )


    spot = order_data[
        "spots"
    ][wavelength]


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


    x_centroid.append(
        np.mean(x)
    )

    y_centroid.append(
        np.mean(y)
    )

    valid_wavelengths.append(
        wavelength
    )


x_centroid = np.asarray(
    x_centroid
)

y_centroid = np.asarray(
    y_centroid
)

valid_wavelengths = np.asarray(
    valid_wavelengths
)


##############################################################
# Construct coordinate along spectral trace
##############################################################

dx = np.diff(
    x_centroid
)

dy = np.diff(
    y_centroid
)


ds = np.sqrt(
    dx**2
    + dy**2
)


s = np.concatenate(
    (
        [0.0],
        np.cumsum(ds)
    )
)


##############################################################
# Local linear dispersion
##############################################################

linear_dispersion = np.gradient(
    s,
    valid_wavelengths
)


##############################################################
# Reciprocal linear dispersion
##############################################################

reciprocal_dispersion_um = (
    1.0
    / linear_dispersion
)


reciprocal_dispersion_nm = (
    1000.0
    * reciprocal_dispersion_um
)


##############################################################
# Print experiment summary
##############################################################

print()

print(
    "EMAR single-order dispersion experiment"
)

print(
    "--------------------------------------"
)

print(
    f"Order:             m = {target_order}"
)

print(
    f"Status:            {status}"
)

print(
    f"Base config:       {base_config}"
)

print(
    f"Image surface:     {surf}"
)

print(
    f"Pupil rays:        {n_rays}"
)

print(
    f"Wavelength points: {len(valid_wavelengths)}"
)

print(
    f"Wavelength range:  "
    f"{valid_wavelengths[0]:.6f} - "
    f"{valid_wavelengths[-1]:.6f} um"
)

print()

print(
    "Mean linear dispersion:"
)

print(
    f"    {np.mean(linear_dispersion):.6f} mm/um"
)

print()

print(
    "Mean reciprocal dispersion:"
)

print(
    f"    {np.mean(reciprocal_dispersion_nm):.6f} nm/mm"
)


##############################################################
# Plot X centroid versus wavelength
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    valid_wavelengths,
    x_centroid
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mu\mathrm{m})$",
    fontsize=16
)

ax.set_ylabel(
    "X centroid (mm)",
    fontsize=16
)

ax.tick_params(
    axis="both",
    labelsize=16
)

plt.tight_layout()
plt.show()


##############################################################
# Plot Y centroid versus wavelength
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    valid_wavelengths,
    y_centroid
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mu\mathrm{m})$",
    fontsize=16
)

ax.set_ylabel(
    "Y centroid (mm)",
    fontsize=16
)


ax.tick_params(
    axis="both",
    labelsize=16
)

plt.tight_layout()
plt.show()


##############################################################
# Plot spectral coordinate versus wavelength
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    valid_wavelengths,
    s
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mu\mathrm{m})$",
    fontsize=16
)

ax.set_ylabel(
    "Spectral coordinate (s) [mm]",
    fontsize=16
)


ax.tick_params(
    axis="both",
    labelsize=16
)

plt.tight_layout()
plt.show()


##############################################################
# Plot linear dispersion
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    valid_wavelengths,
    linear_dispersion
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mu\mathrm{m})$",
    fontsize=16
)

ax.set_ylabel(
    r"$ds/d\lambda\;(\mathrm{mm}/\mu\mathrm{m})$",
    fontsize=16
)

ax.tick_params(
    axis="both",
    labelsize=16
)

plt.tight_layout()
plt.show()

