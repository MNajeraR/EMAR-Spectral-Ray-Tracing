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
    x_centroid,
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mu\mathrm{m})$",
    fontsize=18
)

ax.set_ylabel(
    r"$X_{\mathrm{c}}\;(\mathrm{mm})$",
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


##############################################################
# Plot Y centroid versus wavelength
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    valid_wavelengths,
    y_centroid,
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mu\mathrm{m})$",
    fontsize=18
)

ax.set_ylabel(
    r"$Y_{\mathrm{c}}\;(\mathrm{mm})$",
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


##############################################################
# Plot spectral coordinate versus wavelength
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    valid_wavelengths,
    s,
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mu\mathrm{m})$",
    fontsize=18
)

ax.set_ylabel(
    r"$\mathrm{Spectral coordinate}\;(\mathrm{s})\;[\mathrm{mm}]$",
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


##############################################################
# Plot  linear dispersion
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    valid_wavelengths*1000.0,
    reciprocal_dispersion_nm,
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mathrm{nm})$",
    fontsize=18
)

ax.set_ylabel(
    r"$d\lambda/ds\;(\mathrm{nm}/\mathrm{mm})$",
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


##############################################################
# Spot projection onto the local spectral direction
##############################################################

# Compute the wavelength derivatives of the spectral centroid.
# Together, dx/dlambda and dy/dlambda define the local tangent
# direction of the spectral trace on the image plane.
dx_dlambda = np.gradient(
    x_centroid,
    valid_wavelengths
)

dy_dlambda = np.gradient(
    y_centroid,
    valid_wavelengths
)


# Normalize the centroid derivative vector by the local linear
# dispersion ds/dlambda. The resulting components (tx, ty)
# define a unit tangent vector along the local direction of
# spectral dispersion.
tx = (
    dx_dlambda
    / linear_dispersion
)

ty = (
    dy_dlambda
    / linear_dispersion
)


# Project each monochromatic spot onto its local spectral
# direction. For every wavelength, the ray coordinates are
# first expressed relative to the corresponding centroid and
# then projected onto the unit tangent vector.
#
# The resulting coordinate u measures the displacement of each
# ray along the local dispersion direction, in mm.
u = []

for i in range(len(valid_wavelengths)):

    wavelength = valid_wavelengths[i]

    spot = order_data[
        "spots"
    ][float(wavelength)]

    x = np.asarray(
        spot["x"],
        dtype=float
    )

    y = np.asarray(
        spot["y"],
        dtype=float
    )


    # Monochromatic spot centroid.
    xc = x_centroid[i]
    yc = y_centroid[i]


    # Ray positions relative to the spot centroid.
    delta_x = x - xc
    delta_y = y - yc


    # Projection onto the local spectral direction:
    #
    #     u = delta_r . t_hat
    #
    # where delta_r = (delta_x, delta_y) and
    # t_hat = (tx, ty).
    u_i = (
        delta_x * tx[i]
        + delta_y * ty[i]
    )

    u.append(
        u_i
    )


# Array shape:
#
#     (number of wavelengths, number of pupil rays)
#
# Each row therefore represents the projected monochromatic
# footprint for one wavelength.
u = np.asarray(
    u
)

##############################################################
# Diagnostic plot of one projected monochromatic spot
##############################################################

plot_index = (
    len(valid_wavelengths)
    // 2
)

wavelength = float(
    valid_wavelengths[plot_index]
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


# Move the monochromatic spot centroid to (0, 0).
delta_x = (
    x - x_centroid[plot_index]
)

delta_y = (
    y - y_centroid[plot_index]
)


fig, ax = plt.subplots(
    figsize=(7, 7)
)

ax.scatter(
    delta_x * 1000,
    delta_y * 1000,
    s=8,
    color="black"
)

ax.set_xlabel(
    r"$\Delta X\;(\mu\mathrm{m})$",
    fontsize=18
)

ax.set_ylabel(
    r"$\Delta Y\;(\mu\mathrm{m})$",
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

ax.set_aspect(
    "equal"
)

plt.tight_layout()
plt.show()

##############################################################
# Effective spectral image width
##############################################################

# In the absence of an explicit entrance slit in the Zemax
# model, the full aberrational extent of each monochromatic
# footprint along the local dispersion direction is adopted
# as an effective slit-image width:
#
#     omega' = u_max - u_min
#
# This corresponds to the maximum geometrical width of the
# monochromatic image in the spectral direction.
omega_prime = (
    np.max(u, axis=1)
    - np.min(u, axis=1)
)

print()

print(
    "Mean geometrical width:"
)
print(
    f"    {np.mean(omega_prime)*1000:.6f} um"
)


##############################################################
# Spectral purity
##############################################################

# Following the adopted spectral-purity criterion, Delta lambda
# is the wavelength interval for which the spectral displacement
# equals the effective image width omega':
#
#     delta l = omega'
#
# Since d(lambda)/ds is the reciprocal linear dispersion,
#
#     Delta lambda = omega' * d(lambda)/ds
#
# omega_prime is in mm and reciprocal_dispersion_nm is in nm/mm,
# therefore spectral_purity_nm is obtained in nm.
spectral_purity_nm = (
    omega_prime
    * reciprocal_dispersion_nm
)

print()

print(
    "Mean spectral purity:"
)
print(
    f"    {np.mean(spectral_purity_nm):.6f} nm"
)

##############################################################
# Plot spectral purity
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)


ax.plot(
    valid_wavelengths*1000,
    spectral_purity_nm*1000, 
    color="black",
    linewidth=2.0
)


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mathrm{nm})$",
    fontsize=18
)

ax.set_ylabel(
    r"$\delta\lambda\;(\mathrm{pm})$",
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

##############################################################
# Plot effective spectral image width
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.plot(
    valid_wavelengths * 1000.0,
    omega_prime * 1000.0, 
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

plt.tight_layout()
plt.show()

##############################################################
# Resolving power
##############################################################

wavelengths_nm = (
    valid_wavelengths
    * 1000.0
)

resolving_power = (
    wavelengths_nm
    / spectral_purity_nm
)

exponent = int(
    np.floor(
        np.log10(
            np.max(resolving_power)
        )
    )
)

scale = 10**exponent

print()

print(
    "Max Resolving power:"
)
print(
    f"    {np.max(resolving_power):.2f}"
)

##############################################################
# Plot resolving power
##############################################################

fig, ax = plt.subplots(
    figsize=(8, 5)
)

ax.plot(
    wavelengths_nm,
    resolving_power / scale,
    color="black",
    linewidth=2.0
)

ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\mathrm{nm})$",
    fontsize=18
)

ax.set_ylabel(
    rf"$R\;(10^{{{exponent}}})$",
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