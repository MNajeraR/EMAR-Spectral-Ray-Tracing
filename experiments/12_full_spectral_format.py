# -*- coding: utf-8 -*-

"""
12_full_spectral_format.py
==========================

Full EMAR echelle spectral analysis for orders m = 60 ... 144.

For every diffraction order:

1. Trace the central fiber position over the wavelength range.
2. Calculate the spectral trace and reciprocal linear dispersion.
3. Trace the complete 100 um circular fiber ONLY at the two
   wavelength extremes of the order.
4. Measure:

       omega'_blue = Delta X(lambda_min)
       omega'_red  = Delta X(lambda_max)

5. Calculate resolving power only at those two extreme points:

       delta_lambda = omega' * d(lambda)/dx
       R = lambda / delta_lambda

6. Plot:
       - Linear dispersion for all orders
       - Resolving power using only the two extreme points
       - Full echelle spectral format with fiber-image width

Color convention:

       low order  -> red
       high order -> blue

This follows the physical wavelength ordering of the echelle.

Author
------
Morgan Rhaí Nájera Roa
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
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
    prefix="",
    bar_length=40
):

    fraction = current / total

    filled = int(
        bar_length * fraction
    )

    bar = (
        "=" * filled
        + "-" * (bar_length - filled)
    )

    print(
        f"\r{prefix}"
        f"[{bar}] "
        f"{100.0*fraction:6.2f}% "
        f"({current}/{total})",
        end="",
        flush=True
    )

    if current == total:
        print()


##############################################################
# Orders
##############################################################

orders = np.arange(
    60,
    145
)

n_orders = len(
    orders
)


##############################################################
# Optical surfaces
##############################################################

surf_detector = 54


##############################################################
# Fiber
##############################################################

fiber_diameter_um = 100.0

n_fiber_rings = 5
n_pupil_rings = 5


##############################################################
# Calibrated fiber field coordinates
##############################################################

hx_max = 0.286480
hy_max = 0.285400


##############################################################
# MCE rows
##############################################################

row_order = 1
row_wave1 = 2


##############################################################
# Number of wavelengths used to describe each spectral order
##############################################################

n_wavelengths = 100


##############################################################
# Zemax model
##############################################################

zemax_file = (
    PROJECT_DIR
    / "Zemax"
    / "WP - Con prismas diseñados - camara - theoretical slit.zmx"
)


##############################################################
# Output files
##############################################################

results_dir = (
    PROJECT_DIR
    / "Results"
)

results_dir.mkdir(
    parents=True,
    exist_ok=True
)


spectral_file = (
    results_dir
    / "full_spectral_format.csv"
)


resolution_file = (
    results_dir
    / "full_fiber_resolution_extremes.csv"
)


##############################################################
# Fiber sampling
##############################################################

fiber_x_norm, fiber_y_norm = (
    emar_utils.sample_hexapolar_pupil(
        n_rings=n_fiber_rings,
        include_center=True
    )
)


hx_fiber = (
    hx_max
    * fiber_x_norm
)

hy_fiber = (
    hy_max
    * fiber_y_norm
)


##############################################################
# Pupil sampling
##############################################################

px, py = (
    emar_utils.sample_hexapolar_pupil(
        n_rings=n_pupil_rings,
        include_center=True
    )
)


##############################################################
# Central ray for spectral trace
##############################################################

px_center = np.array([
    0.0
])

py_center = np.array([
    0.0
])


##############################################################
# Storage
##############################################################

spectral_rows = []

resolution_rows = []


##############################################################
# Helper: trace complete fiber image
##############################################################

def trace_fiber_image(
    ln,
    wavelength_um,
    order,
    base_config
):

    all_x = []
    all_y = []


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
            target_order=int(order),
            base_config=int(base_config),
            row_order=row_order,
            row_wave1=row_wave1,
            hx=float(hx),
            hy=float(hy)
        )


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


        if len(x) > 0:

            all_x.extend(
                x
            )

            all_y.extend(
                y
            )


    all_x = np.asarray(
        all_x,
        dtype=float
    )

    all_y = np.asarray(
        all_y,
        dtype=float
    )


    if len(all_x) == 0:

        return {
            "xc": np.nan,
            "yc": np.nan,
            "xmin": np.nan,
            "xmax": np.nan,
            "ymin": np.nan,
            "ymax": np.nan,
            "dx_um": np.nan,
            "dy_um": np.nan,
            "n_rays": 0
        }


    xmin = np.min(all_x)
    xmax = np.max(all_x)

    ymin = np.min(all_y)
    ymax = np.max(all_y)


    return {

        "xc":
            np.mean(all_x),

        "yc":
            np.mean(all_y),

        "xmin":
            xmin,

        "xmax":
            xmax,

        "ymin":
            ymin,

        "ymax":
            ymax,

        "dx_um":
            (xmax - xmin) * 1000.0,

        "dy_um":
            (ymax - ymin) * 1000.0,

        "n_rays":
            len(all_x)
    }


##############################################################
# Open Zemax
##############################################################

ln = pyz.createLink()

ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Main calculation
##############################################################

try:

    print()

    print(
        "EMAR full spectral-format analysis"
    )

    print(
        "----------------------------------"
    )

    print(
        f"Orders:            "
        f"{orders.min()} - {orders.max()}"
    )

    print(
        f"Number of orders:  {n_orders}"
    )

    print(
        f"Fiber points:      {len(hx_fiber)}"
    )

    print(
        f"Pupil rays:        {len(px)}"
    )

    print()


    ##########################################################
    # Loop over diffraction orders
    ##########################################################

    for order_index, order in enumerate(
        orders
    ):

        ######################################################
        # Configuration
        ######################################################

        base_config = (
            (int(order) - 60) // 7 + 1
        )


        ######################################################
        # IMPORTANT
        #
        # Obtain wavelength range for this order.
        #
        # The existing EMAR utilities should provide the
        # wavelength range associated with each diffraction
        # order/configuration.
        ######################################################

        wavelengths_um = (
            emar_utils.get_order_wavelengths(
                ln=ln,
                target_order=int(order),
                base_config=int(base_config),
                n_wavelengths=n_wavelengths,
                row_order=row_order,
                row_wave1=row_wave1
            )
        )


        wavelengths_um = np.asarray(
            wavelengths_um,
            dtype=float
        )


        wavelengths_nm = (
            wavelengths_um
            * 1000.0
        )


        ######################################################
        # Trace spectral centerline
        ######################################################

        trace = emar_utils.trace_all_order(
            ln=ln,
            wavelengths=wavelengths_um,
            px=px_center,
            py=py_center,
            surf=surf_detector,
            target_order=int(order),
            base_config=int(base_config),
            row_order=row_order,
            row_wave1=row_wave1,
            hx=0.0,
            hy=0.0
        )


        ######################################################
        # Extract centroid coordinates
        ######################################################

        x_centroid = []
        y_centroid = []


        spots = list(
            trace["spots"].values()
        )


        if len(spots) != len(wavelengths_um):

            raise RuntimeError(
                f"Order {order}: "
                f"{len(spots)} spectral points returned "
                f"for {len(wavelengths_um)} wavelengths."
            )


        for spot in spots:

            x = np.asarray(
                spot["x"],
                dtype=float
            )

            y = np.asarray(
                spot["y"],
                dtype=float
            )


            if len(x) == 0:

                x_centroid.append(
                    np.nan
                )

                y_centroid.append(
                    np.nan
                )

            else:

                x_centroid.append(
                    np.mean(x)
                )

                y_centroid.append(
                    np.mean(y)
                )


        x_centroid = np.asarray(
            x_centroid,
            dtype=float
        )

        y_centroid = np.asarray(
            y_centroid,
            dtype=float
        )


        ######################################################
        # Reciprocal dispersion along X
        ######################################################

        # dx/dlambda:
        #
        #       mm / nm

        dx_dlambda = np.gradient(
            x_centroid,
            wavelengths_nm
        )


        # d(lambda)/dx:
        #
        #       nm / mm

        reciprocal_dispersion = (
            1.0
            / np.abs(dx_dlambda)
        )


        ######################################################
        # Convert to Angstrom/mm
        ######################################################

        linear_dispersion_A_per_mm = (
            reciprocal_dispersion
            * 10.0
        )


        ######################################################
        # Save spectral trace
        ######################################################

        for i in range(
            len(wavelengths_nm)
        ):

            spectral_rows.append(
                {
                    "order":
                        int(order),

                    "wavelength_nm":
                        wavelengths_nm[i],

                    "wavelength_A":
                        wavelengths_nm[i] * 10.0,

                    "x_centroid_mm":
                        x_centroid[i],

                    "y_centroid_mm":
                        y_centroid[i],

                    "reciprocal_dispersion_nm_per_mm":
                        reciprocal_dispersion[i],

                    "linear_dispersion_A_per_mm":
                        linear_dispersion_A_per_mm[i]
                }
            )


        ######################################################
        # Fiber image at BLUE edge
        ######################################################

        blue_index = 0

        blue_fiber = trace_fiber_image(
            ln=ln,
            wavelength_um=wavelengths_um[
                blue_index
            ],
            order=order,
            base_config=base_config
        )


        ######################################################
        # Fiber image at RED edge
        ######################################################

        red_index = (
            len(wavelengths_um)
            - 1
        )

        red_fiber = trace_fiber_image(
            ln=ln,
            wavelength_um=wavelengths_um[
                red_index
            ],
            order=order,
            base_config=base_config
        )


        ######################################################
        # Calculate R at the two extremes
        ######################################################

        for label, index, fiber in [

            (
                "blue",
                blue_index,
                blue_fiber
            ),

            (
                "red",
                red_index,
                red_fiber
            )
        ]:

            omega_prime_um = (
                fiber["dx_um"]
            )


            omega_prime_mm = (
                omega_prime_um
                / 1000.0
            )


            delta_lambda_nm = (
                omega_prime_mm
                * reciprocal_dispersion[
                    index
                ]
            )


            resolving_power = (
                wavelengths_nm[
                    index
                ]
                / delta_lambda_nm
            )


            resolution_rows.append(
                {
                    "order":
                        int(order),

                    "edge":
                        label,

                    "wavelength_nm":
                        wavelengths_nm[index],

                    "wavelength_A":
                        wavelengths_nm[index] * 10.0,

                    "x_centroid_mm":
                        fiber["xc"],

                    "y_centroid_mm":
                        fiber["yc"],

                    "x_min_mm":
                        fiber["xmin"],

                    "x_max_mm":
                        fiber["xmax"],

                    "y_min_mm":
                        fiber["ymin"],

                    "y_max_mm":
                        fiber["ymax"],

                    "delta_x_um":
                        fiber["dx_um"],

                    "delta_y_um":
                        fiber["dy_um"],

                    "omega_prime_um":
                        omega_prime_um,

                    "reciprocal_dispersion_nm_per_mm":
                        reciprocal_dispersion[index],

                    "delta_lambda_nm":
                        delta_lambda_nm,

                    "resolving_power":
                        resolving_power,

                    "valid_ray_count":
                        fiber["n_rays"]
                }
            )


        ######################################################
        # Progress
        ######################################################

        print_progress(
            current=order_index + 1,
            total=n_orders,
            prefix=f"m={order:3d}  "
        )


finally:

    pyz.closeLink()


##############################################################
# Convert results to DataFrames
##############################################################

spectral_df = pd.DataFrame(
    spectral_rows
)

resolution_df = pd.DataFrame(
    resolution_rows
)


##############################################################
# Save CSV files
##############################################################

spectral_df.to_csv(
    spectral_file,
    index=False
)


resolution_df.to_csv(
    resolution_file,
    index=False
)


##############################################################
# Color map
#
# Low order  -> red
# High order -> blue
##############################################################

cmap = plt.cm.turbo_r

norm = Normalize(
    vmin=orders.min(),
    vmax=orders.max()
)


##############################################################
# Plot 1
#
# Linear dispersion for all orders
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


for order in orders:

    data = spectral_df[
        spectral_df["order"] == order
    ]

    color = cmap(
        norm(order)
    )


    ax.plot(
        data["wavelength_A"],
        data["linear_dispersion_A_per_mm"],
        color=color,
        linewidth=1.4
    )


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\AA)$",
    fontsize=18
)

ax.set_ylabel(
    r"$\mathrm{Linear\ dispersion}\;(\AA/\mathrm{mm})$",
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


sm = plt.cm.ScalarMappable(
    norm=norm,
    cmap=cmap
)

sm.set_array([])


cbar = fig.colorbar(
    sm,
    ax=ax,
    pad=0.02
)

cbar.set_label(
    r"Diffraction order, $m$",
    fontsize=18
)

cbar.ax.tick_params(
    labelsize=14
)


fig.tight_layout()

plt.show()


##############################################################
# Plot 2
#
# Resolving power
#
# ONLY the two extreme fiber measurements of each order
# are used.
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


for order in orders:

    data = resolution_df[
        resolution_df["order"] == order
    ].sort_values(
        "wavelength_A"
    )


    color = cmap(
        norm(order)
    )


    ax.plot(
        data["wavelength_A"],
        data["resolving_power"],
        color=color,
        linewidth=1.5
    )


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\AA)$",
    fontsize=18
)

ax.set_ylabel(
    r"$\mathrm{Resolving\ power},\ R=\lambda/\Delta\lambda$",
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


sm = plt.cm.ScalarMappable(
    norm=norm,
    cmap=cmap
)

sm.set_array([])


cbar = fig.colorbar(
    sm,
    ax=ax,
    pad=0.02
)

cbar.set_label(
    r"Diffraction order, $m$",
    fontsize=18
)

cbar.ax.tick_params(
    labelsize=14
)


fig.tight_layout()

plt.show()


##############################################################
# Plot 3
#
# Echelle spectral format
# including projected fiber-image size
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 7)
)


for order in orders:

    spectral_order = spectral_df[
        spectral_df["order"] == order
    ].sort_values(
        "wavelength_nm"
    )


    fiber_order = resolution_df[
        resolution_df["order"] == order
    ].sort_values(
        "wavelength_nm"
    )


    color = cmap(
        norm(order)
    )


    ##########################################################
    # Spectral centerline
    ##########################################################

    x = spectral_order[
        "x_centroid_mm"
    ].to_numpy()

    y = spectral_order[
        "y_centroid_mm"
    ].to_numpy()


    ax.plot(
        x,
        y,
        color=color,
        linewidth=1.3
    )


    ##########################################################
    # Fiber width at both order extremes
    #
    # Draw the actual measured rectangular envelope:
    #
    #       Delta X x Delta Y
    #
    # around each extreme.
    ##########################################################

    for _, row in fiber_order.iterrows():

        xmin = row[
            "x_min_mm"
        ]

        xmax = row[
            "x_max_mm"
        ]

        ymin = row[
            "y_min_mm"
        ]

        ymax = row[
            "y_max_mm"
        ]


        ######################################################
        # Horizontal extent
        ######################################################

        ax.plot(
            [
                xmin,
                xmax
            ],
            [
                row["y_centroid_mm"],
                row["y_centroid_mm"]
            ],
            color=color,
            linewidth=2.0
        )


        ######################################################
        # Vertical extent
        ######################################################

        ax.plot(
            [
                row["x_centroid_mm"],
                row["x_centroid_mm"]
            ],
            [
                ymin,
                ymax
            ],
            color=color,
            linewidth=2.0
        )


ax.set_xlabel(
    r"$X\;(\mathrm{mm})$",
    fontsize=18
)

ax.set_ylabel(
    r"$Y\;(\mathrm{mm})$",
    fontsize=18
)


ax.set_aspect(
    "equal",
    adjustable="box"
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


sm = plt.cm.ScalarMappable(
    norm=norm,
    cmap=cmap
)

sm.set_array([])


cbar = fig.colorbar(
    sm,
    ax=ax,
    pad=0.02
)

cbar.set_label(
    r"Diffraction order, $m$",
    fontsize=18
)

cbar.ax.tick_params(
    labelsize=14
)


fig.tight_layout()

plt.show()


##############################################################
# Summary
##############################################################

print()

print(
    "Full EMAR analysis completed"
)

print(
    "----------------------------"
)


print(
    f"Orders:              "
    f"{orders.min()} - {orders.max()}"
)

print(
    f"Spectral points:     "
    f"{len(spectral_df)}"
)

print(
    f"Resolution points:   "
    f"{len(resolution_df)}"
)


print()


print(
    f"R minimum:           "
    f"{resolution_df['resolving_power'].min():.0f}"
)

print(
    f"R maximum:           "
    f"{resolution_df['resolving_power'].max():.0f}"
)

print(
    f"R mean:              "
    f"{resolution_df['resolving_power'].mean():.0f}"
)


print()


print(
    "Saved:"
)

print(
    spectral_file
)

print(
    resolution_file
)
