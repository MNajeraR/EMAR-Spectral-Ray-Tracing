"""
Experiment 13
=============

Five-field circular-fiber spectral resolution analysis for EMAR.

This experiment repeats the full-order fiber analysis using only five
representative field points of the 100 um circular fiber:

    center, +X, -X, +Y, -Y

The pupil sampling remains unchanged at 91 hexapolar rays.

The purpose is to compare this reduced sampling against the robust
91-field fiber calculation while keeping the optical and spectral
analysis unchanged.

Orders:
    m = 60 ... 144

Fiber:
    100 um diameter
    5 representative fields

Pupil:
    5-ring hexapolar sampling
    91 rays

Spectral sampling:
    100 wavelengths per order

Detector:
    Surface 54
"""

##############################################################
# Imports
##############################################################

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyzdde.zdde as pyz
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


##############################################################
# Project paths
##############################################################

PROJECT_DIR = Path(__file__).resolve().parents[1]

UTILS_DIR = PROJECT_DIR / "utils"

if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))


import emar_utils


##############################################################
# Analysis parameters
##############################################################

surf = 54

row_order = 1
row_wave1 = 2

n_wavelengths = 100

fiber_diameter_um = 100.0

hx_max = 0.286480
hy_max = 0.285400


##############################################################
# Five representative fiber fields
##############################################################
#
#              +Y
#               |
#               3
#               |
#      -X  2 ---0--- 1  +X
#               |
#               4
#               |
#              -Y
#
##############################################################

fiber_hx = np.array(
    [
        0.0,
        +hx_max,
        -hx_max,
        0.0,
        0.0,
    ],
    dtype=float,
)

fiber_hy = np.array(
    [
        0.0,
        0.0,
        0.0,
        +hy_max,
        -hy_max,
    ],
    dtype=float,
)

n_fiber_points = len(fiber_hx)


##############################################################
# Pupil sampling
##############################################################

pupil_px, pupil_py = emar_utils.sample_hexapolar_pupil(
    n_rings=5,
    include_center=True,
)

n_pupil_rays = len(pupil_px)


##############################################################
# Input order table
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
    encoding="utf-8",
)


##############################################################
# Zemax model
##############################################################

zemax_file = (
    emar_utils.ZEMAX_DIR
    / "WP - Con prismas diseñados - camara - theoretical slit.zmx"
)


##############################################################
# Output files
##############################################################

output_csv = (
    emar_utils.RESULTS_DIR
    / "full_orders_fiber_resolution_5fields.csv"
)

output_dispersion = (
    emar_utils.RESULTS_DIR
    / "linear_dispersion_5fields.png"
)

output_resolution = (
    emar_utils.RESULTS_DIR
    / "resolving_power_5fields.png"
)

output_format = (
    emar_utils.RESULTS_DIR
    / "spectral_format_5fields.png"
)

output_envelope = (
    emar_utils.RESULTS_DIR
    / "spectral_format_fiber_envelope_5fields.png"
)


##############################################################
# Print experiment information
##############################################################

print()
print("=" * 70)
print("EMAR FIVE-FIELD FIBER ANALYSIS")
print("=" * 70)

print(
    f"Fiber diameter      : {fiber_diameter_um:.1f} um"
)

print(
    f"Fiber fields        : {n_fiber_points}"
)

print(
    f"Pupil rays          : {n_pupil_rays}"
)

print(
    f"Rays / wavelength   : "
    f"{n_fiber_points * n_pupil_rays}"
)

print(
    f"Wavelengths / order : {n_wavelengths}"
)

print(
    f"Orders              : "
    f"{len(order_table)}"
)

print(
    f"Total ray traces    : "
    f"{n_fiber_points * n_pupil_rays * n_wavelengths * len(order_table):,}"
)

print("=" * 70)
print()


##############################################################
# Connect to Zemax
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
# Storage
##############################################################

all_rows = []

total_calls = (
    len(order_table)
    * n_fiber_points
)

completed_calls = 0

start_time = time.time()


##############################################################
# Trace all diffraction orders
##############################################################

try:

    for order_index, row in enumerate(
        order_table,
        start=1,
    ):

        order = int(
            row["order"]
        )

        status = row["status"]


        ######################################################
        # Dense wavelength sampling
        ######################################################

        reference_wavelengths = np.asarray(
            [
                row[f"wave_{j}"]
                for j in range(1, 12)
            ],
            dtype=float,
        )

        wavelengths = np.linspace(
            reference_wavelengths[0],
            reference_wavelengths[-1],
            n_wavelengths,
        )


        ######################################################
        # Zemax template configuration
        ######################################################

        base_config = (
            (order - 60) // 7
        ) + 1


        ######################################################
        # Storage for this order
        ######################################################

        x_by_wavelength = [
            []
            for _ in wavelengths
        ]

        y_by_wavelength = [
            []
            for _ in wavelengths
        ]


        ######################################################
        # Trace five fiber fields
        ######################################################

        for field_index, (hx, hy) in enumerate(
            zip(
                fiber_hx,
                fiber_hy,
            ),
            start=1,
        ):

            result = emar_utils.trace_all_order(
                ln=ln,
                base_config=base_config,
                target_order=order,
                wavelengths=wavelengths,
                surf=surf,
                px=pupil_px,
                py=pupil_py,
                row_order=row_order,
                row_wave1=row_wave1,
                hx=float(hx),
                hy=float(hy),
            )


            ##################################################
            # Collect footprints wavelength by wavelength
            ##################################################

            spots = result[
                "spots"
            ]

            if len(spots) != len(wavelengths):

                raise RuntimeError(
                    f"Unexpected number of wavelengths "
                    f"for order m={order}: "
                    f"{len(spots)} != {len(wavelengths)}"
                )


            for wavelength_index, wavelength in enumerate(
                wavelengths
            ):

                spot = spots[
                    float(wavelength)
                ]

                x = np.asarray(
                    spot["x"],
                    dtype=float,
                )

                y = np.asarray(
                    spot["y"],
                    dtype=float,
                )

                if len(x) > 0:

                    x_by_wavelength[
                        wavelength_index
                    ].append(
                        x
                    )

                    y_by_wavelength[
                        wavelength_index
                    ].append(
                        y
                    )


            ##################################################
            # Global progress
            ##################################################

            completed_calls += 1

            fraction = (
                completed_calls
                / total_calls
            )

            bar_length = 30

            filled = int(
                fraction
                * bar_length
            )

            bar = (
                "=" * filled
                + "-"
                * (
                    bar_length
                    - filled
                )
            )

            print(
                "\r"
                f"[{bar}] "
                f"{100.0 * fraction:6.2f}% "
                f"| m={order:3d} "
                f"| fiber "
                f"{field_index}/{n_fiber_points}",
                end="",
                flush=True,
            )


        ######################################################
        # Analyze detector fiber image
        ######################################################

        x_centroid = np.full(
            len(wavelengths),
            np.nan,
        )

        y_centroid = np.full(
            len(wavelengths),
            np.nan,
        )

        x_min = np.full(
            len(wavelengths),
            np.nan,
        )

        x_max = np.full(
            len(wavelengths),
            np.nan,
        )

        y_min = np.full(
            len(wavelengths),
            np.nan,
        )

        y_max = np.full(
            len(wavelengths),
            np.nan,
        )

        delta_x_um = np.full(
            len(wavelengths),
            np.nan,
        )

        delta_y_um = np.full(
            len(wavelengths),
            np.nan,
        )

        n_valid_rays = np.zeros(
            len(wavelengths),
            dtype=int,
        )


        for i in range(
            len(wavelengths)
        ):

            if not x_by_wavelength[i]:
                continue

            x_all = np.concatenate(
                x_by_wavelength[i]
            )

            y_all = np.concatenate(
                y_by_wavelength[i]
            )


            x_centroid[i] = np.mean(
                x_all
            )

            y_centroid[i] = np.mean(
                y_all
            )


            x_min[i] = np.min(
                x_all
            )

            x_max[i] = np.max(
                x_all
            )

            y_min[i] = np.min(
                y_all
            )

            y_max[i] = np.max(
                y_all
            )


            delta_x_um[i] = (
                x_max[i]
                - x_min[i]
            ) * 1000.0

            delta_y_um[i] = (
                y_max[i]
                - y_min[i]
            ) * 1000.0


            n_valid_rays[i] = len(
                x_all
            )


        ######################################################
        # Wavelength units
        ######################################################

        wavelengths_nm = (
            wavelengths
            * 1000.0
        )


        ######################################################
        # Reciprocal linear dispersion
        #
        # d(lambda) / dX
        #
        # wavelength : nm
        # X          : mm
        #
        # result     : nm/mm
        ######################################################

        dx_dlambda = np.gradient(
            x_centroid,
            wavelengths_nm,
        )

        reciprocal_dispersion = (
            1.0
            / np.abs(
                dx_dlambda
            )
        )


        ######################################################
        # Convert to Angstrom / mm
        ######################################################

        linear_dispersion_A_per_mm = (
            reciprocal_dispersion
            * 10.0
        )


        ######################################################
        # Projected fiber width
        ######################################################

        omega_prime_um = (
            delta_x_um.copy()
        )

        omega_prime_mm = (
            omega_prime_um
            / 1000.0
        )


        ######################################################
        # Spectral purity
        #
        # delta_lambda =
        # omega' * d(lambda)/dX
        ######################################################

        delta_lambda_nm = (
            omega_prime_mm
            * reciprocal_dispersion
        )


        ######################################################
        # Resolving power
        ######################################################

        resolving_power = (
            wavelengths_nm
            / delta_lambda_nm
        )


        ######################################################
        # Store rows
        ######################################################

        for i, wavelength in enumerate(
            wavelengths
        ):

            all_rows.append(
                {
                    "order": order,
                    "status": status,
                    "wavelength_um": wavelength,
                    "wavelength_nm": wavelengths_nm[i],
                    "x_centroid_mm": x_centroid[i],
                    "y_centroid_mm": y_centroid[i],
                    "x_min_mm": x_min[i],
                    "x_max_mm": x_max[i],
                    "y_min_mm": y_min[i],
                    "y_max_mm": y_max[i],
                    "delta_x_um": delta_x_um[i],
                    "delta_y_um": delta_y_um[i],
                    "omega_prime_um": omega_prime_um[i],
                    "reciprocal_dispersion_nm_per_mm":
                        reciprocal_dispersion[i],
                    "linear_dispersion_A_per_mm":
                        linear_dispersion_A_per_mm[i],
                    "delta_lambda_nm":
                        delta_lambda_nm[i],
                    "resolving_power":
                        resolving_power[i],
                    "n_valid_rays":
                        n_valid_rays[i],
                }
            )


        ######################################################
        # Save after every completed order
        ######################################################

        df_partial = pd.DataFrame(
            all_rows
        )

        df_partial.to_csv(
            output_csv,
            index=False,
        )


        print(
            f"\nCompleted m={order} "
            f"({order_index}/{len(order_table)})"
        )


##############################################################
# Always close Zemax
##############################################################

finally:

    ln.close()


##############################################################
# Final dataframe
##############################################################

df = pd.DataFrame(
    all_rows
)


##############################################################
# Runtime
##############################################################

elapsed = (
    time.time()
    - start_time
)

print()
print("=" * 70)

print(
    f"Finished in "
    f"{elapsed / 60.0:.2f} min"
)

print(
    f"Saved: {output_csv}"
)

print("=" * 70)
print()


##############################################################
# Plot settings
##############################################################

orders = np.sort(
    df["order"].unique()
)

cmap = plt.cm.turbo_r

norm = Normalize(
    vmin=orders.min(),
    vmax=orders.max(),
)


##############################################################
# 1. Linear dispersion
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)

for order in orders:

    subset = df[
        df["order"] == order
    ]

    ax.plot(
        subset["wavelength_nm"] * 10.0,
        subset["linear_dispersion_A_per_mm"],
        linewidth=1.2,
        color=cmap(
            norm(order)
        ),
    )


ax.set_xlabel(
    r"Wavelength [$\AA$]",
    fontsize=18,
)

ax.set_ylabel(
    r"Linear dispersion [$\AA$/mm]",
    fontsize=18,
)

ax.set_title(
    "EMAR Linear Dispersion",
    fontsize=18,
)

ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5,
)

ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0,
)

ax.minorticks_on()


sm = ScalarMappable(
    norm=norm,
    cmap=cmap,
)

sm.set_array([])

cbar = fig.colorbar(
    sm,
    ax=ax,
)

cbar.set_label(
    "Echelle diffraction order m",
    fontsize=16,
)

cbar.ax.tick_params(
    labelsize=14
)


fig.tight_layout()

fig.savefig(
    output_dispersion,
    dpi=300,
    bbox_inches="tight",
)


##############################################################
# 2. Resolving power
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)

for order in orders:

    subset = df[
        df["order"] == order
    ]

    ax.plot(
        subset["wavelength_nm"] * 10.0,
        subset["resolving_power"],
        linewidth=1.2,
        color=cmap(
            norm(order)
        ),
    )


ax.set_xlabel(
    r"Wavelength [$\AA$]",
    fontsize=18,
)

ax.set_ylabel(
    "Resolving power R",
    fontsize=18,
)

ax.set_title(
    "EMAR Spectral Resolving Power",
    fontsize=18,
)

ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5,
)

ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0,
)

ax.minorticks_on()


sm = ScalarMappable(
    norm=norm,
    cmap=cmap,
)

sm.set_array([])

cbar = fig.colorbar(
    sm,
    ax=ax,
)

cbar.set_label(
    "Echelle diffraction order m",
    fontsize=16,
)

cbar.ax.tick_params(
    labelsize=14
)


fig.tight_layout()

fig.savefig(
    output_resolution,
    dpi=300,
    bbox_inches="tight",
)


##############################################################
# 3. Echelle spectral format - centroids
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 7)
)

for order in orders:

    subset = df[
        df["order"] == order
    ]

    ax.plot(
        subset["x_centroid_mm"],
        subset["y_centroid_mm"],
        linewidth=1.2,
        color=cmap(
            norm(order)
        ),
    )


ax.set_xlabel(
    "X centroid [mm]",
    fontsize=18,
)

ax.set_ylabel(
    "Y centroid [mm]",
    fontsize=18,
)

ax.set_title(
    "EMAR Echelle Spectral Format",
    fontsize=18,
)

ax.set_aspect(
    "equal",
    adjustable="box",
)

ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5,
)

ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0,
)

ax.minorticks_on()


sm = ScalarMappable(
    norm=norm,
    cmap=cmap,
)

sm.set_array([])

cbar = fig.colorbar(
    sm,
    ax=ax,
)

cbar.set_label(
    "Echelle diffraction order m",
    fontsize=16,
)

cbar.ax.tick_params(
    labelsize=14
)


fig.tight_layout()

fig.savefig(
    output_format,
    dpi=300,
    bbox_inches="tight",
)


##############################################################
# 4. Echelle spectral format with fiber envelope
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 7)
)

for order in orders:

    subset = df[
        df["order"] == order
    ]

    color = cmap(
        norm(order)
    )

    ax.fill_between(
        subset["x_centroid_mm"],
        subset["y_min_mm"],
        subset["y_max_mm"],
        color=color,
        alpha=0.30,
        linewidth=0.0,
    )

    ax.plot(
        subset["x_centroid_mm"],
        subset["y_centroid_mm"],
        color=color,
        linewidth=1.0,
    )


ax.set_xlabel(
    "X centroid [mm]",
    fontsize=18,
)

ax.set_ylabel(
    "Y [mm]",
    fontsize=18,
)

ax.set_title(
    "EMAR Echelle Format with 100 um Fiber Envelope",
    fontsize=18,
)

ax.set_aspect(
    "equal",
    adjustable="box",
)

ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5,
)

ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0,
)

ax.minorticks_on()


sm = ScalarMappable(
    norm=norm,
    cmap=cmap,
)

sm.set_array([])

cbar = fig.colorbar(
    sm,
    ax=ax,
)

cbar.set_label(
    "Echelle diffraction order m",
    fontsize=16,
)

cbar.ax.tick_params(
    labelsize=14
)


fig.tight_layout()

fig.savefig(
    output_envelope,
    dpi=300,
    bbox_inches="tight",
)


##############################################################
# Show plots
##############################################################

plt.show()