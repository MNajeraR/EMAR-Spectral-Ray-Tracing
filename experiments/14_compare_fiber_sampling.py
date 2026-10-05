# -*- coding: utf-8 -*-

"""
14_compare_fiber_sampling.py
============================

Comparison between the robust 91-field fiber sampling and the
reduced five-field fiber sampling for the complete EMAR echelle
spectral format.

No Zemax ray tracing is performed.

The script reads:

    Results/full_orders_fiber_resolution.csv
    Results/full_orders_fiber_resolution_5fields.csv

and compares, wavelength by wavelength and order by order:

    omega' = Delta X
    resolving power R

The following differences are calculated:

    Delta omega' = omega'_5 - omega'_91

    Relative omega' difference =
        100 * (omega'_5 - omega'_91) / omega'_91

    Delta R = R_5 - R_91

    Relative R difference =
        100 * (R_5 - R_91) / R_91

Outputs:

    Results/fiber_sampling_comparison.csv
    Results/fiber_sampling_statistics_by_order.csv

    Results/omega_comparison_5_vs_91.png
    Results/omega_relative_difference_5_vs_91.png
    Results/resolving_power_comparison_5_vs_91.png
    Results/resolving_power_relative_difference_5_vs_91.png

Author
------
Morgan Rhaí Nájera Roa
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import numpy as np
import pandas as pd


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
# Input files
##############################################################

file_91 = (
    emar_utils.RESULTS_DIR
    / "full_orders_fiber_resolution.csv"
)


file_5 = (
    emar_utils.RESULTS_DIR
    / "full_orders_fiber_resolution_5fields.csv"
)


##############################################################
# Output files
##############################################################

output_comparison = (
    emar_utils.RESULTS_DIR
    / "fiber_sampling_comparison.csv"
)


output_statistics = (
    emar_utils.RESULTS_DIR
    / "fiber_sampling_statistics_by_order.csv"
)


output_omega = (
    emar_utils.RESULTS_DIR
    / "omega_comparison_5_vs_91.png"
)


output_omega_relative = (
    emar_utils.RESULTS_DIR
    / "omega_relative_difference_5_vs_91.png"
)


output_R = (
    emar_utils.RESULTS_DIR
    / "resolving_power_comparison_5_vs_91.png"
)


output_R_relative = (
    emar_utils.RESULTS_DIR
    / "resolving_power_relative_difference_5_vs_91.png"
)


##############################################################
# Load data
##############################################################

df_91 = pd.read_csv(
    file_91
)


df_5 = pd.read_csv(
    file_5
)


##############################################################
# Keep required quantities
##############################################################

df_91 = df_91[
    [
        "order",
        "wavelength_nm",
        "omega_prime_um",
        "resolving_power"
    ]
].copy()


df_5 = df_5[
    [
        "order",
        "wavelength_nm",
        "omega_prime_um",
        "resolving_power"
    ]
].copy()


##############################################################
# Rename quantities according to sampling
##############################################################

df_91 = df_91.rename(
    columns={
        "omega_prime_um":
            "omega_91_um",

        "resolving_power":
            "R_91"
    }
)


df_5 = df_5.rename(
    columns={
        "omega_prime_um":
            "omega_5_um",

        "resolving_power":
            "R_5"
    }
)


##############################################################
# Match order and wavelength
##############################################################

comparison = pd.merge(
    df_91,
    df_5,
    on=[
        "order",
        "wavelength_nm"
    ],
    how="inner",
    validate="one_to_one"
)


##############################################################
# Wavelength in Angstrom
##############################################################

comparison[
    "wavelength_A"
] = (
    comparison[
        "wavelength_nm"
    ]
    * 10.0
)


##############################################################
# Absolute omega difference
##############################################################

comparison[
    "delta_omega_um"
] = (
    comparison[
        "omega_5_um"
    ]
    - comparison[
        "omega_91_um"
    ]
)


##############################################################
# Relative omega difference
##############################################################

comparison[
    "relative_omega_percent"
] = (
    100.0
    * comparison[
        "delta_omega_um"
    ]
    / comparison[
        "omega_91_um"
    ]
)


##############################################################
# Absolute resolving-power difference
##############################################################

comparison[
    "delta_R"
] = (
    comparison[
        "R_5"
    ]
    - comparison[
        "R_91"
    ]
)


##############################################################
# Relative resolving-power difference
##############################################################

comparison[
    "relative_R_percent"
] = (
    100.0
    * comparison[
        "delta_R"
    ]
    / comparison[
        "R_91"
    ]
)


##############################################################
# Sort
##############################################################

comparison = comparison.sort_values(
    [
        "order",
        "wavelength_nm"
    ]
).reset_index(
    drop=True
)


##############################################################
# Save complete comparison
##############################################################

comparison.to_csv(
    output_comparison,
    index=False
)


##############################################################
# Statistics by diffraction order
##############################################################

statistics = (
    comparison
    .groupby(
        "order"
    )
    .agg(
        omega_difference_mean_um=(
            "delta_omega_um",
            "mean"
        ),

        omega_difference_min_um=(
            "delta_omega_um",
            "min"
        ),

        omega_difference_max_um=(
            "delta_omega_um",
            "max"
        ),

        omega_relative_mean_percent=(
            "relative_omega_percent",
            "mean"
        ),

        omega_relative_median_percent=(
            "relative_omega_percent",
            "median"
        ),

        omega_relative_std_percent=(
            "relative_omega_percent",
            "std"
        ),

        omega_relative_min_percent=(
            "relative_omega_percent",
            "min"
        ),

        omega_relative_max_percent=(
            "relative_omega_percent",
            "max"
        ),

        R_difference_mean=(
            "delta_R",
            "mean"
        ),

        R_difference_min=(
            "delta_R",
            "min"
        ),

        R_difference_max=(
            "delta_R",
            "max"
        ),

        R_relative_mean_percent=(
            "relative_R_percent",
            "mean"
        ),

        R_relative_median_percent=(
            "relative_R_percent",
            "median"
        ),

        R_relative_std_percent=(
            "relative_R_percent",
            "std"
        ),

        R_relative_min_percent=(
            "relative_R_percent",
            "min"
        ),

        R_relative_max_percent=(
            "relative_R_percent",
            "max"
        )
    )
    .reset_index()
)


##############################################################
# Save statistics
##############################################################

statistics.to_csv(
    output_statistics,
    index=False
)


##############################################################
# Global statistics
##############################################################

print()

print(
    "Fiber sampling comparison"
)

print(
    "-------------------------"
)

print(
    f"Matched points: "
    f"{len(comparison):,}"
)

print()


print(
    "omega' relative difference [%]"
)

print(
    f"Mean:    "
    f"{comparison['relative_omega_percent'].mean():+.4f}"
)

print(
    f"Median:  "
    f"{comparison['relative_omega_percent'].median():+.4f}"
)

print(
    f"Std:     "
    f"{comparison['relative_omega_percent'].std():.4f}"
)

print(
    f"Minimum: "
    f"{comparison['relative_omega_percent'].min():+.4f}"
)

print(
    f"Maximum: "
    f"{comparison['relative_omega_percent'].max():+.4f}"
)


print()


print(
    "Resolving-power relative difference [%]"
)

print(
    f"Mean:    "
    f"{comparison['relative_R_percent'].mean():+.4f}"
)

print(
    f"Median:  "
    f"{comparison['relative_R_percent'].median():+.4f}"
)

print(
    f"Std:     "
    f"{comparison['relative_R_percent'].std():.4f}"
)

print(
    f"Minimum: "
    f"{comparison['relative_R_percent'].min():+.4f}"
)

print(
    f"Maximum: "
    f"{comparison['relative_R_percent'].max():+.4f}"
)

print()


##############################################################
# Plot settings
##############################################################

orders = np.sort(
    comparison[
        "order"
    ].unique()
)


cmap = plt.cm.turbo_r


norm = Normalize(
    vmin=orders.min(),
    vmax=orders.max()
)


##############################################################
# 1. omega' comparison
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


for order in orders:

    subset = comparison[
        comparison[
            "order"
        ] == order
    ]


    color = cmap(
        norm(order)
    )


    ax.plot(
        subset[
            "wavelength_A"
        ],
        subset[
            "omega_91_um"
        ],
        color=color,
        linewidth=1.5,
        alpha=0.45
    )


    ax.plot(
        subset[
            "wavelength_A"
        ],
        subset[
            "omega_5_um"
        ],
        color=color,
        linewidth=1.0
    )


ax.set_xlabel(
    r"Wavelength [$\AA$]",
    fontsize=18
)


ax.set_ylabel(
    r"Projected fiber width $\omega'$ [$\mu$m]",
    fontsize=18
)


ax.set_title(
    r"Projected Fiber Width: 5 vs. 91 Fields",
    fontsize=18
)


ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5
)


ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0
)


ax.minorticks_on()


fig.tight_layout()


fig.savefig(
    output_omega,
    dpi=300,
    bbox_inches="tight"
)


##############################################################
# 2. Relative omega difference
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


for order in orders:

    subset = comparison[
        comparison[
            "order"
        ] == order
    ]


    ax.plot(
        subset[
            "wavelength_A"
        ],
        subset[
            "relative_omega_percent"
        ],
        linewidth=1.2,
        color=cmap(
            norm(order)
        )
    )


ax.axhline(
    0.0,
    color="black",
    linewidth=1.0,
    linestyle="--"
)


ax.set_xlabel(
    r"Wavelength [$\AA$]",
    fontsize=18
)


ax.set_ylabel(
    r"$100(\omega'_5-\omega'_{91})/\omega'_{91}$ [%]",
    fontsize=18
)


ax.set_title(
    r"Relative Difference in Projected Fiber Width",
    fontsize=18
)


ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5
)


ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0
)


ax.minorticks_on()


fig.tight_layout()


fig.savefig(
    output_omega_relative,
    dpi=300,
    bbox_inches="tight"
)


##############################################################
# 3. Resolving-power comparison
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


for order in orders:

    subset = comparison[
        comparison[
            "order"
        ] == order
    ]


    color = cmap(
        norm(order)
    )


    ax.plot(
        subset[
            "wavelength_A"
        ],
        subset[
            "R_91"
        ],
        color=color,
        linewidth=1.5,
        alpha=0.45
    )


    ax.plot(
        subset[
            "wavelength_A"
        ],
        subset[
            "R_5"
        ],
        color=color,
        linewidth=1.0
    )


ax.set_xlabel(
    r"Wavelength [$\AA$]",
    fontsize=18
)


ax.set_ylabel(
    "Resolving power R",
    fontsize=18
)


ax.set_title(
    "Spectral Resolving Power: 5 vs. 91 Fields",
    fontsize=18
)


ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5
)


ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0
)


ax.minorticks_on()


fig.tight_layout()


fig.savefig(
    output_R,
    dpi=300,
    bbox_inches="tight"
)


##############################################################
# 4. Relative resolving-power difference
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


for order in orders:

    subset = comparison[
        comparison[
            "order"
        ] == order
    ]


    ax.plot(
        subset[
            "wavelength_A"
        ],
        subset[
            "relative_R_percent"
        ],
        linewidth=1.2,
        color=cmap(
            norm(order)
        )
    )


ax.axhline(
    0.0,
    color="black",
    linewidth=1.0,
    linestyle="--"
)


ax.set_xlabel(
    r"Wavelength [$\AA$]",
    fontsize=18
)


ax.set_ylabel(
    r"$100(R_5-R_{91})/R_{91}$ [%]",
    fontsize=18
)


ax.set_title(
    "Relative Difference in Spectral Resolving Power",
    fontsize=18
)


ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5
)


ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0
)


ax.minorticks_on()


fig.tight_layout()


fig.savefig(
    output_R_relative,
    dpi=300,
    bbox_inches="tight"
)


##############################################################
# Show plots
##############################################################

plt.show()