# -*- coding: utf-8 -*-
"""
EMAR - Complete Echelle Order Reconstruction
============================================

This experiment validates the reconstruction of the complete echelle
spectral format from diffraction order m = 60 to m = 144.

The original Zemax model contains 13 representative configurations:

    m = 60, 67, 74, 81, 88, 95, 102,
        109, 116, 123, 130, 137, 144

The remaining intermediate orders are reconstructed using the wavelength
table stored in:

    Results/echelle_orders_60_144.csv

For this validation experiment, only the 11 reference wavelengths
associated with each diffraction order are traced.

Original orders are traced using their native Zemax configurations.
Intermediate orders are reconstructed by temporarily modifying the
echelle diffraction order and WAVE 1 MCE entries.

The resulting spectral centroids are traced to Surface 54, corresponding
to the final image plane.

The purpose of this experiment is to verify that the complete sequence
of 85 echelle orders produces a continuous and physically consistent
spectral format before generating a densely sampled version.

Author
------
Morgan R. Najera

Date
----
September 2026
"""


##############################################################
# Standard-library imports
##############################################################

from pathlib import Path
import sys
import csv


##############################################################
# EMAR project path
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
# External packages
##############################################################

import numpy as np
import matplotlib.pyplot as plt
import pyzdde.zdde as pyz


##############################################################
# EMAR utilities
##############################################################

from utils import emar_utils


##############################################################
# Analysis parameters
##############################################################

# Final image plane.

surf = 54


# Pupil sampling.

n_rays = 100
seed = 123


# Multi-Configuration Editor rows.

row_order = 1
row_wave1 = 2


##############################################################
# Original Zemax configurations and diffraction orders
##############################################################

# Mapping:
#
#     diffraction order -> Zemax configuration

original_configs = {
    60: 1,
    67: 2,
    74: 3,
    81: 4,
    88: 5,
    95: 6,
    102: 7,
    109: 8,
    116: 9,
    123: 10,
    130: 11,
    137: 12,
    144: 13
}


##############################################################
# Complete order sequence
##############################################################

orders = np.arange(
    60,
    145
)


##############################################################
# Generate common pupil sampling
##############################################################

px, py = emar_utils.sample_random_pupil(
    n_rays=n_rays,
    seed=seed
)


##############################################################
# File locations
##############################################################

zemax_file = (
    emar_utils.ZEMAX_DIR
    / "WP - Con prismas diseñados - camara.zmx"
)


orders_file = (
    emar_utils.RESULTS_DIR
    / "echelle_orders_60_144.csv"
)


##############################################################
# Read reconstructed wavelength table
##############################################################

wavelengths_by_order = {}




with open(
    orders_file,
    "r",
    encoding="utf-8"
) as f:

    reader = csv.DictReader(
        f
    )


    for row in reader:

        order = int(
            row["order"]
        )


        wavelengths = np.array(
            [
                float(
                    row[f"wave_{j}"]
                )
                for j in range(
                    1,
                    12
                )
            ]
        )


        wavelengths_by_order[
            order
        ] = wavelengths



##############################################################
# Connect to Zemax
##############################################################

ln = pyz.createLink()


if ln is None:

    raise RuntimeError(
        "Could not establish a PyZDDE connection to Zemax OpticStudio."
    )


ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Trace complete order sequence
##############################################################

results = {}


print(
    "\nTracing complete echelle order sequence..."
)

print(
    f"Orders: {orders[0]} - {orders[-1]}"
)

print(
    f"Total orders: {len(orders)}"
)

print(
    f"Wavelengths per order: 11"
)

print()


for i, order in enumerate(
    orders,
    start=1
):

    order = int(
        order
    )

    ##########################################################
    # Original Zemax order
    ##########################################################

    if order in original_configs:

        config = original_configs[
            order
        ]

        results[
            order
        ] = emar_utils.trace_original_order(
            ln=ln,
            config=config,
            order=order,
            surf=surf,
            px=px,
            py=py
        )

    ##########################################################
    # Reconstructed intermediate order
    ##########################################################

    else:

        lower_original_orders = [
            original_order
            for original_order in original_configs
            if original_order < order
        ]

        base_order = max(
            lower_original_orders
        )

        base_config = original_configs[
            base_order
        ]

        results[
            order
        ] = emar_utils.trace_reconstructed_order(
            ln=ln,
            base_config=base_config,
            target_order=order,
            wavelengths=wavelengths_by_order[
                order
            ],
            surf=surf,
            px=px,
            py=py,
            row_order=row_order,
            row_wave1=row_wave1
        )

    ##########################################################
    # Progress bar
    ##########################################################

    percentage = 100.0 * i / len(orders)

    bar_length = 30

    filled_length = int(
        bar_length * i / len(orders)
    )

    bar = (
        "#" * filled_length
        + "-" * (bar_length - filled_length)
    )

    print(
        f"\r[{bar}] "
        f"{percentage:6.2f}% | "
        f"m = {order}",
        end="",
        flush=True
    )

##############################################################
# Close Zemax
##############################################################

ln.close()


##############################################################
# Calculate spectral centroids
##############################################################

centroids = {}


for order in orders:

    order = int(
        order
    )


    result = results[
        order
    ]


    xc = []
    yc = []


    for wavelength in result[
        "wavelengths"
    ]:

        spot = result[
            "spots"
        ][
            float(wavelength)
        ]


        x = spot[
            "x"
        ]

        y = spot[
            "y"
        ]


        if len(x) == 0:

            xc.append(
                np.nan
            )

            yc.append(
                np.nan
            )


        else:

            xc.append(
                np.mean(x)
            )

            yc.append(
                np.mean(y)
            )


    centroids[
        order
    ] = {
        "xc": np.asarray(
            xc
        ),
        "yc": np.asarray(
            yc
        )
    }


##############################################################
# Plot complete spectral format
##############################################################

fig, ax = plt.subplots(
    figsize=(12, 10)
)


##############################################################
# Colormap
##############################################################

cmap = plt.get_cmap(
    "turbo_r"
)


norm = plt.Normalize(
    orders.min(),
    orders.max()
)


##############################################################
# Plot each diffraction order
##############################################################

for order in orders:

    order = int(
        order
    )


    ax.plot(
        centroids[order]["xc"],
        centroids[order]["yc"],
        "-",
        linewidth=1.1,
        color=cmap(
            norm(order)
        )
    )

##############################################################
# Plot formatting
##############################################################

ax.set_xlabel(
    "X centroid [mm]",
    fontsize=14
)


ax.set_ylabel(
    "Y centroid [mm]",
    fontsize=14
)


ax.set_title(
    "Complete reconstructed echelle format "
    "(m = 60-144) - Surface 54",
    fontsize=15
)


ax.tick_params(
    axis="both",
    labelsize=12
)


##############################################################
# Colorbar
##############################################################

sm = plt.cm.ScalarMappable(
    norm=norm,
    cmap=cmap
)


sm.set_array(
    []
)


cbar = fig.colorbar(
    sm,
    ax=ax
)


cbar.set_label(
    "Echelle diffraction order m",
    fontsize=12
)


##############################################################
# Final layout
##############################################################

fig.tight_layout()


##############################################################
# Save figure
##############################################################

output_file = (
    emar_utils.RESULTS_DIR
    / "Complete_echelle_orders_60_144_surf54.png"
)


fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)


plt.show()


print(
    f"\nFigure saved to:\n{output_file}"
)