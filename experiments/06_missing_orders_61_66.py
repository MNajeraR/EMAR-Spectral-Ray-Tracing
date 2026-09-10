# -*- coding: utf-8 -*-
"""
EMAR - Reconstructed Missing Echelle Orders 61-66
=================================================

This experiment extends the validation performed for the reconstructed
intermediate order m = 61.

The original Zemax model contains:

    Configuration 1 -> m = 60
    Configuration 2 -> m = 67

The reconstructed wavelength table stored in

    Results/echelle_orders_60_144.csv

contains the wavelength values corresponding to every integer order
between m = 60 and m = 144.

In this experiment, the missing orders

    m = 61, 62, 63, 64, 65, 66

are reconstructed using Configuration 1 as a temporary template.

The spectral traces of all reconstructed orders are compared with the
original m = 60 and m = 67 orders at the final image plane, Surface 54.

For each reconstructed order:

1. Configuration 1 is activated.
2. The original diffraction order and WAVE 1 MCE entries are stored.
3. PRAM(Surface 15, Parameter 2) is temporarily changed.
4. The 11 reconstructed wavelengths are read from the CSV.
5. Each wavelength is temporarily assigned to WAVE 1.
6. The same normalized pupil sampling is traced to Surface 54.
7. The original MCE entries are restored.

The purpose is to verify that the reconstructed intermediate orders
progressively fill the detector space between the original m = 60 and
m = 67 spectral traces.

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

# Original neighboring configurations.

config_60 = 1
config_67 = 2


# Original neighboring orders.

order_lower = 60
order_upper = 67


# Reconstructed missing orders.

missing_orders = np.arange(
    61,
    67
)


# All orders displayed in the experiment.

orders_to_plot = np.arange(
    60,
    68
)


# Final image plane.

surf = 54


##############################################################
# MCE rows
##############################################################

# Row 1:
#     PRAM Surface 15 / Parameter 2
#     Echelle diffraction order
#
# Row 2:
#     WAVE 1

row_order = 1
row_wave1 = 2


##############################################################
# Pupil sampling
##############################################################

n_rays = 100
seed = 123


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
# Helper: trace one original Zemax configuration
##############################################################

def trace_original_order(
    ln,
    config,
    order,
    surf,
    px,
    py
):
    """
    Trace the 11 original wavelength slots of a Zemax configuration.
    """

    ln.zSetConfig(
        config
    )

    ln.zGetUpdate()


    wavelengths = []
    spots = {}


    for wave_num in range(
        1,
        12
    ):

        wavelength = ln.zGetWave(
            wave_num
        ).wavelength


        x, y = emar_utils.trace_pupil(
            ln=ln,
            wave_num=wave_num,
            surf=surf,
            px=px,
            py=py,
            hx=0.0,
            hy=0.0
        )


        wavelengths.append(
            wavelength
        )


        spots[
            float(wavelength)
        ] = {
            "x": x.copy(),
            "y": y.copy()
        }


    return {
        "order": order,
        "wavelengths": np.asarray(
            wavelengths
        ),
        "spots": spots
    }


##############################################################
# Helper: trace one reconstructed intermediate order
##############################################################

def trace_reconstructed_order(
    ln,
    base_config,
    target_order,
    wavelengths,
    surf,
    px,
    py,
    row_order,
    row_wave1
):
    """
    Trace one reconstructed echelle order by temporarily modifying
    Configuration 1.

    The original diffraction order and WAVE 1 MCE entries are always
    restored before returning.
    """

    ##########################################################
    # Activate base configuration
    ##########################################################

    ln.zSetConfig(
        base_config
    )

    ln.zGetUpdate()


    ##########################################################
    # Save original MCE entries
    ##########################################################

    original_order = ln.zGetMulticon(
        base_config,
        row_order
    )


    original_wave1 = ln.zGetMulticon(
        base_config,
        row_wave1
    )


    ##########################################################
    # Storage
    ##########################################################

    spots = {}


    ##########################################################
    # Temporary modification
    ##########################################################

    try:

        ######################################################
        # Set reconstructed diffraction order
        ######################################################

        ln.zSetMulticon(
            base_config,
            row_order,
            float(target_order),
            original_order.status,
            original_order.pickupRow,
            original_order.pickupConfig,
            original_order.scale,
            original_order.offset
        )


        ln.zGetUpdate()


        ######################################################
        # Verify order
        ######################################################

        order_zemax = ln.zGetMulticon(
            base_config,
            row_order
        ).value


        print(
            f"\nOrder m = {target_order}"
        )

        print(
            f"Zemax PRAM = {order_zemax:.0f}"
        )


        ######################################################
        # Trace reconstructed wavelengths
        ######################################################

        for wavelength in wavelengths:

            ln.zSetMulticon(
                base_config,
                row_wave1,
                float(wavelength),
                original_wave1.status,
                original_wave1.pickupRow,
                original_wave1.pickupConfig,
                original_wave1.scale,
                original_wave1.offset
            )


            ln.zGetUpdate()


            wave_zemax = ln.zGetWave(
                1
            ).wavelength


            x, y = emar_utils.trace_pupil(
                ln=ln,
                wave_num=1,
                surf=surf,
                px=px,
                py=py,
                hx=0.0,
                hy=0.0
            )


            spots[
                float(wavelength)
            ] = {
                "x": x.copy(),
                "y": y.copy()
            }


            if len(x) != len(px):

                print(
                    f"  lambda = {wave_zemax:.9f} um | "
                    f"{len(x)}/{len(px)} valid rays"
                )


    ##########################################################
    # Restore original configuration
    ##########################################################

    finally:

        ln.zSetMulticon(
            base_config,
            row_wave1,
            original_wave1.value,
            original_wave1.status,
            original_wave1.pickupRow,
            original_wave1.pickupConfig,
            original_wave1.scale,
            original_wave1.offset
        )


        ln.zSetMulticon(
            base_config,
            row_order,
            original_order.value,
            original_order.status,
            original_order.pickupRow,
            original_order.pickupConfig,
            original_order.scale,
            original_order.offset
        )


        ln.zGetUpdate()


    return {
        "order": target_order,
        "wavelengths": wavelengths.copy(),
        "spots": spots
    }


##############################################################
# Trace original lower order: m = 60
##############################################################

print(
    "\nTracing original order m = 60..."
)


results = {}


results[
    order_lower
] = trace_original_order(
    ln=ln,
    config=config_60,
    order=order_lower,
    surf=surf,
    px=px,
    py=py
)


##############################################################
# Trace reconstructed orders: m = 61 ... 66
##############################################################

for order in missing_orders:

    print(
        f"\nTracing reconstructed order m = {order}..."
    )


    results[
        int(order)
    ] = trace_reconstructed_order(
        ln=ln,
        base_config=config_60,
        target_order=int(order),
        wavelengths=wavelengths_by_order[int(order)],
        surf=surf,
        px=px,
        py=py,
        row_order=row_order,
        row_wave1=row_wave1
    )


##############################################################
# Trace original upper order: m = 67
##############################################################

print(
    "\nTracing original order m = 67..."
)


results[
    order_upper
] = trace_original_order(
    ln=ln,
    config=config_67,
    order=order_upper,
    surf=surf,
    px=px,
    py=py
)


##############################################################
# Final restoration check
##############################################################

ln.zSetConfig(
    config_60
)

ln.zGetUpdate()


restored_order = ln.zGetMulticon(
    config_60,
    row_order
).value


restored_wave1 = ln.zGetMulticon(
    config_60,
    row_wave1
).value


print(
    "\nFinal restoration check:"
)

print(
    f"Order restored : {restored_order:.0f}"
)

print(
    f"WAVE 1 restored: {restored_wave1:.9f} um"
)


##############################################################
# Close Zemax
##############################################################

ln.close()


##############################################################
# Calculate spectral centroids
##############################################################

centroids = {}


for order in orders_to_plot:

    result = results[
        int(order)
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
        int(order)
    ] = {
        "xc": np.asarray(xc),
        "yc": np.asarray(yc)
    }


##############################################################
# Plot spectral traces
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 8)
)


cmap = plt.get_cmap(
    "viridis"
)


norm = plt.Normalize(
    order_lower,
    order_upper
)


for order in orders_to_plot:

    order = int(order)


    ax.plot(
        centroids[order]["xc"],
        centroids[order]["yc"],
        linewidth=2.5,
        color=cmap(norm(order)),
        label=f"m = {order}"
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
    "Reconstructed echelle orders 60-67 - Surface 54",
    fontsize=15
)


ax.legend(
    fontsize=9,
    ncol=2
)


ax.tick_params(
    axis="both",
    labelsize=12
)


fig.tight_layout()


##############################################################
# Save figure
##############################################################

output_file = (
    emar_utils.RESULTS_DIR
    / "Reconstructed_orders_60_67_surf54.png"
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

