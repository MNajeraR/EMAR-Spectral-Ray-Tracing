# -*- coding: utf-8 -*-
"""
EMAR - First Reconstructed Intermediate Echelle Order
=====================================================

This experiment tests whether an intermediate echelle diffraction order
that is not explicitly present in the original Zemax Multi-Configuration
Editor (MCE) can be reproduced by temporarily modifying an existing
configuration.

The original Zemax model contains:

    Configuration 1 -> m = 60
    Configuration 2 -> m = 67

The reconstructed wavelength table stored in

    Results/echelle_orders_60_144.csv

contains the wavelength values corresponding to every integer order
between m = 60 and m = 144.

In this experiment, order m = 61 is generated using Configuration 1
as a temporary template.

The test compares three spectral orders:

    m = 60  -> original Zemax configuration
    m = 61  -> reconstructed intermediate order
    m = 67  -> original Zemax configuration

For the reconstructed order m = 61:

1. Configuration 1 is activated.
2. The original diffraction order and WAVE 1 MCE entries are stored.
3. PRAM(Surface 15, Parameter 2) is temporarily changed from 60 to 61.
4. The 11 reconstructed wavelengths for m = 61 are read from the CSV.
5. Each wavelength is temporarily assigned to WAVE 1.
6. The model is updated with zGetUpdate().
7. The same normalized pupil sampling is traced to surface 54.
8. The resulting ray footprints are stored.
9. The original diffraction order and WAVE 1 values are restored.

The purpose of the experiment is to verify that the reconstructed order
appears geometrically between the neighboring original echelle orders
before extending the method to the complete m = 60 ... 144 sequence.

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

# Original configurations used as reference.

config_60 = 1
config_67 = 2


# Corresponding original diffraction orders.

order_60 = 60
order_test = 61
order_67 = 67


# Final image plane.

surf = 54


# Multi-Configuration Editor rows.
#
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

# For this first validation experiment, use the same 100 random
# pupil coordinates for all three orders.

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

# Store wavelengths using diffraction order as dictionary key:
#
#     wavelengths_by_order[60]
#     wavelengths_by_order[61]
#     ...
#     wavelengths_by_order[144]

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
# Extract wavelengths required for this experiment
##############################################################

waves_60 = wavelengths_by_order[
    order_60
]

waves_61 = wavelengths_by_order[
    order_test
]

waves_67 = wavelengths_by_order[
    order_67
]


print(
    "\nWavelength ranges:"
)

print(
    f"m = {order_60}: "
    f"{waves_60[0]:.6f} - {waves_60[-1]:.6f} um"
)

print(
    f"m = {order_test}: "
    f"{waves_61[0]:.6f} - {waves_61[-1]:.6f} um"
)

print(
    f"m = {order_67}: "
    f"{waves_67[0]:.6f} - {waves_67[-1]:.6f} um"
)


##############################################################
# Connect to Zemax
##############################################################

ln = pyz.createLink()


ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Helper: trace original configuration wavelengths
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
    Trace the 11 original wavelength slots of one Zemax configuration.

    No MCE values are modified by this function.
    """

    ln.zSetConfig(
        config
    )

    ln.zGetUpdate()


    wavelengths = []

    spots_by_wave = {}


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


        spots_by_wave[
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
        "spots": spots_by_wave
    }


##############################################################
# Trace original order m = 60
##############################################################

print(
    "\nTracing original order m = 60..."
)


result_60 = trace_original_order(
    ln=ln,
    config=config_60,
    order=order_60,
    surf=surf,
    px=px,
    py=py
)


##############################################################
# Trace reconstructed order m = 61
##############################################################

print(
    "\nTracing reconstructed order m = 61..."
)


##############################################################
# Activate Configuration 1
##############################################################

# Configuration 1 is used as the temporary template because its
# original diffraction order, m = 60, is the nearest lower order
# to the reconstructed m = 61 order.

ln.zSetConfig(
    config_60
)

ln.zGetUpdate()


##############################################################
# Save original MCE entries
##############################################################

# Store the complete MCE information rather than only the numerical
# values so status, pickup, scale, and offset can also be restored.

original_order = ln.zGetMulticon(
    config_60,
    row_order
)


original_wave1 = ln.zGetMulticon(
    config_60,
    row_wave1
)


##############################################################
# Storage for reconstructed order
##############################################################

spots_61 = {}


##############################################################
# Temporarily modify Configuration 1
##############################################################

try:

    ##########################################################
    # Change diffraction order: 60 -> 61
    ##########################################################

    ln.zSetMulticon(
        config_60,
        row_order,
        float(order_test),
        original_order.status,
        original_order.pickupRow,
        original_order.pickupConfig,
        original_order.scale,
        original_order.offset
    )


    ln.zGetUpdate()


    ##########################################################
    # Verify diffraction order
    ##########################################################

    order_zemax = ln.zGetMulticon(
        config_60,
        row_order
    ).value


    print(
        f"Requested diffraction order : {order_test}"
    )

    print(
        f"Zemax diffraction order     : "
        f"{order_zemax:.0f}"
    )


    ##########################################################
    # Loop over reconstructed wavelengths
    ##########################################################

    for j, wavelength in enumerate(
        waves_61,
        start=1
    ):

        ######################################################
        # Replace WAVE 1
        ######################################################

        ln.zSetMulticon(
            config_60,
            row_wave1,
            float(wavelength),
            original_wave1.status,
            original_wave1.pickupRow,
            original_wave1.pickupConfig,
            original_wave1.scale,
            original_wave1.offset
        )


        ######################################################
        # Update Zemax
        ######################################################

        ln.zGetUpdate()


        ######################################################
        # Verify wavelength
        ######################################################

        wave_zemax = ln.zGetWave(
            1
        ).wavelength


        ######################################################
        # Trace pupil
        ######################################################

        x, y = emar_utils.trace_pupil(
            ln=ln,
            wave_num=1,
            surf=surf,
            px=px,
            py=py,
            hx=0.0,
            hy=0.0
        )


        ######################################################
        # Store footprint
        ######################################################

        spots_61[
            float(wavelength)
        ] = {
            "x": x.copy(),
            "y": y.copy()
        }


        ######################################################
        # Diagnostic output
        ######################################################

        print(
            f"WAVE {j:2d} | "
            f"Requested = {wavelength:.9f} um | "
            f"Zemax = {wave_zemax:.9f} um | "
            f"{len(x)}/{len(px)} valid rays"
        )


##############################################################
# Always restore original Zemax configuration
##############################################################

finally:

    ##########################################################
    # Restore diffraction order
    ##########################################################

    ln.zSetMulticon(
        config_60,
        row_order,
        original_order.value,
        original_order.status,
        original_order.pickupRow,
        original_order.pickupConfig,
        original_order.scale,
        original_order.offset
    )


    ##########################################################
    # Restore WAVE 1
    ##########################################################

    ln.zSetMulticon(
        config_60,
        row_wave1,
        original_wave1.value,
        original_wave1.status,
        original_wave1.pickupRow,
        original_wave1.pickupConfig,
        original_wave1.scale,
        original_wave1.offset
    )


    ##########################################################
    # Return Zemax to original state
    ##########################################################

    ln.zGetUpdate()


##############################################################
# Save reconstructed order
##############################################################

result_61 = {
    "order": order_test,
    "wavelengths": waves_61.copy(),
    "spots": spots_61
}


##############################################################
# Trace original order m = 67
##############################################################

print(
    "\nTracing original order m = 67..."
)


result_67 = trace_original_order(
    ln=ln,
    config=config_67,
    order=order_67,
    surf=surf,
    px=px,
    py=py
)


##############################################################
# Verify restoration
##############################################################

# Return to Configuration 1 and check that the original diffraction
# order and WAVE 1 value were recovered.

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
    "\nRestoration check:"
)

print(
    f"Order restored : "
    f"{restored_order:.0f}"
)

print(
    f"WAVE 1 restored: "
    f"{restored_wave1:.9f} um"
)


##############################################################
# Close Zemax
##############################################################

ln.close()


##############################################################
# Group results
##############################################################

results = {
    order_60: result_60,
    order_test: result_61,
    order_67: result_67
}


##############################################################
# Calculate centroid traces
##############################################################

centroids = {}


for order, result in results.items():

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
        "xc": np.asarray(xc),
        "yc": np.asarray(yc)
    }


##############################################################
# Plot spectral traces
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 8)
)


colors = {
    60: "tab:red",
    61: "black",
    67: "tab:orange"
}


for order in [
    60,
    61,
    67
]:

    ax.plot(
        centroids[order]["xc"],
        centroids[order]["yc"],
        linewidth=2.0,
        color=colors[order],
        label=f"Order m = {order}"
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
    f"Intermediate order validation - Surface {surf}",
    fontsize=15
)


ax.legend(
    fontsize=11
)


ax.tick_params(
    axis="both",
    labelsize=12
)


fig.tight_layout()

plt.show()

