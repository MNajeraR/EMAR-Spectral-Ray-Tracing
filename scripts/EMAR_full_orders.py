# -*- coding: utf-8 -*-

"""
EMAR_full_orders.py
===================

Dense ray tracing of the complete EMAR echelle spectral format.

The complete echelle-order definition is read directly from:

    Results/echelle_orders_60_144.csv

The table contains:

    order
    status
    wave_1 ... wave_11

for every diffraction order from m=60 to m=144.

For each order, a dense wavelength grid is generated between the
first and last reference wavelengths stored in the table. The order
is then traced through Zemax using the nearest lower original
configuration as the optical template.

The same normalized pupil sample is used for every wavelength and
every diffraction order.

Author
------
Morgan Rhaí Nájera Roa
"""

import sys
from pathlib import Path

import numpy as np
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
from utils import emar_plots


##############################################################
# Ray-tracing parameters
##############################################################

surf = 54

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
# Load Zemax file
##############################################################

ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Trace complete dense echelle format
##############################################################

results = {}


try:

    for i, row in enumerate(
        order_table,
        start=1
    ):

        ######################################################
        # Order information from reconstructed table
        ######################################################

        order = int(
            row["order"]
        )


        status = row[
            "status"
        ]


        ######################################################
        # Reference wavelengths from CSV
        ######################################################

        reference_wavelengths = np.asarray(
            [
                row[f"wave_{j}"]
                for j in range(1, 12)
            ],
            dtype=float
        )


        ######################################################
        # Dense wavelength sampling
        ######################################################

        wavelengths_dense = np.linspace(
            reference_wavelengths[0],
            reference_wavelengths[-1],
            n_wavelengths
        )


        ######################################################
        # Select nearest lower original Zemax configuration
        ######################################################

        base_config = (
            (order - 60) // 7
            + 1
        )


        ######################################################
        # Trace echelle order
        ######################################################

        results[
            order
        ] = emar_utils.trace_all_order(
            ln=ln,
            base_config=base_config,
            target_order=order,
            wavelengths=wavelengths_dense,
            surf=surf,
            px=px,
            py=py,
            row_order=row_order,
            row_wave1=row_wave1
        )


        ######################################################
        # Preserve status metadata
        ######################################################

        results[
            order
        ]["status"] = status


        ######################################################
        # Progress bar
        ######################################################

        percentage = (
            100.0
            * i
            / len(order_table)
        )


        bar_length = 30


        filled_length = int(
            bar_length
            * i
            / len(order_table)
        )


        bar = (
            "#"
            * filled_length
            + "-"
            * (
                bar_length
                - filled_length
            )
        )


        print(
            f"\r[{bar}] "
            f"{percentage:6.2f}% | "
            f"m = {order}",
            end="",
            flush=True
        )


    print()


finally:

    ln.close()


##############################################################
# Plot complete ray footprint
##############################################################

emar_plots.plot_echelle_orders(
    spots_by_order=results,
    surf=surf,
    view_mode="scatter",
    point_size=1,
    alpha=0.20,
    save=True,
    file_format="png"
)


##############################################################
# Plot spectral centroid traces
##############################################################

emar_plots.plot_echelle_orders(
    spots_by_order=results,
    surf=surf,
    view_mode="spectral_trace",
    linewidth=1.0,
    save=True,
    file_format="png"
)


##############################################################
# Plot detector representation
##############################################################

emar_plots.plot_echelle_orders(
    spots_by_order=results,
    surf=surf,
    view_mode="detector",
    linewidth=1.0,
    save=True,
    file_format="png"
)