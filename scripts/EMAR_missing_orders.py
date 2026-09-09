# -*- coding: utf-8 -*-
"""
EMAR - Reconstruction of Missing Echelle Orders
================================================

This script reconstructs the wavelength sampling for every integer
echelle diffraction order between m = 60 and m = 144.

The original EMAR Zemax model contains 13 representative
multi-configurations associated with the following diffraction orders:

    60, 67, 74, 81, 88, 95, 102,
    109, 116, 123, 130, 137, 144

Each configuration contains 11 reference wavelengths defined in the
Zemax Multi-Configuration Editor (MCE).

The original configurations show that, for a fixed wavelength index j,
the product

    m * lambda_j

is approximately invariant across the 13 configurations.

For each of the 11 wavelength positions, this script therefore computes

    C_j = <m * lambda_j>

where the average is evaluated over all original configurations.

The standard deviation

    sigma_j = std(m * lambda_j)

is also calculated to verify the numerical consistency of this
relationship in the original Zemax model.

Once the reference constants C_j are obtained, the corresponding
wavelengths for every integer diffraction order are reconstructed from

    lambda_j(m) = C_j / m

for

    m = 60, 61, 62, ..., 144.

The resulting wavelength table contains both:

    ORIGINAL orders
        Orders already present in the Zemax multi-configuration model.

    NEW orders
        Intermediate diffraction orders reconstructed from the
        m*lambda relation.

The complete table is saved as

    Results/echelle_orders_60_144.csv

for use by subsequent scripts that will construct and trace the
continuous echelle-order sequence.

IMPORTANT
---------
This script only reads information from the original Zemax
Multi-Configuration Editor.

It does NOT add, delete, or modify Zemax configurations.

The generated wavelength table is an external reconstruction that will
be used later to build the missing configurations explicitly.

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

# This script is located inside:
#
#     EMAR/scripts/
#
# Shared project modules are located inside:
#
#     EMAR/utils/
#
# Therefore, the EMAR project root must be added to Python's
# module-search path before importing the project utilities.

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
import pyzdde.zdde as pyz


##############################################################
# EMAR utilities
##############################################################

from utils import emar_utils


##############################################################
# Zemax optical model
##############################################################

# ZEMAX_DIR is defined centrally in utils/emar_utils.py.
#
# Using the centralized project directory avoids dependence on
# Spyder's current working directory.

zemax_file = (
    emar_utils.ZEMAX_DIR
    / "WP - Con prismas diseñados - camara.zmx"
)


##############################################################
# Original Zemax configurations
##############################################################

# The original optical model contains 13 representative
# multi-configurations.

configs_original = np.arange(
    1,
    14
)


##############################################################
# Multi-Configuration Editor rows
##############################################################

# Row containing the echelle diffraction order.

row_order = 1


# Rows containing the 11 reference wavelengths:
#
#     MCE row  2 -> WAVE 1
#     ...
#     MCE row 12 -> WAVE 11

row_wave_first = 2
row_wave_last = 12


# Number of wavelength reference positions contained in the MCE.

n_waves = (
    row_wave_last
    - row_wave_first
    + 1
)


##############################################################
# Complete echelle-order sequence
##############################################################

# Generate every integer diffraction order from 60 through 144,
# including both limits.
#
# This produces 85 orders in total.

orders_complete = np.arange(
    60,
    145
)


##############################################################
# Connect to Zemax
##############################################################

# Create communication link between Python and Zemax OpticStudio.

ln = pyz.createLink()


# Load the original EMAR optical model.

ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Allocate arrays for original MCE data
##############################################################

n_configs = len(
    configs_original
)


# One diffraction order per original configuration.

orders_original = np.empty(
    n_configs,
    dtype=int
)


# Wavelength matrix:
#
#                 WAVE 1  ...  WAVE 11
#
# Config 1
# Config 2
#   ...
# Config 13
#
# Shape:
#
#     (13 configurations, 11 wavelengths)

waves_original = np.empty(
    (
        n_configs,
        n_waves
    ),
    dtype=float
)


##############################################################
# Read original MCE
##############################################################

# Extract the diffraction order and the 11 reference wavelengths
# directly from every original Zemax configuration.

for i, config in enumerate(
    configs_original
):

    ##########################################################
    # Read diffraction order
    ##########################################################

    orders_original[i] = int(
        ln.zGetMulticon(
            config,
            row_order
        ).value
    )


    ##########################################################
    # Read WAVE 1 ... WAVE 11
    ##########################################################

    for j, row in enumerate(
        range(
            row_wave_first,
            row_wave_last + 1
        )
    ):

        waves_original[
            i,
            j
        ] = ln.zGetMulticon(
            config,
            row
        ).value


##############################################################
# Calculate m*lambda for original configurations
##############################################################

# orders_original has shape:
#
#     (13,)
#
# Adding a new axis converts it to:
#
#     (13, 1)
#
# This allows NumPy broadcasting against waves_original:
#
#     (13, 1) * (13, 11)
#
# resulting in:
#
#     (13, 11)

m_lambda = (
    orders_original[:, None]
    * waves_original
)


##############################################################
# Reference m*lambda values
##############################################################

# For each of the 11 wavelength positions, calculate the mean
# m*lambda value across all 13 original configurations.
#
# The result contains 11 reference constants:
#
#     C_1, C_2, ..., C_11

m_lambda_reference = np.mean(
    m_lambda,
    axis=0
)


##############################################################
# Standard deviation of m*lambda
##############################################################

# Calculate the dispersion of m*lambda across the original
# configurations for every wavelength position.
#
# Values close to zero indicate that the relation
#
#     m * lambda_j = constant
#
# is numerically preserved by the original MCE data.

m_lambda_std = np.std(
    m_lambda,
    axis=0
)


##############################################################
# Verify invariance of m*lambda
##############################################################

print(
    "\n"
    "============================================================\n"
    " Original MCE: m*lambda verification\n"
    "============================================================"
)


for j, (value, std) in enumerate(
    zip(
        m_lambda_reference,
        m_lambda_std
    ),
    start=1
):

    print(
        f"WAVE {j:2d} | "
        f"<m*lambda> = {value:.9f} um | "
        f"std = {std:.3e}"
    )


##############################################################
# Reconstruct wavelengths for all orders
##############################################################

# For every diffraction order m and wavelength reference
# position j, calculate:
#
#                   C_j
#     lambda_j(m) = -----
#                    m
#
# m_lambda_reference has shape:
#
#     (11,)
#
# orders_complete[:, None] has shape:
#
#     (85, 1)
#
# Broadcasting therefore generates a matrix with shape:
#
#     (85 orders, 11 wavelengths)

waves_complete = (
    m_lambda_reference[None, :]
    / orders_complete[:, None]
)


##############################################################
# Identify original and reconstructed orders
##############################################################

# Determine which integer orders between 60 and 144 are absent
# from the original 13-configuration Zemax model.

orders_missing = np.setdiff1d(
    orders_complete,
    orders_original
)


##############################################################
# Print order summary
##############################################################

print(
    "\n"
    "============================================================\n"
    " Echelle order summary\n"
    "============================================================"
)


print(
    "\nOriginal orders:\n",
    orders_original
)


print(
    "\nMissing orders:\n",
    orders_missing
)


print(
    f"\nOriginal configurations : "
    f"{len(orders_original)}"
)


print(
    f"Reconstructed orders    : "
    f"{len(orders_missing)}"
)


print(
    f"Total orders            : "
    f"{len(orders_complete)}"
)


##############################################################
# Close Zemax connection
##############################################################

# All required MCE information has now been extracted.
#
# No modifications were made to the optical model.

ln.close()


##############################################################
# Output file
##############################################################

# RESULTS_DIR is defined centrally in utils/emar_utils.py.
#
# This guarantees that the output is written to:
#
#     EMAR/Results/
#
# independently of the current working directory.

output_file = (
    emar_utils.RESULTS_DIR
    / "echelle_orders_60_144.csv"
)


# Make sure the Results directory exists.

emar_utils.RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True
)


##############################################################
# Write complete wavelength table to CSV
##############################################################

# CSV structure:
#
# order | status | wave_1 | wave_2 | ... | wave_11
#
# status is:
#
#     ORIGINAL
#
# if the diffraction order already exists in the original Zemax
# model, or:
#
#     NEW
#
# if it was reconstructed from the m*lambda relation.

with open(
    output_file,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(
        f
    )


    ##########################################################
    # CSV header
    ##########################################################

    header = [
        "order",
        "status"
    ] + [
        f"wave_{j}"
        for j in range(
            1,
            n_waves + 1
        )
    ]


    writer.writerow(
        header
    )


    ##########################################################
    # CSV data
    ##########################################################

    for m, wavelengths in zip(
        orders_complete,
        waves_complete
    ):

        # Mark whether the order comes directly from the original
        # Zemax MCE or was reconstructed.

        status = (
            "ORIGINAL"
            if m in orders_original
            else "NEW"
        )


        # Store the diffraction order, its status, and the
        # corresponding 11 wavelengths.

        row = [
            int(m),
            status
        ] + [
            f"{wavelength:.12f}"
            for wavelength in wavelengths
        ]


        writer.writerow(
            row
        )


##############################################################
# Confirmation
##############################################################

print(
    "\n"
    "============================================================\n"
    " Output\n"
    "============================================================"
)


print(
    "\nComplete echelle-order table saved to:"
)


print(
    output_file.resolve()
)