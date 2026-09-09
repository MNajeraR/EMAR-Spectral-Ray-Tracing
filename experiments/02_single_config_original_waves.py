# -*- coding: utf-8 -*-
"""
EMAR - Single Configuration with Original Zemax Wavelengths
============================================================

This script provides a first inspection of the spectral behavior of the
EMAR optical model using the wavelengths already defined in Zemax.

A single Zemax multi-configuration is selected and the chief ray is traced
for each of the 11 wavelengths defined in the wavelength table.

The ray is propagated up to a selected surface in order to inspect how its
position changes as a function of wavelength.

For the current example:

    Configuration = 2
    Surface       = 26

Surface 26 corresponds to the region of the optical system where the
spectral dispersion produced by the Echelle grating can be inspected.

Only the chief ray is traced in this script:

    Hx = 0
    Hy = 0
    Px = 0
    Py = 0

Therefore, the resulting coordinates represent the position of the chief
ray for each wavelength and not the complete geometrical footprint of the
beam.

The script produces three diagnostic plots:

1. Chief-ray positions in the XY plane.
2. X position as a function of wavelength.
3. Y position as a function of wavelength.

This script represents the first step in the EMAR spectral-tracing
development. Later scripts extend this analysis by:

- sampling the complete pupil,
- generating intermediate wavelengths,
- modifying the Zemax multi-configuration wavelength values,
- and tracing multiple configurations.

Author
------
Morgan R. Najera

Date
----
September 2026
"""


##############################################################
# Imports
##############################################################

import matplotlib.pyplot as plt
import pyzdde.zdde as pyz
# EMAR project path
from pathlib import Path
import sys
PROJECT_DIR = Path(__file__).resolve().parent.parent
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(
        0,
        str(PROJECT_DIR)
    )
from utils.emar_utils import ZEMAX_DIR


##############################################################
# Analysis parameters
##############################################################

# Zemax multi-configuration to analyze.
config = 2

# Surface where ray coordinates will be evaluated.
surf = 26

# Number of wavelengths already defined in the Zemax model.
n_waves = 11


##############################################################
# Zemax model
##############################################################

# Optical model location.
#
# ZEMAX_DIR is defined in utils/emar_utils.py and points to the
# Zemax directory located at the root of the EMAR project.

file_path = (
    ZEMAX_DIR
    / "WP - Con prismas diseñados - camara.zmx"
)


##############################################################
# Connect to Zemax
##############################################################

# Create the PyZDDE communication link with Zemax.

ln = pyz.createLink()


# Load the EMAR optical model.

ln.zLoadFile(
    str(file_path)
)


##############################################################
# Select Zemax configuration
##############################################################

# Activate the requested column of the Multi-Configuration Editor.

ln.zSetConfig(
    config
)

# Update the optical system after changing configuration.

ln.zGetUpdate()


##############################################################
# Storage arrays
##############################################################

# These lists store the wavelength and the corresponding chief-ray
# coordinates at the selected surface.

wavelengths = []

x_position = []

y_position = []


##############################################################
# Trace the original Zemax wavelengths
##############################################################

# The Zemax model contains 11 wavelengths.
#
# Here we use them exactly as they are defined in the optical model.
# No wavelength values are modified in this script.

for waveNum in range(
    1,
    n_waves + 1
):

    ##########################################################
    # Read wavelength from Zemax
    ##########################################################

    # zGetWave() returns the wavelength and its corresponding
    # wavelength weight.

    wave_data = ln.zGetWave(
        waveNum
    )

    wavelength = wave_data.wavelength


    ##########################################################
    # Trace chief ray
    ##########################################################

    # zGetTrace() traces a sequential ray through the optical
    # system up to the requested surface.
    #
    # Here:
    #
    #     Hx = Hy = 0
    #
    # corresponds to the central field, while:
    #
    #     Px = Py = 0
    #
    # selects the center of the normalized pupil. Therefore,
    # this is the chief ray for the central field.

    ret = ln.zGetTrace(
        waveNum,
        0,
        surf,
        0.0,       # Hx
        0.0,       # Hy
        0.0,       # Px
        0.0        # Py
    )


    ##########################################################
    # Ray-trace status
    ##########################################################

    error = ret[0]

    vignette = ret[1]


    ##########################################################
    # Store valid ray
    ##########################################################

    if error == 0:

        # Ray coordinates at the requested surface.

        x = ret[2]

        y = ret[3]


        # Store wavelength and position.

        wavelengths.append(
            wavelength
        )

        x_position.append(
            x
        )

        y_position.append(
            y
        )


        ######################################################
        # Diagnostic output
        ######################################################

        print(
            f"Wave {waveNum:2d} | "
            f"{wavelength:.6f} um | "
            f"X = {x:.6f} mm | "
            f"Y = {y:.6f} mm | "
            f"Vig = {vignette}"
        )


##############################################################
# Close Zemax connection
##############################################################

ln.close()


##############################################################
# Plot 1: Chief-ray positions in the XY plane
##############################################################

# This plot shows the geometrical trajectory followed by the
# chief-ray positions as wavelength changes.

plt.figure(
    figsize=(8, 8)
)

sc = plt.scatter(
    x_position,
    y_position,
    c=wavelengths,
    s=60
)

plt.xlabel(
    "X [mm]",
    fontsize=14
)

plt.ylabel(
    "Y [mm]",
    fontsize=14
)

plt.axis(
    "equal"
)


##############################################################
# Wavelength colorbar
##############################################################

cbar = plt.colorbar(
    sc
)

cbar.set_label(
    r"Wavelength [$\mu$m]"
)


plt.title(
    f"Configuration {config} - Surface {surf}"
)

plt.tight_layout()

plt.show()


##############################################################
# Plot 2: X position vs wavelength
##############################################################

# This plot isolates the spectral displacement along X.
#
# It is particularly useful for identifying the primary
# dispersion direction at the selected surface.

plt.figure(
    figsize=(9, 7)
)

plt.plot(
    wavelengths,
    x_position,
    "o-"
)

plt.xlabel(
    r"Wavelength [$\mu$m]",
    fontsize=14
)

plt.ylabel(
    "X position [mm]",
    fontsize=14
)

plt.tight_layout()

plt.show()


##############################################################
# Plot 3: Y position vs wavelength
##############################################################

# This plot shows the wavelength-dependent displacement along Y.
#
# Comparing this curve with the X-position curve helps determine
# whether the spectral trace is purely linear along X or also
# contains a wavelength-dependent curvature in Y.

plt.figure(
    figsize=(9, 7)
)

plt.plot(
    wavelengths,
    y_position,
    "o-"
)

plt.xlabel(
    r"Wavelength [$\mu$m]",
    fontsize=14
)

plt.ylabel(
    "Y position [mm]",
    fontsize=14
)

plt.tight_layout()

plt.show()