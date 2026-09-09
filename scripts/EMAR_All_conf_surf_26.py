# -*- coding: utf-8 -*-
"""
EMAR - All Original Configurations at Surface 26
================================================

This script traces the complete set of original EMAR Zemax
multi-configurations to surface 26.

The original optical model contains 13 representative configurations,
each associated with a different echelle diffraction order.

For every configuration:

1. The wavelength interval is read directly from the Zemax
   Multi-Configuration Editor (MCE).
2. A denser spectral sampling is generated between the first and
   last wavelength of the configuration.
3. The same normalized pupil sampling is used for every wavelength
   and every configuration.
4. Rays are propagated through the optical system up to surface 26.
5. The complete X-Y ray footprints are stored for later visualization.

Surface 26
----------
Surface 26 is used to inspect the spectral dispersion produced by the
echelle stage before the following optical block.

At this surface, the different configurations are expected to remain
geometrically registered in a similar region while each order spans its
corresponding wavelength interval.

The script can generate three different representations:

1. Scientific footprint
   Displays all pupil-ray intersections for every wavelength.

2. Spectral trace
   Calculates the centroid of each wavelength footprint and connects
   consecutive wavelength centroids.

3. Detector representation
   Uses the same centroid traces, but displays them as white continuous
   lines on a black background.

Reusable ray-tracing functions are implemented in:

    utils/emar_utils.py

Visualization functions are implemented in:

    utils/emar_plots.py

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

import pyzdde.zdde as pyz


##############################################################
# EMAR utilities
##############################################################

from utils import emar_utils
from utils import emar_plots


##############################################################
# Analysis parameters
##############################################################

# Original Zemax configurations.
#
# The current EMAR model contains 13 representative echelle
# configurations.

configs = range(
    1,
    14
)


# Surface where ray coordinates will be evaluated.

surf = 26


# Number of wavelengths generated between WAVE 1 and WAVE 11
# for every configuration.

n_wavelengths = 100


# Number of normalized pupil rays traced for every wavelength.

n_rays = 100


# Pupil-sampling method.
#
# Available options currently implemented in emar_utils:
#
#     "random"
#     "hexapolar"

pupil_sampling = "random"


# Random-number seed used only for random pupil sampling.
#
# A fixed seed guarantees that the same pupil coordinates are
# generated every time the experiment is executed.

seed = 123


# Number of hexapolar rings.
#
# This parameter is ignored when pupil_sampling="random".
#
# n_rings = 5 produces 91 pupil rays.

n_rings = 5


##############################################################
# Zemax optical model
##############################################################

# ZEMAX_DIR is defined centrally in utils/emar_utils.py.
#
# This avoids dependence on Spyder's current working directory.

zemax_file = (
    emar_utils.ZEMAX_DIR
    / "WP - Con prismas diseñados - camara.zmx"
)


##############################################################
# Connect to Zemax
##############################################################

# Create communication link between Python and Zemax OpticStudio.

ln = pyz.createLink()


# Load the EMAR optical model.

ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Trace all original configurations
##############################################################

# trace_configurations() performs the complete wavelength and pupil
# sampling internally.
#
# For each configuration it:
#
#     - activates the MCE configuration,
#     - reads WAVE 1 and WAVE 11,
#     - generates n_wavelengths intermediate wavelengths,
#     - temporarily modifies WAVE 1,
#     - updates Zemax,
#     - traces the complete pupil,
#     - stores the X-Y footprints,
#     - and restores the original MCE wavelength afterwards.
#
# The same pupil sampling is reused for all wavelengths and
# configurations.

spots_by_config = emar_utils.trace_configurations(
    ln=ln,
    configs=configs,
    surf=surf,
    n_wavelengths=n_wavelengths,
    n_rays=n_rays,
    pupil_sampling=pupil_sampling,
    n_rings=n_rings,
    seed=seed
)


##############################################################
# 1. Scientific footprint
##############################################################

# Display every valid pupil-ray intersection for every wavelength
# and configuration.
#
# The configurations follow the spectral ordering:
#
#     Config 1  -> red
#     ...
#     Config 13 -> blue

emar_plots.plot_all_footprints(
    spots_by_config=spots_by_config,
    surf=surf,
    color_mode="wavelength",
    view_mode="scatter"
)


##############################################################
# 2. Spectral centroid traces
##############################################################

# For every wavelength footprint:
#
#     Xc = mean(X)
#     Yc = mean(Y)
#
# Consecutive wavelength centroids are connected to form the
# continuous spectral trace of each configuration.
#
# The configuration-dependent spectral colors are preserved.

emar_plots.plot_all_footprints(
    spots_by_config=spots_by_config,
    surf=surf,
    color_mode="wavelength",
    view_mode="spectral_trace",
    linewidth=2.5
)

##############################################################
# Close Zemax connection
##############################################################

# All ray-tracing results are already stored in Python, so the
# communication link can be safely closed.

ln.close()