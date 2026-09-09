# -*- coding: utf-8 -*-
"""
EMAR - Two Configurations with 100 Wavelengths
==============================================

This script extends the single-configuration spectral analysis to two
Zemax multi-configurations.

For the current example:

    Configurations      = 1 and 2
    Surface             = 26
    Spectral samples    = 100 wavelengths/configuration
    Pupil samples       = 100 rays/wavelength

Each configuration represents a different echelle diffraction order
with its own wavelength interval defined in the Zemax
Multi-Configuration Editor (MCE).

For every configuration:

1. The configuration is activated in Zemax.
2. WAVE 1 and WAVE 11 are read from the MCE.
3. 100 wavelengths are generated between those limits.
4. WAVE 1 is temporarily replaced by each sampled wavelength.
5. Zemax is updated.
6. The actual wavelength used by Zemax is verified.
7. The same normalized pupil sampling is traced to the selected surface.
8. The resulting X-Y ray footprint is stored.
9. The original WAVE 1 value is restored.

The two configurations are then compared through:

1. Complete X-Y ray footprints.
2. Spectral centroid traces.

This script represents the transition from the single-order analysis
to the multi-configuration analysis later generalized to all 13
original configurations.

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



import numpy as np
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

from utils import emar_utils
from utils import emar_plots


##############################################################
# Analysis parameters
##############################################################

# Zemax configurations to compare.

configs = [
    1,
    2
]


# Surface where ray coordinates will be evaluated.

surf = 26


# Number of wavelengths generated within each configuration.

n_wavelengths = 100


# Number of normalized pupil rays traced for every wavelength.

n_rays = 100


# Random seed used for the pupil sampling.
#
# The same pupil coordinates are used for both configurations so that
# differences in the resulting footprints originate from the optical
# system rather than from different random pupil realizations.

seed = 123


# MCE rows defining the wavelength interval.

row_wave1 = 2
row_wave11 = 12


##############################################################
# Zemax optical model
##############################################################

file_path = (
    emar_utils.ZEMAX_DIR
    / "WP - Con prismas diseñados - camara.zmx"
)


##############################################################
# Connect to Zemax
##############################################################

ln = pyz.createLink()

ln.zLoadFile(
    str(file_path)
)


##############################################################
# Generate pupil sampling
##############################################################

# Generate one random pupil distribution and reuse it for both
# configurations and all wavelengths.

px, py = emar_utils.sample_random_pupil(
    n_rays=n_rays,
    seed=seed
)


##############################################################
# Store results for all configurations
##############################################################

spots_by_config = {}


##############################################################
# Loop over configurations
##############################################################

for config in configs:

    ##########################################################
    # Activate configuration
    ##########################################################

    ln.zSetConfig(
        config
    )

    ln.zGetUpdate()


    ##########################################################
    # Read spectral limits from the MCE
    ##########################################################

    wave_min = ln.zGetMulticon(
        config,
        row_wave1
    ).value


    wave_max = ln.zGetMulticon(
        config,
        row_wave11
    ).value


    ##########################################################
    # Generate wavelength sampling
    ##########################################################

    wavelengths_um = np.linspace(
        wave_min,
        wave_max,
        n_wavelengths
    )


    print(
        f"\nConfiguration {config}: "
        f"{wave_min:.6f} - {wave_max:.6f} um"
    )


    ##########################################################
    # Save original WAVE 1 MCE data
    ##########################################################

    # WAVE 1 will be modified temporarily.
    #
    # The complete MCE cell is stored so that its value and solve
    # properties can be restored after the configuration is traced.

    original_wave1 = ln.zGetMulticon(
        config,
        row_wave1
    )


    ##########################################################
    # Store all wavelength footprints for this configuration
    ##########################################################

    spots_by_wave = {}


    ##########################################################
    # Loop over wavelengths
    ##########################################################

    try:

        for wavelength in wavelengths_um:

            ##################################################
            # Modify WAVE 1 directly in the MCE
            ##################################################

            ln.zSetMulticon(
                config,
                row_wave1,
                float(wavelength),
                original_wave1.status,
                original_wave1.pickupRow,
                original_wave1.pickupConfig,
                original_wave1.scale,
                original_wave1.offset
            )


            ##################################################
            # Update Zemax
            ##################################################

            ln.zGetUpdate()


            ##################################################
            # Verify actual wavelength used by Zemax
            ##################################################

            wave_zemax = ln.zGetWave(
                1
            ).wavelength


            ##################################################
            # Trace complete pupil
            ##################################################

            x, y = emar_utils.trace_pupil(
                ln=ln,
                wave_num=1,
                surf=surf,
                px=px,
                py=py,
                hx=0.0,
                hy=0.0
            )


            ##################################################
            # Store wavelength footprint
            ##################################################

            spots_by_wave[
                float(wavelength)
            ] = {
                "x": x.copy(),
                "y": y.copy()
            }


            ##################################################
            # Diagnostic output
            ##################################################

            # This output is intentionally kept in this experimental
            # script because one of its purposes is to verify that the
            # wavelength requested through the MCE is the wavelength
            # actually used by Zemax.

            print(
                f"Config {config:2d} | "
                f"Requested = {wavelength:.6f} um | "
                f"Zemax = {wave_zemax:.6f} um | "
                f"{len(x)}/{n_rays} valid rays"
            )


    finally:

        ######################################################
        # Restore original WAVE 1
        ######################################################

        # The original MCE value is restored even if an exception
        # occurs during the wavelength loop.

        ln.zSetMulticon(
            config,
            row_wave1,
            original_wave1.value,
            original_wave1.status,
            original_wave1.pickupRow,
            original_wave1.pickupConfig,
            original_wave1.scale,
            original_wave1.offset
        )

        ln.zGetUpdate()


    ##########################################################
    # Store complete configuration
    ##########################################################

    spots_by_config[
        config
    ] = {
        "wavelengths": wavelengths_um.copy(),
        "spots": spots_by_wave
    }


##############################################################
# Close Zemax connection
##############################################################

ln.close()


##############################################################
# Calculate spectral centroids
##############################################################

# Convert the complete pupil footprints into one centroid position
# (Xc, Yc) for every wavelength and configuration.

centroids_by_config = (
    emar_utils.calculate_centroids(
        spots_by_config
    )
)


##############################################################
# Plot complete footprints
##############################################################

# Each configuration is assigned a different color.
#
# The complete pupil footprints are shown for all 100 wavelengths
# belonging to each configuration.

emar_plots.plot_all_footprints(
    spots_by_config=spots_by_config,
    surf=surf,
    point_size=3,
    alpha=0.6,
    save=False
)


##############################################################
# Plot spectral centroid traces
##############################################################

# The centroid trace provides a compact representation of the spectral
# trajectory followed by each configuration.

emar_plots.plot_spectral_traces(
    centroids_by_config=centroids_by_config,
    surf=surf,
    linewidth=1.5,
    save=False
)