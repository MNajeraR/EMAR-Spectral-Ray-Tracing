# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
emar_utils.py
=============

Utility functions for ray-tracing studies of the EMAR white-pupil
echelle spectrograph using Zemax OpticStudio through PyZDDE.

This module contains reusable functions for:

1. Generating a reproducible random sampling of the normalized pupil.
2. Tracing a set of pupil rays for a selected wavelength and surface.
3. Reading the spectral limits of a Zemax multi-configuration.
4. Temporarily modifying one wavelength in the Multi-Configuration
   Editor (MCE) to generate a denser wavelength sampling.
5. Tracing several configurations using the same pupil sampling.
6. Calculating the centroid of the ray footprint for each wavelength.
7. Displaying the progress of long ray-tracing calculations.

The functions in this module do not define a particular experiment.
Parameters such as the Zemax file, configurations, target surface,
number of wavelengths, and number of rays are specified in the
main scripts.

The current Zemax model represents different echelle diffraction
orders using different configurations. Within each configuration,
WAVE 1 and WAVE 11 define the spectral limits used here for the
continuous wavelength sampling.

Author
------
Morgan Rhaí Nájera Roa
"""

import numpy as np
from pathlib import Path


##############################################################
# EMAR project directories
##############################################################

PROJECT_DIR = Path(__file__).resolve().parent.parent

ZEMAX_DIR = (
    PROJECT_DIR
    / "Zemax"
)

RESULTS_DIR = (
    PROJECT_DIR
    / "Results"
)


##############################################################
# Random pupil sampling
##############################################################

def sample_random_pupil(
    n_rays,
    seed=123
):
    """
    Generate random ray coordinates over a normalized circular pupil.

    The Zemax sequential ray-tracing function zGetTrace() uses normalized
    pupil coordinates (Px, Py). Valid points satisfy

        Px^2 + Py^2 <= 1.

    A uniform random distribution in radius would artificially place too
    many rays near the center of the pupil. To obtain a uniform sampling
    in pupil AREA, the radial coordinate is generated as

        r = sqrt(u),

    where u is uniformly distributed between 0 and 1.

    The same random seed can be used for every wavelength and
    configuration. This is useful because all wavelengths are then traced
    through exactly the same pupil coordinates, making chromatic
    comparisons easier to interpret.

    Parameters
    ----------
    n_rays : int
        Number of rays to generate over the normalized pupil.

    seed : int, optional
        Seed used by NumPy's random-number generator.
        The default is 123.

    Returns
    -------
    px : ndarray
        Normalized pupil X coordinates.

    py : ndarray
        Normalized pupil Y coordinates.
    """

    # Reproducible random-number generator
    rng = np.random.default_rng(seed)

    # Uniform random variables
    u = rng.random(n_rays)
    theta = 2.0 * np.pi * rng.random(n_rays)

    # sqrt(u) gives a uniform density over the pupil area
    r = np.sqrt(u)

    # Convert polar pupil coordinates to Cartesian coordinates
    px = r * np.cos(theta)
    py = r * np.sin(theta)

    return px, py

##############################################################
# Hexapolar pupil sampling
##############################################################

def sample_hexapolar_pupil(
    n_rings,
    include_center=True
):
    """
    Generate a deterministic hexapolar sampling of the normalized pupil.

    The pupil is divided into concentric rings. For ring ``i``, the
    normalized radial coordinate is

        r_i = i / n_rings

    and the number of azimuthal samples is

        N_i = 6 * i.

    Therefore, the number of rays increases with pupil radius:

        Ring 1 ->  6 rays
        Ring 2 -> 12 rays
        Ring 3 -> 18 rays
        ...

    This produces a structured approximately uniform sampling of the
    circular pupil that is useful for geometrical ray-tracing studies.

    Unlike random pupil sampling, the hexapolar pattern is completely
    deterministic and therefore does not require a random seed.

    Parameters
    ----------
    n_rings : int
        Number of concentric pupil rings.

    include_center : bool, optional
        If True, include one ray at the pupil center (Px, Py) = (0, 0).
        Default is True.

    Returns
    -------
    px : ndarray
        Normalized pupil X coordinates.

    py : ndarray
        Normalized pupil Y coordinates.

    Notes
    -----
    The total number of rays is

        N = 1 + 3 * n_rings * (n_rings + 1)

    when the central ray is included.
    """

    px = []
    py = []


    # --------------------------------------------------------
    # Central pupil ray
    # --------------------------------------------------------

    if include_center:

        px.append(0.0)
        py.append(0.0)


    # --------------------------------------------------------
    # Generate concentric hexapolar rings
    # --------------------------------------------------------

    for ring in range(
        1,
        n_rings + 1
    ):

        # Normalized radius of this ring
        r = (
            ring
            / n_rings
        )

        # Number of azimuthal rays in this ring
        n_theta = (
            6 * ring
        )

        # Equally spaced azimuthal coordinates
        theta = np.linspace(
            0.0,
            2.0 * np.pi,
            n_theta,
            endpoint=False
        )

        # Convert polar coordinates to normalized pupil coordinates
        px_ring = (
            r
            * np.cos(theta)
        )

        py_ring = (
            r
            * np.sin(theta)
        )

        px.extend(
            px_ring
        )

        py.extend(
            py_ring
        )


    return (
        np.asarray(px),
        np.asarray(py)
    )

##############################################################
# Trace one wavelength through the pupil
##############################################################

def trace_pupil(
    ln,
    wave_num,
    surf,
    px,
    py,
    hx=0.0,
    hy=0.0
):
    """
    Trace a set of pupil rays to a selected Zemax surface.

    For one wavelength and one field point, this function traces all
    normalized pupil coordinates supplied in ``px`` and ``py``.

    In the EMAR studies performed so far, the field is normally the
    central field:

        Hx = 0
        Hy = 0

    while Px and Py are varied to sample the complete pupil.

    The returned X and Y coordinates therefore describe the ray
    footprint produced by that wavelength at the selected surface.

    Examples of its interpretation are:

    - At an image plane, the coordinates form a spot/line image.
    - At or near a pupil plane, they describe the footprint of the beam.
    - By repeating the calculation for many wavelengths, the chromatic
      evolution of the footprint can be studied.

    Parameters
    ----------
    ln : PyZDDE link
        Active connection to Zemax OpticStudio.

    wave_num : int
        Zemax wavelength number used by zGetTrace().
        The current EMAR scripts normally use WAVE 1 and modify its value
        through the Multi-Configuration Editor.

    surf : int
        Zemax surface number at which the ray coordinates are requested.

    px, py : array_like
        Normalized pupil coordinates.

    hx, hy : float, optional
        Normalized field coordinates.
        The default is the central field: (0, 0).

    Returns
    -------
    x : ndarray
        X coordinates [mm] of all valid rays at the selected surface.

    y : ndarray
        Y coordinates [mm] of all valid rays at the selected surface.
    """

    x = []
    y = []

    # Trace one ray for every sampled pupil coordinate
    for pxi, pyi in zip(px, py):

        ret = ln.zGetTrace(
            wave_num,
            0,          # Sequential ray-tracing mode
            surf,
            hx,
            hy,
            float(pxi),
            float(pyi)
        )

        # PyZDDE returns error and vignette flags in the first entries
        error = ret[0]
        vignette = ret[1]

        # Store only rays successfully propagated to the surface
        if error == 0 and vignette == 0:

            x.append(ret[2])
            y.append(ret[3])

    return (
        np.asarray(x),
        np.asarray(y)
    )


##############################################################
# Read spectral limits from one Zemax configuration
##############################################################

def get_config_spectral_range(
    ln,
    config,
    row_wave_first=2,
    row_wave_last=12
):
    """
    Read the wavelength limits assigned to one MCE configuration.

    In the current EMAR Zemax model:

        MCE row 1  -> echelle diffraction-order parameter
        MCE row 2  -> WAVE 1
        ...
        MCE row 12 -> WAVE 11

    WAVE 1 and WAVE 11 are therefore used as the lower and upper
    wavelength limits of the configuration.

    The original model contains 11 wavelength samples per configuration.
    For the numerical experiments, these limits are used to construct a
    denser wavelength grid with ``numpy.linspace``.

    Parameters
    ----------
    ln : PyZDDE link
        Active Zemax connection.

    config : int
        Zemax configuration number.

    row_wave_first : int, optional
        MCE row containing the first wavelength.
        Default is row 2.

    row_wave_last : int, optional
        MCE row containing the last wavelength.
        Default is row 12.

    Returns
    -------
    wave_min : float
        First wavelength of the configuration [micrometers].

    wave_max : float
        Last wavelength of the configuration [micrometers].
    """

    wave_min = ln.zGetMulticon(
        config,
        row_wave_first
    ).value

    wave_max = ln.zGetMulticon(
        config,
        row_wave_last
    ).value

    return (
        wave_min,
        wave_max
    )


##############################################################
# Console progress bar
##############################################################

def print_progress_bar(
    completed,
    total,
    config,
    config_index,
    n_configs,
    wavelength_index,
    n_wavelengths,
    bar_length=30
):
    """
    Display the global progress of a multi-configuration ray trace.

    The same console line is continuously overwritten using ``\\r``.
    This avoids printing one line for every wavelength, which becomes
    impractical when hundreds or thousands of wavelength sets are traced.

    The percentage corresponds to the complete experiment, not only to
    the current configuration.

    Parameters
    ----------
    completed : int
        Number of wavelength sets already traced.

    total : int
        Total number of wavelength sets in the experiment.

    config : int
        Current Zemax configuration identifier.

    config_index : int
        Sequential position of the current configuration in the requested
        list.

    n_configs : int
        Total number of configurations.

    wavelength_index : int
        Current wavelength position inside the configuration.

    n_wavelengths : int
        Number of sampled wavelengths per configuration.

    bar_length : int, optional
        Number of characters used by the progress bar.
    """

    fraction = completed / total
    percent = 100.0 * fraction

    filled = int(
        bar_length * fraction
    )

    bar = (
        "#" * filled
        + "-" * (bar_length - filled)
    )

    print(
        f"\r"
        f"[{bar}] "
        f"{percent:6.2f}% | "
        f"Config {config_index:2d}/{n_configs:2d} "
        f"(ID {config:2d}) | "
        f"Wave {wavelength_index:3d}/{n_wavelengths:3d}",
        end="",
        flush=True
    )


##############################################################
# Trace one complete Zemax configuration
##############################################################

def trace_configuration(
    ln,
    config,
    surf,
    px,
    py,
    n_wavelengths=100,
    row_wave_first=2,
    row_wave_last=12,
    wave_num=1,
    hx=0.0,
    hy=0.0,
    progress_callback=None
):
    """
    Trace a dense wavelength sampling for one Zemax configuration.

    Each EMAR configuration represents one selected echelle diffraction
    order and contains 11 wavelengths in the Multi-Configuration Editor.

    Instead of tracing only those 11 discrete values, this function:

    1. Activates the requested configuration.
    2. Reads WAVE 1 and WAVE 11 from the MCE.
    3. Generates ``n_wavelengths`` values between those limits.
    4. Temporarily replaces the value of WAVE 1 in the MCE.
    5. Updates Zemax.
    6. Traces the complete pupil for that wavelength.
    7. Stores the resulting X-Y footprint.
    8. Repeats the process over the entire spectral interval.
    9. Restores the original WAVE 1 value when finished.

    IMPORTANT
    ---------
    The wavelength must be modified directly in the MCE.

    A simple call such as

        ln.zSetWave(...)

    is not sufficient for this model because the active
    Multi-Configuration operand overwrites the global wavelength value
    when Zemax is updated.

    The ``try/finally`` block guarantees that the original MCE value is
    restored even if the ray tracing fails during the calculation.

    Parameters
    ----------
    ln : PyZDDE link
        Active connection to Zemax.

    config : int
        Configuration to trace.

    surf : int
        Target Zemax surface.

    px, py : array_like
        Normalized pupil coordinates shared by all wavelengths.

    n_wavelengths : int, optional
        Number of wavelengths sampled between WAVE 1 and WAVE 11.

    row_wave_first, row_wave_last : int, optional
        MCE rows containing the spectral limits.

    wave_num : int, optional
        Zemax wavelength number passed to zGetTrace().
        Default is 1.

    hx, hy : float, optional
        Normalized field coordinates.

    progress_callback : callable or None
        Optional function called after each wavelength is traced.

    Returns
    -------
    result : dict
        Dictionary containing:

        ``wavelengths``
            Sampled wavelength array.

        ``spots``
            X-Y ray coordinates for every wavelength.

        ``wave_min``
            Original lower spectral limit.

        ``wave_max``
            Original upper spectral limit.
    """

    # --------------------------------------------------------
    # Activate the requested Zemax configuration
    # --------------------------------------------------------

    ln.zSetConfig(config)
    ln.zGetUpdate()


    # --------------------------------------------------------
    # Read the spectral interval defined in the MCE
    # --------------------------------------------------------

    wave_min, wave_max = get_config_spectral_range(
        ln=ln,
        config=config,
        row_wave_first=row_wave_first,
        row_wave_last=row_wave_last
    )


    # --------------------------------------------------------
    # Create a denser continuous wavelength sampling
    # --------------------------------------------------------

    wavelengths = np.linspace(
        wave_min,
        wave_max,
        n_wavelengths
    )


    # --------------------------------------------------------
    # Save the complete original MCE information for WAVE 1
    #
    # Besides the value itself, this stores status, pickup,
    # scale, and offset information so the cell can later be
    # restored exactly.
    # --------------------------------------------------------

    original_wave = ln.zGetMulticon(
        config,
        row_wave_first
    )


    # Dictionary used to store one footprint per wavelength
    spots_by_wave = {}


    try:

        # ----------------------------------------------------
        # Wavelength loop
        # ----------------------------------------------------

        for wave_index, wavelength in enumerate(
            wavelengths,
            start=1
        ):

            # ------------------------------------------------
            # Replace WAVE 1 directly in the active MCE
            # ------------------------------------------------

            ln.zSetMulticon(
                config,
                row_wave_first,
                float(wavelength),
                original_wave.status,
                original_wave.pickupRow,
                original_wave.pickupConfig,
                original_wave.scale,
                original_wave.offset
            )

            # Recalculate the optical system using the new λ
            ln.zGetUpdate()


            # ------------------------------------------------
            # Trace all pupil rays to the requested surface
            # ------------------------------------------------

            x, y = trace_pupil(
                ln=ln,
                wave_num=wave_num,
                surf=surf,
                px=px,
                py=py,
                hx=hx,
                hy=hy
            )


            # ------------------------------------------------
            # Store the complete ray footprint
            # ------------------------------------------------

            spots_by_wave[
                float(wavelength)
            ] = {
                "x": x.copy(),
                "y": y.copy()
            }


            # ------------------------------------------------
            # Update console progress
            # ------------------------------------------------

            if progress_callback is not None:

                progress_callback(
                    wave_index
                )


    finally:

        # ----------------------------------------------------
        # Restore the original value of WAVE 1
        #
        # This is performed even if an exception occurs.
        # ----------------------------------------------------

        ln.zSetMulticon(
            config,
            row_wave_first,
            original_wave.value,
            original_wave.status,
            original_wave.pickupRow,
            original_wave.pickupConfig,
            original_wave.scale,
            original_wave.offset
        )

        ln.zGetUpdate()


    return {
        "wavelengths": wavelengths.copy(),
        "spots": spots_by_wave,
        "wave_min": wave_min,
        "wave_max": wave_max
    }


##############################################################
# Trace several Zemax configurations
##############################################################

def trace_configurations(
    ln,
    configs,
    surf,
    n_wavelengths=100,
    n_rays=100,
    pupil_sampling="random",
    n_rings=5,
    seed=123,
    row_wave_first=2,
    row_wave_last=12,
    wave_num=1,
    hx=0.0,
    hy=0.0
):
    
    """
    Trace multiple EMAR configurations to the same Zemax surface.

    The same pupil sampling is used for every wavelength and every
    configuration. This is essential for meaningful comparisons because
    differences between configurations then originate from the optical
    system rather than from different random pupil samples.

    For the current model, the original 13 configurations represent
    selected echelle diffraction orders.

    Parameters
    ----------
    ln : PyZDDE link
        Active Zemax connection.

    configs : iterable of int
        Configuration numbers to trace.

    surf : int
        Target Zemax surface.

    n_wavelengths : int, optional
        Number of wavelengths sampled per configuration.

    n_rays : int, optional
        Number of pupil rays traced per wavelength.

    seed : int, optional
        Random pupil seed.

    row_wave_first, row_wave_last : int, optional
        MCE wavelength-limit rows.

    wave_num : int, optional
        Zemax wavelength number used by zGetTrace().

    hx, hy : float, optional
        Normalized field coordinates.

    Returns
    -------
    spots_by_config : dict
        Nested dictionary containing the wavelengths and footprints
        obtained for every configuration.
    """

    # Convert range/array/etc. to an explicit Python list
    configs = list(configs)

    n_configs = len(configs)


    # --------------------------------------------------------
    # Generate the pupil only ONCE
    #
    # Every wavelength and configuration uses the same rays.
    # --------------------------------------------------------

    # --------------------------------------------------------
    # Generate pupil sampling
    # --------------------------------------------------------
    
    if pupil_sampling == "random":
    
        px, py = sample_random_pupil(
            n_rays=n_rays,
            seed=seed
        )
    
    
    elif pupil_sampling == "hexapolar":
    
        px, py = sample_hexapolar_pupil(
            n_rings=n_rings,
            include_center=True
        )
    
    
    else:
    
        raise ValueError(
            "Unknown pupil sampling method: "
            f"{pupil_sampling}. "
            "Available methods are 'random' and 'hexapolar'."
        )
    
    
    # Actual number of rays generated
    n_rays_actual = len(px)


    # Total wavelength sets to be traced
    total_steps = (
        n_configs
        * n_wavelengths
    )

    completed_steps = 0


    # Main output dictionary
    spots_by_config = {}


    # --------------------------------------------------------
    # Configuration loop
    # --------------------------------------------------------

    for config_index, config in enumerate(
        configs,
        start=1
    ):

        # ----------------------------------------------------
        # Callback used by trace_configuration() to update
        # the global progress bar.
        # ----------------------------------------------------

        def update_progress(
            wavelength_index,
            config_index=config_index,
            config=config
        ):

            nonlocal completed_steps

            completed_steps += 1

            print_progress_bar(
                completed=completed_steps,
                total=total_steps,
                config=config,
                config_index=config_index,
                n_configs=n_configs,
                wavelength_index=wavelength_index,
                n_wavelengths=n_wavelengths
            )


        # ----------------------------------------------------
        # Trace the complete wavelength range of this config
        # ----------------------------------------------------

        result = trace_configuration(
            ln=ln,
            config=config,
            surf=surf,
            px=px,
            py=py,
            n_wavelengths=n_wavelengths,
            row_wave_first=row_wave_first,
            row_wave_last=row_wave_last,
            wave_num=wave_num,
            hx=hx,
            hy=hy,
            progress_callback=update_progress
        )


        # Store the complete result
        spots_by_config[
            config
        ] = result


    # Move console output to a new line after progress bar
    print()


    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print(
            f"\nTrace completed:"
            f"\n  Surface               : {surf}"
            f"\n  Configurations        : {n_configs}"
            f"\n  Pupil sampling        : {pupil_sampling}"
            f"\n  Wavelengths/config    : {n_wavelengths}"
            f"\n  Rays/wavelength       : {n_rays_actual}"
            f"\n  Total wavelength sets : {total_steps}"
            f"\n  Total ray traces      : "
            f"{total_steps * n_rays_actual}"
        )


    return spots_by_config


##############################################################
# Calculate spectral centroids
##############################################################

def calculate_centroids(
    spots_by_config
):
    """
    Calculate the centroid of every wavelength footprint.

    For each wavelength, the centroid is defined as

        Xc = mean(X)
        Yc = mean(Y)

    where X and Y contain all valid pupil rays reaching the selected
    surface.

    The centroids provide a compact representation of the spectral trace
    and are useful for studying:

    - the main dispersion direction,
    - order curvature,
    - geometrical registration between configurations,
    - separation of echelle orders after cross dispersion.

    Parameters
    ----------
    spots_by_config : dict
        Output generated by ``trace_configurations()``.

    Returns
    -------
    centroids_by_config : dict
        Dictionary containing wavelength, X centroid, and Y centroid
        arrays for every configuration.
    """

    centroids_by_config = {}


    # --------------------------------------------------------
    # Loop over configurations
    # --------------------------------------------------------

    for config, config_data in spots_by_config.items():

        wavelengths = config_data[
            "wavelengths"
        ]

        spots_by_wave = config_data[
            "spots"
        ]

        xc = []
        yc = []


        # ----------------------------------------------------
        # Calculate centroid for every wavelength
        # ----------------------------------------------------

        for wavelength in wavelengths:

            spot = spots_by_wave[
                float(wavelength)
            ]

            x = spot["x"]
            y = spot["y"]


            # If no valid rays reach the selected surface,
            # store NaN rather than producing an invalid mean.
            if len(x) == 0:

                xc.append(np.nan)
                yc.append(np.nan)

            else:

                xc.append(
                    np.mean(x)
                )

                yc.append(
                    np.mean(y)
                )


        # ----------------------------------------------------
        # Store centroid arrays
        # ----------------------------------------------------

        centroids_by_config[
            config
        ] = {
            "wavelengths": np.asarray(
                wavelengths
            ),
            "xc": np.asarray(xc),
            "yc": np.asarray(yc)
        }


    return centroids_by_config

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


        x, y = trace_pupil(
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


            x, y = trace_pupil(
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



