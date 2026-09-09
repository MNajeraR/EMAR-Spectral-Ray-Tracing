# -*- coding: utf-8 -*-

"""
emar_plots.py
=============

Plotting utilities for the EMAR white-pupil echelle spectrograph
ray-tracing studies.

This module contains reusable visualization functions for:

1. Plotting the complete X-Y ray footprints of several Zemax
   configurations at a selected surface.
2. Displaying continuous spectral centroid traces.
3. Displaying a detector-like representation of the spectral orders.
4. Visualizing normalized pupil sampling.
5. Saving figures to the project ``Results`` directory in PNG,
   PDF, or both formats.

Three spectral visualization modes are available:

Scientific footprint
--------------------
``view_mode="scatter"``

Displays all valid ray intersections for every sampled wavelength.

For the original EMAR configuration ordering:

    Config 1  -> red  -> longer wavelengths
    ...
    Config 13 -> blue -> shorter wavelengths


Spectral trace
--------------
``view_mode="spectral_trace"``

For every wavelength, the centroid of the complete pupil footprint is
calculated:

    Xc = mean(X)
    Yc = mean(Y)

The centroids are connected to form continuous spectral traces.

Each configuration retains its spectral color and the figure uses a
white background.


Detector view
-------------
``view_mode="detector"``

Uses the same centroid traces as the spectral-trace representation, but
displays all configurations as white lines on a black background.

This mode intentionally removes wavelength-dependent visible color and
is intended as a simplified representation of the spatial spectral
format recorded by a detector.

The plotting functions receive data already calculated by
``emar_utils.py`` and do not communicate directly with Zemax.

Author
------
Morgan Rhaí Nájera Roa
"""

from pathlib import Path

from utils.emar_utils import RESULTS_DIR

import matplotlib.pyplot as plt
import numpy as np


##############################################################
# Internal figure-saving utility
##############################################################

def save_figure(
    fig,
    filename,
    results_dir=RESULTS_DIR,
    file_format="png",
    dpi=300
):
    """
    Save a Matplotlib figure in the Results directory.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to save.

    filename : str
        Base filename without extension.

    results_dir : str or pathlib.Path, optional
        Directory where figures will be stored.

    file_format : {"png", "pdf", "both"}, optional
        Output format.

    dpi : int, optional
        Resolution used for raster output.

    Returns
    -------
    saved_files : list of pathlib.Path
        Paths of the files created.
    """

    results_dir = Path(
        results_dir
    )

    results_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    file_format = file_format.lower()


    if file_format == "png":

        formats = [
            "png"
        ]


    elif file_format == "pdf":

        formats = [
            "pdf"
        ]


    elif file_format == "both":

        formats = [
            "png",
            "pdf"
        ]


    else:

        raise ValueError(
            "Unknown file format: "
            f"{file_format}. "
            "Available formats are 'png', 'pdf', or 'both'."
        )


    saved_files = []


    for fmt in formats:

        output_file = (
            results_dir
            / f"{filename}.{fmt}"
        )

        fig.savefig(
            output_file,
            dpi=dpi,
            bbox_inches="tight"
        )

        saved_files.append(
            output_file
        )


    for output_file in saved_files:

        print(
            f"Figure saved: "
            f"{output_file.resolve()}"
        )


    return saved_files


##############################################################
# Configuration color mapping
##############################################################

def get_config_color(
    config,
    configs,
    cmap_name="coolwarm"
):
    """
    Assign a continuous spectral color to one configuration.

    The configuration ordering is reversed with respect to the
    colormap normalization so that:

        Config 1  -> red
        Config N  -> blue

    This follows the original EMAR model, where lower configuration
    numbers correspond to longer wavelengths.

    Parameters
    ----------
    config : int
        Configuration identifier.

    configs : iterable
        Complete list of configurations being plotted.

    cmap_name : str, optional
        Matplotlib colormap name.

    Returns
    -------
    color : tuple
        RGBA color assigned to the configuration.
    """

    configs = np.asarray(
        configs,
        dtype=float
    )

    config_min = np.min(
        configs
    )

    config_max = np.max(
        configs
    )


    if config_max == config_min:

        normalized = 0.5


    else:

        normalized = (
            config_max - config
        ) / (
            config_max - config_min
        )


    cmap = plt.get_cmap(
        cmap_name
    )


    return cmap(
        normalized
    )


##############################################################
# Calculate centroid trace from one configuration
##############################################################

def get_centroid_trace(
    wavelengths,
    spots_by_wave
):
    """
    Calculate the X-Y centroid trajectory of one configuration.

    For every sampled wavelength:

        Xc = mean(X)
        Yc = mean(Y)

    Wavelengths for which no valid rays reach the selected surface are
    skipped.

    Parameters
    ----------
    wavelengths : array_like
        Sampled wavelengths.

    spots_by_wave : dict
        Dictionary containing X-Y ray footprints for every wavelength.

    Returns
    -------
    xc : ndarray
        X centroid coordinates.

    yc : ndarray
        Y centroid coordinates.
    """

    xc = []
    yc = []


    for wavelength in wavelengths:

        spot = spots_by_wave[
            float(wavelength)
        ]

        x = spot[
            "x"
        ]

        y = spot[
            "y"
        ]


        if len(x) == 0:

            continue


        xc.append(
            np.mean(x)
        )

        yc.append(
            np.mean(y)
        )


    return (
        np.asarray(xc),
        np.asarray(yc)
    )


##############################################################
# Plot spectral results
##############################################################

def plot_all_footprints(
    spots_by_config,
    surf,
    figsize=(11, 8),
    point_size=2,
    alpha=0.35,
    color_mode="wavelength",
    view_mode="scatter",
    linewidth=1.5,
    detector_alpha=1.0,
    save=False,
    results_dir=RESULTS_DIR,
    file_format="png",
    dpi=300,
    show=True
):
    """
    Plot EMAR spectral results using one of three visualization modes.

    Available modes
    ---------------

    ``view_mode="scatter"``
        Plot the complete ray footprint for every wavelength.

        Each configuration is represented with its spectral color.

    ``view_mode="spectral_trace"``
        Calculate the centroid of each wavelength footprint and connect
        consecutive centroid positions.

        Each configuration retains its spectral color and is displayed
        on a white background.

    ``view_mode="detector"``
        Calculate the same centroid traces, but display all configurations
        as white lines on a black background.

    Parameters
    ----------
    spots_by_config : dict
        Output generated by ``emar_utils.trace_configurations()``.

    surf : int
        Zemax surface represented by the plot.

    figsize : tuple, optional
        Figure size in inches.

    point_size : float, optional
        Marker size used in scatter mode.

    alpha : float, optional
        Marker transparency used in scatter mode.

    color_mode : {"wavelength", "categorical"}, optional
        Color convention used in scatter and spectral-trace modes.

    view_mode : {"scatter", "spectral_trace", "detector"}, optional
        Visualization mode.

    linewidth : float, optional
        Width of spectral traces.

    detector_alpha : float, optional
        Opacity of white detector traces.

    save : bool, optional
        Save figure if True.

    results_dir : str or pathlib.Path, optional
        Output directory.

    file_format : {"png", "pdf", "both"}, optional
        Output format.

    dpi : int, optional
        Raster resolution.

    show : bool, optional
        Display figure if True.

    Returns
    -------
    fig : matplotlib.figure.Figure

    ax : matplotlib.axes.Axes
    """

    ##########################################################
    # Validate visualization mode
    ##########################################################

    valid_modes = [
        "scatter",
        "spectral_trace",
        "detector"
    ]


    if view_mode not in valid_modes:

        raise ValueError(
            "Unknown view_mode: "
            f"{view_mode}. "
            "Available modes are "
            "'scatter', 'spectral_trace', and 'detector'."
        )


    ##########################################################
    # Create figure
    ##########################################################

    fig, ax = plt.subplots(
        figsize=figsize
    )


    ##########################################################
    # Configuration identifiers
    ##########################################################

    configs = list(
        spots_by_config.keys()
    )


    ##########################################################
    # Detector appearance
    ##########################################################

    if view_mode == "detector":

        ax.set_facecolor(
            "black"
        )

        fig.patch.set_facecolor(
            "black"
        )


    ##########################################################
    # Loop over configurations
    ##########################################################

    for i, config in enumerate(
        configs
    ):

        wavelengths = spots_by_config[
            config
        ]["wavelengths"]

        spots_by_wave = spots_by_config[
            config
        ]["spots"]


        ######################################################
        # Configuration color
        ######################################################

        if view_mode == "detector":

            color = "white"


        elif color_mode == "wavelength":

            color = get_config_color(
                config=config,
                configs=configs,
                cmap_name="coolwarm"
            )


        elif color_mode == "categorical":

            color = plt.cm.tab20(
                i % 20
            )


        else:

            raise ValueError(
                "Unknown color_mode: "
                f"{color_mode}. "
                "Available modes are "
                "'wavelength' and 'categorical'."
            )


        ######################################################
        # Scatter footprint
        ######################################################

        if view_mode == "scatter":

            for wavelength in wavelengths:

                spot = spots_by_wave[
                    float(wavelength)
                ]

                x = spot[
                    "x"
                ]

                y = spot[
                    "y"
                ]


                ax.scatter(
                    x,
                    y,
                    s=point_size,
                    color=color,
                    alpha=alpha,
                    linewidths=0
                )


            # One legend entry per configuration
            ax.scatter(
                [],
                [],
                color=color,
                s=30,
                label=f"Config {config}"
            )


        ######################################################
        # Spectral trace
        ######################################################

        elif view_mode == "spectral_trace":

            xc, yc = get_centroid_trace(
                wavelengths=wavelengths,
                spots_by_wave=spots_by_wave
            )


            ax.plot(
                xc,
                yc,
                "-",
                color=color,
                linewidth=linewidth,
                label=f"Config {config}"
            )


        ######################################################
        # Detector trace
        ######################################################

        elif view_mode == "detector":

            xc, yc = get_centroid_trace(
                wavelengths=wavelengths,
                spots_by_wave=spots_by_wave
            )


            ax.plot(
                xc,
                yc,
                "-",
                color="white",
                linewidth=linewidth,
                alpha=detector_alpha
            )


    ##########################################################
    # Axis labels
    ##########################################################

    if view_mode == "spectral_trace":

        ax.set_xlabel(
            "X centroid [mm]",
            fontsize=14
        )

        ax.set_ylabel(
            "Y centroid [mm]",
            fontsize=14
        )


    else:

        ax.set_xlabel(
            "X [mm]",
            fontsize=14
        )

        ax.set_ylabel(
            "Y [mm]",
            fontsize=14
        )


    ##########################################################
    # Title
    ##########################################################

    if view_mode == "scatter":

        title = (
            f"All configurations - Surface {surf}"
        )


    elif view_mode == "spectral_trace":

        title = (
            f"Spectral traces - Surface {surf}"
        )


    else:

        title = (
            f"Detector spectral trace - Surface {surf}"
        )


    ax.set_title(
        title,
        fontsize=15
    )


    ##########################################################
    # Scientific-view formatting
    ##########################################################

    if view_mode in [
        "scatter",
        "spectral_trace"
    ]:

        ax.tick_params(
            axis="both",
            labelsize=12
        )


        ax.legend(
            fontsize=9,
            ncol=2
        )


    ##########################################################
    # Detector formatting
    ##########################################################

    elif view_mode == "detector":

        ax.xaxis.label.set_color(
            "white"
        )

        ax.yaxis.label.set_color(
            "white"
        )

        ax.title.set_color(
            "white"
        )


        ax.tick_params(
            axis="both",
            colors="white",
            labelsize=12
        )


        for spine in ax.spines.values():

            spine.set_color(
                "white"
            )


    ##########################################################
    # Layout
    ##########################################################

    fig.tight_layout()


    ##########################################################
    # Save figure
    ##########################################################

    if save:

        if view_mode == "scatter":

            filename = (
                f"Footprint_surf{surf}"
            )


        elif view_mode == "spectral_trace":

            filename = (
                f"Spectral_trace_surf{surf}"
            )


        else:

            filename = (
                f"Detector_spectral_trace_surf{surf}"
            )


        save_figure(
            fig=fig,
            filename=filename,
            results_dir=results_dir,
            file_format=file_format,
            dpi=dpi
        )


    ##########################################################
    # Display
    ##########################################################

    if show:

        plt.show()


    return fig, ax


##############################################################
# Plot normalized pupil sampling
##############################################################

def plot_pupil_sampling(
    px,
    py,
    title="Pupil sampling",
    figsize=(7, 7),
    point_size=30,
    save=False,
    results_dir=RESULTS_DIR,
    file_format="png",
    filename="Pupil_sampling",
    dpi=300,
    show=True
):
    """
    Visualize normalized pupil coordinates used for ray tracing.

    The normalized circular pupil satisfies

        Px^2 + Py^2 <= 1.

    This function is useful for comparing pupil-sampling strategies
    such as random area-uniform and hexapolar sampling.

    Parameters
    ----------
    px, py : array_like
        Normalized pupil coordinates.

    title : str, optional
        Figure title.

    figsize : tuple, optional
        Figure size.

    point_size : float, optional
        Scatter-marker size.

    save : bool, optional
        Save figure if True.

    results_dir : str or pathlib.Path, optional
        Output directory.

    file_format : {"png", "pdf", "both"}, optional
        Output format.

    filename : str, optional
        Base filename used when saving.

    dpi : int, optional
        Raster resolution.

    show : bool, optional
        Display figure if True.

    Returns
    -------
    fig, ax
    """

    fig, ax = plt.subplots(
        figsize=figsize
    )

    ##########################################################
    # Normalized pupil boundary
    ##########################################################

    theta = np.linspace(
        0.0,
        2.0 * np.pi,
        500
    )


    ax.plot(
        np.cos(theta),
        np.sin(theta),
        "-",
        color="k",
        linewidth=2
    )
    
    ##########################################################
    # Plot pupil coordinates
    ##########################################################

    ax.scatter(
        px,
        py,
        marker = "x",
        color="#0072B2",
        linewidths=1.5,
        s=point_size
    )

    ##########################################################
    # Labels
    ##########################################################

    ax.set_xlabel(
        r"$\rho_x$",
        fontsize=16
    )

    ax.set_ylabel(
        r"$\rho_y$",
        fontsize=16
    )


    ax.set_title(
        title,
        fontsize=14
    )


    ##########################################################
    # Preserve circular pupil geometry
    ##########################################################

    ax.set_aspect(
        "equal"
    )

    ax.set_xlim(
        -1.05,
        1.05
    )

    ax.set_ylim(
        -1.05,
        1.05
    )
    
    ##########################################################
    # Remove axis ticks
    ##########################################################
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    
    ##########################################################
    # Normalized scale bar
    ##########################################################
    
    # bar_x0 = -1.
    # bar_x1 = -0.0
    # bar_y  =  0.0
    
    # ax.plot(
    #     [bar_x0, bar_x1],
    #     [bar_y, bar_y],
    #     color="k",
    #     linewidth=2
    # )
    
    # ax.plot(
    #     [bar_x0, bar_x0],
    #     [bar_y - 0.025, bar_y + 0.025],
    #     color="k",
    #     linewidth=2
    # )
    
    # ax.plot(
    #     [bar_x1, bar_x1],
    #     [bar_y - 0.025, bar_y + 0.025],
    #     color="k",
    #     linewidth=2
    # )
    
    # ax.text(
    #     (bar_x0 + bar_x1) / 2,
    #     bar_y + 0.06,
    #     r"$0.5$",
    #     ha="center",
    #     va="bottom",
    #     fontsize=12
    # )


    ##########################################################
    # Grid
    ##########################################################

    # ax.grid(
        # alpha=0.3
    # )


    fig.tight_layout()


    ##########################################################
    # Save
    ##########################################################

    if save:

        save_figure(
            fig=fig,
            filename=filename,
            results_dir=results_dir,
            file_format=file_format,
            dpi=dpi
        )


    ##########################################################
    # Display
    ##########################################################

    if show:

        plt.show()


    return fig, ax