# -*- coding: utf-8 -*-

"""
12_full_fiber_spectrum.py
=========================

Full EMAR echelle spectral analysis including a 100 um
circular fiber.

Orders:
    m = 60 ... 144

For every diffraction order:

1. Read the wavelength limits from:
       Results/echelle_orders_60_144.csv

2. Generate 100 wavelengths between wave_1 and wave_11.

3. Sample the 100 um circular fiber with a deterministic
   5-ring hexapolar distribution:
       91 fiber points.

4. For every fiber point, trace all 100 wavelengths using
   a deterministic 5-ring hexapolar pupil:
       91 pupil rays.

5. For every wavelength measure:
       Xc, Yc
       Xmin, Xmax
       Ymin, Ymax
       Delta X
       Delta Y

6. Calculate the reciprocal linear dispersion along X:
       d(lambda)/dX

7. Maintain the geometrical definition:
       omega' = Delta X

8. Calculate:
       Delta lambda = omega' * d(lambda)/dX
       R = lambda / Delta lambda

9. Save each completed order immediately.

10. Produce:
       - Linear dispersion for all orders
       - Resolving power for all orders
       - Echelle centroid spectral format
       - Echelle format including the projected fiber size

Color convention:
       m = 60  -> red
       m = 144 -> blue

Author
------
Morgan Rhaí Nájera Roa
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import numpy as np
import pandas as pd
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


##############################################################
# Optical surfaces
##############################################################

surf_detector = 54


##############################################################
# Spectral sampling
##############################################################

n_wavelengths = 100


##############################################################
# Fiber geometry
##############################################################

fiber_diameter_um = 100.0


##############################################################
# Hexapolar sampling
##############################################################

n_fiber_rings = 5

n_pupil_rings = 5


##############################################################
# Calibrated field coordinates for the 100 um fiber
##############################################################

hx_max = 0.286480

hy_max = 0.285400


##############################################################
# MCE rows
##############################################################

row_order = 1

row_wave1 = 2


##############################################################
# Files
##############################################################

orders_file = (
    emar_utils.RESULTS_DIR
    / "echelle_orders_60_144.csv"
)


zemax_file = (
    emar_utils.ZEMAX_DIR
    / "WP - Con prismas diseñados - camara - theoretical slit.zmx"
)


output_file = (
    emar_utils.RESULTS_DIR
    / "full_orders_fiber_resolution.csv"
)


##############################################################
# Load complete echelle-order table
##############################################################

order_table = np.genfromtxt(
    orders_file,
    delimiter=",",
    names=True,
    dtype=None,
    encoding="utf-8"
)


##############################################################
# Orders
##############################################################

orders = np.asarray(
    order_table["order"],
    dtype=int
)


n_orders = len(
    orders
)


##############################################################
# Fiber sampling
##############################################################

fiber_x_norm, fiber_y_norm = (
    emar_utils.sample_hexapolar_pupil(
        n_rings=n_fiber_rings,
        include_center=True
    )
)


hx_fiber = (
    hx_max
    * fiber_x_norm
)


hy_fiber = (
    hy_max
    * fiber_y_norm
)


##############################################################
# Pupil sampling
##############################################################

px, py = (
    emar_utils.sample_hexapolar_pupil(
        n_rings=n_pupil_rings,
        include_center=True
    )
)


##############################################################
# Sampling information
##############################################################

n_fiber_points = len(
    hx_fiber
)


n_pupil_rays = len(
    px
)


n_rays_per_wavelength = (
    n_fiber_points
    * n_pupil_rays
)


n_rays_per_order = (
    n_rays_per_wavelength
    * n_wavelengths
)


n_total_rays = (
    n_rays_per_order
    * n_orders
)


##############################################################
# Progress bar
##############################################################

def print_progress(
    current,
    total,
    order,
    fiber_point,
    fiber_total,
    bar_length=40
):

    """
    Global progress bar over the complete 85-order calculation.
    """

    fraction = (
        current
        / total
    )


    filled_length = int(
        bar_length
        * fraction
    )


    bar = (
        "="
        * filled_length
        + "-"
        * (
            bar_length
            - filled_length
        )
    )


    percentage = (
        100.0
        * fraction
    )


    print(
        f"\r[{bar}] "
        f"{percentage:6.2f}% | "
        f"m={order:3d} | "
        f"fiber {fiber_point:2d}/{fiber_total}",
        end="",
        flush=True
    )


    if current == total:

        print()


##############################################################
# Storage
##############################################################

all_order_results = []


##############################################################
# Print experiment information
##############################################################

print()

print(
    "EMAR full 100 um fiber spectral analysis"
)

print(
    "----------------------------------------"
)


print(
    f"Orders:                    "
    f"{orders.min()} - {orders.max()}"
)


print(
    f"Number of orders:          "
    f"{n_orders}"
)


print(
    f"Wavelengths/order:         "
    f"{n_wavelengths}"
)


print()

print(
    f"Fiber diameter:            "
    f"{fiber_diameter_um:.1f} um"
)


print(
    f"Fiber rings:               "
    f"{n_fiber_rings}"
)


print(
    f"Fiber points:              "
    f"{n_fiber_points}"
)


print(
    f"Pupil rings:               "
    f"{n_pupil_rings}"
)


print(
    f"Pupil rays:                "
    f"{n_pupil_rays}"
)


print()

print(
    f"Rays/wavelength:           "
    f"{n_rays_per_wavelength:,}"
)


print(
    f"Rays/order:                "
    f"{n_rays_per_order:,}"
)


print(
    f"Maximum total rays:        "
    f"{n_total_rays:,}"
)


print()

print(
    "Tracing..."
)

print()


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
# Load Zemax model
##############################################################

ln.zLoadFile(
    str(zemax_file)
)


##############################################################
# Global progress counter
##############################################################

total_calls = (
    n_orders
    * n_fiber_points
)


completed_calls = 0


##############################################################
# Main calculation
##############################################################

try:

    ##########################################################
    # Loop over diffraction orders
    ##########################################################

    for order_index, row in enumerate(
        order_table,
        start=1
    ):

        ######################################################
        # Order information
        ######################################################

        order = int(
            row["order"]
        )


        status = row[
            "status"
        ]


        ######################################################
        # Reference wavelengths
        ######################################################

        reference_wavelengths = np.asarray(
            [
                row[
                    f"wave_{j}"
                ]
                for j in range(
                    1,
                    12
                )
            ],
            dtype=float
        )


        ######################################################
        # Dense wavelength sampling
        ######################################################

        wavelengths_um = np.linspace(
            reference_wavelengths[0],
            reference_wavelengths[-1],
            n_wavelengths
        )


        wavelengths_nm = (
            wavelengths_um
            * 1000.0
        )


        wavelengths_A = (
            wavelengths_nm
            * 10.0
        )


        ######################################################
        # Select nearest lower original Zemax configuration
        ######################################################

        base_config = (
            (order - 60) // 7
            + 1
        )


        ######################################################
        # Storage for complete fiber image at each wavelength
        ######################################################

        fiber_x_by_wave = [
            []
            for _ in range(
                n_wavelengths
            )
        ]


        fiber_y_by_wave = [
            []
            for _ in range(
                n_wavelengths
            )
        ]


        ######################################################
        # Loop over fiber points
        #
        # One trace_all_order call traces ALL wavelengths
        # for one fiber position.
        ######################################################

        for fiber_index, (hx, hy) in enumerate(
            zip(
                hx_fiber,
                hy_fiber
            ),
            start=1
        ):

            result = emar_utils.trace_all_order(
                ln=ln,
                base_config=base_config,
                target_order=order,
                wavelengths=wavelengths_um,
                surf=surf_detector,
                px=px,
                py=py,
                row_order=row_order,
                row_wave1=row_wave1,
                hx=float(hx),
                hy=float(hy)
            )


            ##################################################
            # Retrieve the 100 wavelength footprints
            ##################################################

            spots = list(
                result[
                    "spots"
                ].values()
            )


            if len(spots) != n_wavelengths:

                raise RuntimeError(
                    f"Order m={order}, "
                    f"fiber point {fiber_index}: "
                    f"{len(spots)} spots returned for "
                    f"{n_wavelengths} wavelengths."
                )


            ##################################################
            # Accumulate rays wavelength by wavelength
            ##################################################

            for wave_index, spot in enumerate(
                spots
            ):

                x = np.asarray(
                    spot["x"],
                    dtype=float
                )


                y = np.asarray(
                    spot["y"],
                    dtype=float
                )


                if len(x) > 0:

                    fiber_x_by_wave[
                        wave_index
                    ].extend(
                        x
                    )


                    fiber_y_by_wave[
                        wave_index
                    ].extend(
                        y
                    )


            ##################################################
            # Global progress
            ##################################################

            completed_calls += 1


            print_progress(
                current=completed_calls,
                total=total_calls,
                order=order,
                fiber_point=fiber_index,
                fiber_total=n_fiber_points
            )


        ######################################################
        # Measure fiber image at every wavelength
        ######################################################

        x_centroid = np.full(
            n_wavelengths,
            np.nan
        )


        y_centroid = np.full(
            n_wavelengths,
            np.nan
        )


        x_minimum = np.full(
            n_wavelengths,
            np.nan
        )


        x_maximum = np.full(
            n_wavelengths,
            np.nan
        )


        y_minimum = np.full(
            n_wavelengths,
            np.nan
        )


        y_maximum = np.full(
            n_wavelengths,
            np.nan
        )


        delta_x_um = np.full(
            n_wavelengths,
            np.nan
        )


        delta_y_um = np.full(
            n_wavelengths,
            np.nan
        )


        valid_ray_count = np.zeros(
            n_wavelengths,
            dtype=int
        )


        ######################################################
        # Wavelength loop
        ######################################################

        for wave_index in range(
            n_wavelengths
        ):

            x = np.asarray(
                fiber_x_by_wave[
                    wave_index
                ],
                dtype=float
            )


            y = np.asarray(
                fiber_y_by_wave[
                    wave_index
                ],
                dtype=float
            )


            valid_ray_count[
                wave_index
            ] = len(
                x
            )


            if len(x) == 0:

                continue


            ##################################################
            # Fiber centroid
            ##################################################

            x_centroid[
                wave_index
            ] = np.mean(
                x
            )


            y_centroid[
                wave_index
            ] = np.mean(
                y
            )


            ##################################################
            # Fiber image envelope
            ##################################################

            x_minimum[
                wave_index
            ] = np.min(
                x
            )


            x_maximum[
                wave_index
            ] = np.max(
                x
            )


            y_minimum[
                wave_index
            ] = np.min(
                y
            )


            y_maximum[
                wave_index
            ] = np.max(
                y
            )


            ##################################################
            # Fiber image size
            ##################################################

            delta_x_um[
                wave_index
            ] = (
                x_maximum[
                    wave_index
                ]
                -
                x_minimum[
                    wave_index
                ]
            ) * 1000.0


            delta_y_um[
                wave_index
            ] = (
                y_maximum[
                    wave_index
                ]
                -
                y_minimum[
                    wave_index
                ]
            ) * 1000.0


        ######################################################
        # Check valid spectral points
        ######################################################

        valid = (
            np.isfinite(
                x_centroid
            )
            &
            np.isfinite(
                wavelengths_nm
            )
        )


        if np.count_nonzero(
            valid
        ) < 3:

            raise RuntimeError(
                f"Order m={order}: insufficient valid "
                f"spectral points."
            )


        ######################################################
        # Reciprocal dispersion along detector X
        #
        # dx/dlambda:
        #
        #       mm / nm
        ######################################################

        dx_dlambda = np.full(
            n_wavelengths,
            np.nan
        )


        dx_dlambda[
            valid
        ] = np.gradient(
            x_centroid[
                valid
            ],
            wavelengths_nm[
                valid
            ]
        )


        ######################################################
        # d(lambda)/dX
        #
        #       nm / mm
        ######################################################

        reciprocal_dispersion = np.full(
            n_wavelengths,
            np.nan
        )


        reciprocal_dispersion[
            valid
        ] = (
            1.0
            / np.abs(
                dx_dlambda[
                    valid
                ]
            )
        )


        ######################################################
        # Linear dispersion
        #
        #       Angstrom / mm
        ######################################################

        linear_dispersion_A_per_mm = (
            reciprocal_dispersion
            * 10.0
        )


        ######################################################
        # Our definition:
        #
        #       omega' = Delta X
        ######################################################

        omega_prime_um = (
            delta_x_um.copy()
        )


        omega_prime_mm = (
            omega_prime_um
            / 1000.0
        )


        ######################################################
        # Spectral purity
        #
        #       nm
        ######################################################

        delta_lambda_nm = (
            omega_prime_mm
            * reciprocal_dispersion
        )


        ######################################################
        # Spectral purity
        #
        #       pm
        ######################################################

        delta_lambda_pm = (
            delta_lambda_nm
            * 1000.0
        )


        ######################################################
        # Resolving power
        ######################################################

        resolving_power = (
            wavelengths_nm
            / delta_lambda_nm
        )


        ######################################################
        # Build DataFrame for this order
        ######################################################

        order_df = pd.DataFrame(
            {
                "order":
                    np.full(
                        n_wavelengths,
                        order,
                        dtype=int
                    ),

                "status":
                    np.full(
                        n_wavelengths,
                        status,
                        dtype=object
                    ),

                "wavelength_um":
                    wavelengths_um,

                "wavelength_nm":
                    wavelengths_nm,

                "wavelength_A":
                    wavelengths_A,

                "x_centroid_mm":
                    x_centroid,

                "y_centroid_mm":
                    y_centroid,

                "x_min_mm":
                    x_minimum,

                "x_max_mm":
                    x_maximum,

                "y_min_mm":
                    y_minimum,

                "y_max_mm":
                    y_maximum,

                "delta_x_um":
                    delta_x_um,

                "delta_y_um":
                    delta_y_um,

                "omega_prime_um":
                    omega_prime_um,

                "reciprocal_dispersion_nm_per_mm":
                    reciprocal_dispersion,

                "linear_dispersion_A_per_mm":
                    linear_dispersion_A_per_mm,

                "delta_lambda_nm":
                    delta_lambda_nm,

                "delta_lambda_pm":
                    delta_lambda_pm,

                "resolving_power":
                    resolving_power,

                "valid_ray_count":
                    valid_ray_count
            }
        )


        ######################################################
        # Store in memory
        ######################################################

        all_order_results.append(
            order_df
        )


        ######################################################
        # Save immediately after each completed order
        #
        # This protects the long calculation if Zemax/PyZDDE
        # stops before all 85 orders are finished.
        ######################################################

        current_results = pd.concat(
            all_order_results,
            ignore_index=True
        )


        current_results.to_csv(
            output_file,
            index=False
        )


finally:

    ##########################################################
    # Close Zemax
    ##########################################################

    ln.close()


##############################################################
# Combine complete results
##############################################################

results_df = pd.concat(
    all_order_results,
    ignore_index=True
)


##############################################################
# Final save
##############################################################

results_df.to_csv(
    output_file,
    index=False
)


##############################################################
# Color convention
#
# Low order  -> red
# High order -> blue
##############################################################

cmap = plt.cm.turbo_r


norm = Normalize(
    vmin=orders.min(),
    vmax=orders.max()
)


##############################################################
# Colorbar helper
##############################################################

def add_order_colorbar(
    fig,
    ax
):

    sm = plt.cm.ScalarMappable(
        norm=norm,
        cmap=cmap
    )


    sm.set_array([])


    cbar = fig.colorbar(
        sm,
        ax=ax,
        pad=0.02
    )


    cbar.set_label(
        r"Diffraction order, $m$",
        fontsize=18
    )


    cbar.ax.tick_params(
        labelsize=14
    )


##############################################################
# Common axis style
##############################################################

def style_axis(
    ax
):

    ax.minorticks_on()


    ax.tick_params(
        which="major",
        direction="in",
        top=True,
        right=True,
        length=6,
        width=2.2,
        labelsize=16
    )


    ax.tick_params(
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=3,
        width=1.2
    )


##############################################################
# Plot 1
#
# Linear dispersion for all orders
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


for order in orders:

    data = results_df[
        results_df[
            "order"
        ] == order
    ].sort_values(
        "wavelength_A"
    )


    ax.plot(
        data[
            "wavelength_A"
        ],
        data[
            "linear_dispersion_A_per_mm"
        ],
        color=cmap(
            norm(
                order
            )
        ),
        linewidth=1.4
    )


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\AA)$",
    fontsize=18
)


ax.set_ylabel(
    r"$\mathrm{Linear\ dispersion}\;(\AA/\mathrm{mm})$",
    fontsize=18
)


style_axis(
    ax
)


add_order_colorbar(
    fig,
    ax
)


fig.tight_layout()

plt.show()


##############################################################
# Plot 2
#
# Complete resolving power
#
# omega'(lambda) is independently measured at all
# 100 wavelengths of every order.
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


for order in orders:

    data = results_df[
        results_df[
            "order"
        ] == order
    ].sort_values(
        "wavelength_A"
    )


    ax.plot(
        data[
            "wavelength_A"
        ],
        data[
            "resolving_power"
        ],
        color=cmap(
            norm(
                order
            )
        ),
        linewidth=1.3
    )


ax.set_xlabel(
    r"$\mathrm{Wavelength}\;(\AA)$",
    fontsize=18
)


ax.set_ylabel(
    r"$\mathrm{Resolving\ power},\ R=\lambda/\Delta\lambda$",
    fontsize=18
)


style_axis(
    ax
)


add_order_colorbar(
    fig,
    ax
)


fig.tight_layout()

plt.show()


##############################################################
# Plot 3
#
# Echelle spectral format:
# fiber-image centroids
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 7)
)


for order in orders:

    data = results_df[
        results_df[
            "order"
        ] == order
    ].sort_values(
        "wavelength_nm"
    )


    ax.plot(
        data[
            "x_centroid_mm"
        ],
        data[
            "y_centroid_mm"
        ],
        color=cmap(
            norm(
                order
            )
        ),
        linewidth=1.2
    )


ax.set_xlabel(
    r"$X\;(\mathrm{mm})$",
    fontsize=18
)


ax.set_ylabel(
    r"$Y\;(\mathrm{mm})$",
    fontsize=18
)


ax.set_aspect(
    "equal",
    adjustable="box"
)


style_axis(
    ax
)


add_order_colorbar(
    fig,
    ax
)


fig.tight_layout()

plt.show()


##############################################################
# Plot 4
#
# Echelle spectral format including projected fiber size
#
# The envelope is NOT interpolated.
#
# ymin(lambda) and ymax(lambda) come directly from the
# complete fiber ray tracing at every wavelength.
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 7)
)


for order in orders:

    data = results_df[
        results_df[
            "order"
        ] == order
    ].sort_values(
        "wavelength_nm"
    )


    color = cmap(
        norm(
            order
        )
    )


    ##########################################################
    # Coordinates
    ##########################################################

    xc = data[
        "x_centroid_mm"
    ].to_numpy()


    yc = data[
        "y_centroid_mm"
    ].to_numpy()


    ymin = data[
        "y_min_mm"
    ].to_numpy()


    ymax = data[
        "y_max_mm"
    ].to_numpy()


    ##########################################################
    # Fiber envelope
    ##########################################################

    ax.fill_between(
        xc,
        ymin,
        ymax,
        color=color,
        alpha=0.55,
        linewidth=0.0
    )


    ##########################################################
    # Spectral centroid
    ##########################################################

    ax.plot(
        xc,
        yc,
        color=color,
        linewidth=0.8
    )


ax.set_xlabel(
    r"$X\;(\mathrm{mm})$",
    fontsize=18
)


ax.set_ylabel(
    r"$Y\;(\mathrm{mm})$",
    fontsize=18
)


ax.set_aspect(
    "equal",
    adjustable="box"
)


style_axis(
    ax
)


add_order_colorbar(
    fig,
    ax
)


fig.tight_layout()

plt.show()


##############################################################
# Summary
##############################################################

print()

print(
    "Full EMAR fiber analysis completed"
)

print(
    "----------------------------------"
)


print(
    f"Orders:                    "
    f"{orders.min()} - {orders.max()}"
)


print(
    f"Total spectral points:     "
    f"{len(results_df):,}"
)


print()

print(
    "Fiber image"
)

print(
    "-----------"
)


print(
    f"Delta X range:             "
    f"{results_df['delta_x_um'].min():.3f} - "
    f"{results_df['delta_x_um'].max():.3f} um"
)


print(
    f"Mean Delta X:              "
    f"{results_df['delta_x_um'].mean():.3f} um"
)


print(
    f"Delta Y range:             "
    f"{results_df['delta_y_um'].min():.3f} - "
    f"{results_df['delta_y_um'].max():.3f} um"
)


print(
    f"Mean Delta Y:              "
    f"{results_df['delta_y_um'].mean():.3f} um"
)


print()

print(
    "Linear dispersion"
)

print(
    "-----------------"
)


print(
    f"Range:                     "
    f"{results_df['linear_dispersion_A_per_mm'].min():.3f} - "
    f"{results_df['linear_dispersion_A_per_mm'].max():.3f} A/mm"
)


print(
    f"Mean:                      "
    f"{results_df['linear_dispersion_A_per_mm'].mean():.3f} A/mm"
)


print()

print(
    "Spectral resolution"
)

print(
    "-------------------"
)


print(
    f"Delta lambda range:        "
    f"{results_df['delta_lambda_pm'].min():.3f} - "
    f"{results_df['delta_lambda_pm'].max():.3f} pm"
)


print(
    f"R minimum:                 "
    f"{results_df['resolving_power'].min():.0f}"
)


print(
    f"R maximum:                 "
    f"{results_df['resolving_power'].max():.0f}"
)


print(
    f"R mean:                    "
    f"{results_df['resolving_power'].mean():.0f}"
)


print()

print(
    "Results saved to:"
)

print(
    output_file
)