# -*- coding: utf-8 -*-

"""
15_order_separation.py
======================

Geometrical separation between adjacent EMAR echelle orders using
the projected 100 um fiber envelope obtained with the robust
91-field fiber sampling.

No Zemax ray tracing is performed.

For every pair of adjacent diffraction orders:

    (60, 61), (61, 62), ..., (143, 144)

the fiber envelopes are compared at the same detector X coordinate.

For an upper and lower order, the free geometrical gap is defined as:

    gap(X) = Ymin_upper(X) - Ymax_lower(X)

where:

    Ymin_upper = lower edge of the upper spectral order
    Ymax_lower = upper edge of the lower spectral order

The two envelopes are interpolated over their common X interval.

For every adjacent-order pair the script calculates:

    minimum gap
    maximum gap
    mean gap
    median gap
    standard deviation
    X coordinate of minimum gap
    X coordinate of maximum gap

Positive gap:
    free space between projected fiber envelopes.

Zero gap:
    envelopes touch.

Negative gap:
    geometrical overlap.

Outputs
-------
Results/order_separation_91fields.csv
Results/order_separation_statistics.png
Results/order_minimum_separation.png

Author
------
Morgan Rhaí Nájera Roa
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
# Input file
##############################################################

input_file = (
    emar_utils.RESULTS_DIR
    / "full_orders_fiber_resolution.csv"
)


##############################################################
# Output files
##############################################################

output_table = (
    emar_utils.RESULTS_DIR
    / "order_separation_91fields.csv"
)


output_statistics_plot = (
    emar_utils.RESULTS_DIR
    / "order_separation_statistics.png"
)


output_minimum_plot = (
    emar_utils.RESULTS_DIR
    / "order_minimum_separation.png"
)

output_gap_profiles_plot = (
    emar_utils.RESULTS_DIR
    / "order_separation_vs_detector_x.png"
)

##############################################################
# Interpolation sampling
##############################################################

n_x_samples = 2000


##############################################################
# Load robust 91-field results
##############################################################

df = pd.read_csv(
    input_file
)


##############################################################
# Required columns
##############################################################

required_columns = [
    "order",
    "x_centroid_mm",
    "y_centroid_mm",
    "y_min_mm",
    "y_max_mm"
]


missing_columns = [
    column
    for column in required_columns
    if column not in df.columns
]


if missing_columns:

    raise ValueError(
        "Missing required columns: "
        + ", ".join(
            missing_columns
        )
    )


##############################################################
# Available diffraction orders
##############################################################

orders = np.sort(
    df[
        "order"
    ].unique()
)


print()

print(
    "EMAR adjacent-order separation analysis"
)

print(
    "---------------------------------------"
)

print(
    f"Orders: "
    f"{orders.min()} - {orders.max()}"
)

print(
    f"Number of orders: "
    f"{len(orders)}"
)

print(
    f"Adjacent pairs: "
    f"{len(orders) - 1}"
)

print()


##############################################################
# Storage
##############################################################

results = []

gap_profiles = {}

##############################################################
# Analyze adjacent orders
##############################################################

for i in range(
    len(orders) - 1
):

    order_a = int(
        orders[i]
    )

    order_b = int(
        orders[i + 1]
    )


    ##########################################################
    # Extract both orders
    ##########################################################

    data_a = df[
        df[
            "order"
        ] == order_a
    ].copy()


    data_b = df[
        df[
            "order"
        ] == order_b
    ].copy()


    ##########################################################
    # Sort by detector X
    ##########################################################

    data_a = data_a.sort_values(
        "x_centroid_mm"
    )


    data_b = data_b.sort_values(
        "x_centroid_mm"
    )


    ##########################################################
    # Determine which order is physically upper in Y
    ##########################################################

    mean_y_a = data_a[
        "y_centroid_mm"
    ].mean()


    mean_y_b = data_b[
        "y_centroid_mm"
    ].mean()


    if mean_y_a < mean_y_b:

        lower_order = order_a
        upper_order = order_b

        lower_data = data_a
        upper_data = data_b

    else:

        lower_order = order_b
        upper_order = order_a

        lower_data = data_b
        upper_data = data_a


    ##########################################################
    # Common detector-X interval
    ##########################################################

    x_low = max(
        lower_data[
            "x_centroid_mm"
        ].min(),

        upper_data[
            "x_centroid_mm"
        ].min()
    )


    x_high = min(
        lower_data[
            "x_centroid_mm"
        ].max(),

        upper_data[
            "x_centroid_mm"
        ].max()
    )


    ##########################################################
    # Check overlap in X
    ##########################################################

    if x_high <= x_low:

        print(
            f"Orders {order_a}-{order_b}: "
            "no common X interval"
        )

        continue


    ##########################################################
    # Common X grid
    ##########################################################

    x_common = np.linspace(
        x_low,
        x_high,
        n_x_samples
    )


    ##########################################################
    # Interpolate upper edge of lower order
    ##########################################################

    lower_ymax = np.interp(
        x_common,

        lower_data[
            "x_centroid_mm"
        ],

        lower_data[
            "y_max_mm"
        ]
    )


    ##########################################################
    # Interpolate lower edge of upper order
    ##########################################################

    upper_ymin = np.interp(
        x_common,

        upper_data[
            "x_centroid_mm"
        ],

        upper_data[
            "y_min_mm"
        ]
    )


    ##########################################################
    # Free geometrical gap
    ##########################################################

    gap_mm = (
        upper_ymin
        - lower_ymax
    )


    gap_um = (
        gap_mm
        * 1000.0
    )
    
    gap_profiles[
    (order_a, order_b)
    ] = {
        "x_mm": x_common.copy(),
        "gap_um": gap_um.copy()
    }
        
    ##########################################################
    # Fiber-envelope separation at detector center
    ##########################################################
    
    if x_low <= 0.0 <= x_high:
    
        gap_at_x0_um = float(
            np.interp(
                0.0,
                x_common,
                gap_um
            )
        )
    
    else:
    
        gap_at_x0_um = np.nan    


    ##########################################################
    # Statistics
    ##########################################################

    min_index = np.argmin(
        gap_um
    )


    max_index = np.argmax(
        gap_um
    )


    min_gap_um = gap_um[
        min_index
    ]


    max_gap_um = gap_um[
        max_index
    ]


    mean_gap_um = np.mean(
        gap_um
    )


    median_gap_um = np.median(
        gap_um
    )


    std_gap_um = np.std(
        gap_um
    )


    x_at_min_gap_mm = x_common[
        min_index
    ]


    x_at_max_gap_mm = x_common[
        max_index
    ]


    ##########################################################
    # Save pair result
    ##########################################################

    results.append(
        {
            "order_a":
                order_a,

            "order_b":
                order_b,

            "lower_order":
                lower_order,

            "upper_order":
                upper_order,

            "x_common_min_mm":
                x_low,

            "x_common_max_mm":
                x_high,

            "min_gap_um":
                min_gap_um,

            "max_gap_um":
                max_gap_um,

            "mean_gap_um":
                mean_gap_um,

            "median_gap_um":
                median_gap_um,

            "std_gap_um":
                std_gap_um,

            "x_at_min_gap_mm":
                x_at_min_gap_mm,

            "x_at_max_gap_mm":
                x_at_max_gap_mm,
                
            "gap_at_x0_um":
                gap_at_x0_um    
        }
    )


##############################################################
# Results dataframe
##############################################################

separation = pd.DataFrame(
    results
)


##############################################################
# Save table
##############################################################

separation.to_csv(
    output_table,
    index=False
)


##############################################################
# Global minimum
##############################################################

global_min_index = (
    separation[
        "min_gap_um"
    ].idxmin()
)


global_min = separation.loc[
    global_min_index
]


##############################################################
# Global maximum
##############################################################

global_max_index = (
    separation[
        "max_gap_um"
    ].idxmax()
)


global_max = separation.loc[
    global_max_index
]


##############################################################
# Console summary
##############################################################

print()

print(
    "Global separation statistics"
)

print(
    "----------------------------"
)


print(
    f"Mean of minimum gaps: "
    f"{separation['min_gap_um'].mean():.3f} um"
)


print(
    f"Median of minimum gaps: "
    f"{separation['min_gap_um'].median():.3f} um"
)


print(
    f"Mean of maximum gaps: "
    f"{separation['max_gap_um'].mean():.3f} um"
)


print(
    f"Median of maximum gaps: "
    f"{separation['max_gap_um'].median():.3f} um"
)


print()


print(
    "Minimum separation in complete format"
)

print(
    "-------------------------------------"
)


print(
    f"Orders: "
    f"{int(global_min['order_a'])}-"
    f"{int(global_min['order_b'])}"
)


print(
    f"Gap: "
    f"{global_min['min_gap_um']:.3f} um"
)


print(
    f"Detector X: "
    f"{global_min['x_at_min_gap_mm']:.6f} mm"
)


print()


print(
    "Maximum separation in complete format"
)

print(
    "-------------------------------------"
)


print(
    f"Orders: "
    f"{int(global_max['order_a'])}-"
    f"{int(global_max['order_b'])}"
)


print(
    f"Gap: "
    f"{global_max['max_gap_um']:.3f} um"
)


print(
    f"Detector X: "
    f"{global_max['x_at_max_gap_mm']:.6f} mm"
)


print()

print(
    f"Mean gap at CCD center: "
    f"{separation['gap_at_x0_um'].mean():.3f} um"
)


print(
    f"Median gap at CCD center: "
    f"{separation['gap_at_x0_um'].median():.3f} um"
)


print()

##############################################################
# Check for geometrical overlap
##############################################################

overlap_pairs = separation[
    separation[
        "min_gap_um"
    ] < 0.0
]


if len(
    overlap_pairs
) == 0:

    print(
        "No geometrical overlap between adjacent "
        "fiber envelopes."
    )

else:

    print(
        f"WARNING: {len(overlap_pairs)} adjacent "
        "order pairs show geometrical overlap."
    )


print()


##############################################################
# Minimum, CCD-center and maximum separation
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


x_orders = separation[
    "order_a"
]


##############################################################
# Minimum gap
##############################################################

ax.plot(
    x_orders,
    separation[
        "min_gap_um"
    ],
    marker="o",
    markersize=3,
    linewidth=1.5,
    label="Minimum gap"
)


##############################################################
# Gap at detector center
##############################################################

ax.plot(
    x_orders,
    separation[
        "gap_at_x0_um"
    ],
    marker="o",
    markersize=3,
    linewidth=1.5,
    label=r"Gap at CCD center ($X=0$)"
)


##############################################################
# Maximum gap
##############################################################

ax.plot(
    x_orders,
    separation[
        "max_gap_um"
    ],
    marker="o",
    markersize=3,
    linewidth=1.5,
    label="Maximum gap"
)


##############################################################
# Zero-separation reference
##############################################################

ax.axhline(
    0.0,
    color="black",
    linestyle="--",
    linewidth=1.0
)


##############################################################
# Labels
##############################################################

ax.set_xlabel(
    "Lower diffraction order m",
    fontsize=18
)


ax.set_ylabel(
    r"Fiber-envelope separation [$\mu$m]",
    fontsize=18
)


ax.set_title(
    "EMAR Adjacent-Order Separation",
    fontsize=18
)


##############################################################
# Legend
##############################################################

ax.legend(
    fontsize=13
)


##############################################################
# Tick style
##############################################################

ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5
)


ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0
)


ax.minorticks_on()


##############################################################
# Save
##############################################################

fig.tight_layout()


fig.savefig(
    output_statistics_plot,
    dpi=300,
    bbox_inches="tight"
)


##############################################################
# Minimum separation only
##############################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)


ax.plot(
    separation[
        "order_a"
    ],
    separation[
        "min_gap_um"
    ],
    marker="o",
    markersize=4,
    linewidth=1.5
)


ax.axhline(
    0.0,
    color="black",
    linestyle="--",
    linewidth=1.0
)


##############################################################
# Mark global minimum
##############################################################

ax.scatter(
    global_min[
        "order_a"
    ],
    global_min[
        "min_gap_um"
    ],
    s=70,
    zorder=5
)


ax.annotate(
    (
        f"m={int(global_min['order_a'])}-"
        f"{int(global_min['order_b'])}\n"
        f"{global_min['min_gap_um']:.1f} "
        r"$\mu$m"
    ),

    xy=(
        global_min[
            "order_a"
        ],

        global_min[
            "min_gap_um"
        ]
    ),

    xytext=(
        20,
        20
    ),

    textcoords="offset points",

    fontsize=13,

    arrowprops={
        "arrowstyle": "->"
    }
)


ax.set_xlabel(
    "Lower diffraction order m",
    fontsize=18
)


ax.set_ylabel(
    r"Minimum fiber-envelope separation [$\mu$m]",
    fontsize=18
)


ax.set_title(
    "EMAR Minimum Adjacent-Order Separation",
    fontsize=18
)


ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5
)


ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0
)


ax.minorticks_on()


fig.tight_layout()


fig.savefig(
    output_minimum_plot,
    dpi=300,
    bbox_inches="tight"
)

##############################################################
# Fiber-envelope separation as a function of detector X
##############################################################

representative_pairs = [
    (60, 61),
    (102, 103),
    (143, 144)
]


fig, ax = plt.subplots(
    figsize=(10, 6)
)


for pair in representative_pairs:

    if pair not in gap_profiles:

        print(
            f"WARNING: separation profile "
            f"{pair[0]}-{pair[1]} not available."
        )

        continue


    ##########################################################
    # Retrieve profile
    ##########################################################

    x_profile = gap_profiles[
        pair
    ][
        "x_mm"
    ]


    gap_profile = gap_profiles[
        pair
    ][
        "gap_um"
    ]


    ##########################################################
    # Minimum of this pair
    ##########################################################

    i_min = np.argmin(
        gap_profile
    )


    x_min = x_profile[
        i_min
    ]


    gap_min = gap_profile[
        i_min
    ]


    ##########################################################
    # Plot separation profile
    ##########################################################

    line, = ax.plot(
        x_profile,
        gap_profile,
        linewidth=2.0,
        label=(
            f"m={pair[0]}-{pair[1]}"
        )
    )


    ##########################################################
    # Mark minimum using same line color
    ##########################################################

    ax.scatter(
        x_min,
        gap_min,
        s=55,
        color=line.get_color(),
        zorder=5
    )


##############################################################
# Detector center
##############################################################

ax.axvline(
    0.0,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label="CCD center"
)


##############################################################
# Zero separation reference
##############################################################

ax.axhline(
    0.0,
    color="black",
    linestyle=":",
    linewidth=1.0
)


##############################################################
# Labels
##############################################################

ax.set_xlabel(
    "Detector X [mm]",
    fontsize=18
)


ax.set_ylabel(
    r"Fiber-envelope separation [$\mu$m]",
    fontsize=18
)


ax.set_title(
    "EMAR Adjacent-Order Separation Across the Detector",
    fontsize=18
)


##############################################################
# Legend
##############################################################

ax.legend(
    fontsize=13
)


##############################################################
# Tick style
##############################################################

ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    top=True,
    right=True,
    labelsize=16,
    length=6,
    width=1.5
)


ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    top=True,
    right=True,
    length=3,
    width=1.0
)


ax.minorticks_on()


##############################################################
# Save
##############################################################

fig.tight_layout()


fig.savefig(
    output_gap_profiles_plot,
    dpi=300,
    bbox_inches="tight"
)


##############################################################
# Show plots
##############################################################

plt.show()
