import sys
from pathlib import Path

# ============================================================
# Project root
# ============================================================

project_root = Path(__file__).resolve().parents[1]

if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# ============================================================
# Imports
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyzdde.zdde as pyz

from scipy.spatial import ConvexHull
from shapely.geometry import Polygon

from utils.emar_plots import apply_scientific_style
from utils import emar_utils

# ============================================================
# Input focal-plane calibration
# ============================================================

fiber_diameter_um = 100.0
fiber_radius_um = fiber_diameter_um / 2.0

hx_max = 0.286480
hy_max = 0.285400

# Center-to-center separation between polarimetric images
d_um = 100.0

center_plus_y_um = +d_um / 2.0
center_minus_y_um = -d_um / 2.0


# ============================================================
# Nine representative fiber fields
# ============================================================

angles_deg = np.arange(0.0, 360.0, 45.0)
angles_rad = np.deg2rad(angles_deg)

labels = ["center"] + [
    f"{int(angle)} deg" for angle in angles_deg
]

image_x_um = np.concatenate((
    [0.0],
    fiber_radius_um * np.cos(angles_rad)
))

image_y_offset_um = np.concatenate((
    [0.0],
    fiber_radius_um * np.sin(angles_rad)
))

image_plus_y_um = (
    center_plus_y_um + image_y_offset_um
)

image_minus_y_um = (
    center_minus_y_um + image_y_offset_um
)

# ============================================================
# Physical coordinates -> Zemax field coordinates
# ============================================================

image_hx = (image_x_um / fiber_radius_um) * hx_max

image_plus_hy = (image_plus_y_um / fiber_radius_um) * hy_max

image_minus_hy = (image_minus_y_um / fiber_radius_um) * hy_max


# ============================================================
# Hexapolar pupil sampling
# ============================================================

pupil_px, pupil_py = emar_utils.sample_hexapolar_pupil(
    n_rings=5,
    include_center=True,
)

n_pupil_rays = len(pupil_px)
n_fields = 2 * len(image_hx)

# ============================================================
# Input focal-plane visualization
# ============================================================

fig, ax = plt.subplots(figsize=(7, 7))

circle_plus = plt.Circle(
    (0.0, center_plus_y_um),
    fiber_radius_um,
    fill=False,
    linewidth=2,
    color="k",
)

circle_minus = plt.Circle(
    (0.0, center_minus_y_um),
    fiber_radius_um,
    fill=False,
    linewidth=2,
    color="k",
)

ax.add_patch(circle_plus)
ax.add_patch(circle_minus)

ax.scatter(
    image_x_um,
    image_plus_y_um,
    marker="x",
    linewidths=1.5,
    s=50,
)

ax.scatter(
    image_x_um,
    image_minus_y_um,
    marker="x",
    linewidths=1.5,
    s=50,
)

ax.set_aspect("equal", adjustable="box")

apply_scientific_style(
    ax,
    xlabel=r"X [$\mu$m]",
    ylabel=r"Y [$\mu$m]",
)

margin_um = 25.0

ax.set_xlim(
    -fiber_radius_um - margin_um,
    +fiber_radius_um + margin_um,
)

ax.set_ylim(
    center_minus_y_um - fiber_radius_um - margin_um,
    center_plus_y_um + fiber_radius_um + margin_um,
)

plt.tight_layout()
plt.show()
 
    
# ============================================================
# Zemax reference configuration
# ============================================================

order = 60
base_config = ((order - 60) // 7) + 1

surf = 54
wave_num = 6    
    
# ============================================================
# Zemax connection
# ============================================================

zemax_file = project_root / "zemax" / "EMAR.zmx"

ln = pyz.createLink()

if ln is None:
    raise RuntimeError(
        "Could not establish a PyZDDE connection to Zemax OpticStudio."
    )

ln.zLoadFile(str(zemax_file))

# ============================================================
# Activate reference configuration
# ============================================================

ln.zSetConfig(base_config)
ln.zGetUpdate()

wavelength_um = ln.zGetWave(wave_num).wavelength

print("\nZemax reference configuration")
print("-----------------------------")
print(f"Diffraction order : {order}")
print(f"Configuration     : {base_config}")
print(f"Wavelength slot   : {wave_num}")
print(f"Wavelength        : {wavelength_um:.9f} um")
print(f"Detector surface  : {surf}")

# ============================================================
# Single-ray field propagation test
# ============================================================

px = 0.0
py = 0.0

results = []

for image_name, hy_values in [
    ("Upper", image_plus_hy),
    ("Lower", image_minus_hy),
]:
    for label, hx, hy in zip(labels, image_hx, hy_values):

        ray = ln.zGetTrace(
            waveNum=wave_num,
            mode=0,
            surf=surf,
            hx=hx,
            hy=hy,
            px=px,
            py=py,
        )

        results.append({
            "Image": image_name,
            "Field": label,
            "Hx": hx,
            "Hy": hy,
            "X [mm]": ray[2],
            "Y [mm]": ray[3],
            "Error": ray[0],
            "Vignette": ray[1]
        })


df = pd.DataFrame(results)

# ============================================================
# Central-ray tracing summary
# ============================================================

print("\nCentral-ray validation")
print(f"Fields traced : {len(df)}")
print(f"Errors        : {(df['Error'] != 0).sum()}")
print(f"Vignetted     : {(df['Vignette'] != 0).sum()}")

# ============================================================
# Hexapolar pupil ray tracing
# ============================================================

ray_rows = []

for image_name, hy_values in [
    ("Upper", image_plus_hy),
    ("Lower", image_minus_hy),
]:

    for field_name, hx, hy in zip(
        labels, image_hx, hy_values
    ):

        for pupil_index, (pxi, pyi) in enumerate(
            zip(pupil_px, pupil_py)
        ):

            ray = ln.zGetTrace(
                waveNum=wave_num,
                mode=0,
                surf=surf,
                hx=float(hx),
                hy=float(hy),
                px=float(pxi),
                py=float(pyi),
            )

            ray_rows.append({
                "Image": image_name,
                "Field": field_name,
                "Pupil": pupil_index,
                "Hx": hx,
                "Hy": hy,
                "Px": pxi,
                "Py": pyi,
                "X [mm]": ray[2],
                "Y [mm]": ray[3],
                "Error": ray[0],
                "Vignette": ray[1],
            })


rays_df = pd.DataFrame(ray_rows)

# ============================================================
# Pupil ray-tracing validation
# ============================================================

rays_df["Valid"] = (
    (rays_df["Error"] == 0) &
    (rays_df["Vignette"] == 0) &
    np.isfinite(rays_df["X [mm]"]) &
    np.isfinite(rays_df["Y [mm]"])
)

validation_df = (
    rays_df.groupby(
        ["Image", "Field"],
        sort=False
    )
    .agg(
        Traced=("Pupil", "count"),
        Valid=("Valid", "sum"),
        Errors=("Error", lambda x: (x != 0).sum()),
        Vignetted=("Vignette", lambda x: (x != 0).sum()),
    )
    .reset_index()
)

print("\nHexapolar pupil validation")
print(validation_df.to_string(index=False))

print(f"\nTotal rays traced: {len(rays_df)}")
print(f"Total valid rays: {rays_df['Valid'].sum()}")
assert len(rays_df) == 2 * len(labels) * n_pupil_rays
assert len(labels) == 9 and n_pupil_rays == 91

# ============================================================
# Detector-plane separation validation
# ============================================================

# Image centers
upper_center = df[
    (df["Image"] == "Upper") &
    (df["Field"] == "center")
].iloc[0]

lower_center = df[
    (df["Image"] == "Lower") &
    (df["Field"] == "center")
].iloc[0]

delta_x_center_um = (
    upper_center["X [mm]"] -
    lower_center["X [mm]"]
) * 1000.0

delta_y_center_um = (
    upper_center["Y [mm]"] -
    lower_center["Y [mm]"]
) * 1000.0


# Contact points: Upper -Y and Lower +Y
upper_contact = df[
    (df["Image"] == "Upper") &
    (df["Field"] == "270 deg")
].iloc[0]

lower_contact = df[
    (df["Image"] == "Lower") &
    (df["Field"] == "90 deg")
].iloc[0]

delta_x_contact_um = (
    upper_contact["X [mm]"] -
    lower_contact["X [mm]"]
) * 1000.0

delta_y_contact_um = (
    upper_contact["Y [mm]"] -
    lower_contact["Y [mm]"]
) * 1000.0


# Outer Y-field points: Upper +Y and Lower -Y
upper_ymax = df[
    (df["Image"] == "Upper") &
    (df["Field"] == "90 deg")
].iloc[0]

lower_ymin = df[
    (df["Image"] == "Lower") &
    (df["Field"] == "270 deg")
].iloc[0]

delta_x_outer_um = (
    upper_ymax["X [mm]"] -
    lower_ymin["X [mm]"]
) * 1000.0

delta_y_outer_um = (
    upper_ymax["Y [mm]"] -
    lower_ymin["Y [mm]"]
) * 1000.0


# ============================================================
# Summary table
# ============================================================

separation_df = pd.DataFrame([
    {
        "Reference": "Centers",
        "Delta X [um]": delta_x_center_um,
        "Delta Y [um]": delta_y_center_um,
    },
    {
        "Reference": "Contact points",
        "Delta X [um]": delta_x_contact_um,
        "Delta Y [um]": delta_y_contact_um,
    },
    {
        "Reference": "Outer Y fields",
        "Delta X [um]": delta_x_outer_um,
        "Delta Y [um]": delta_y_outer_um,
    },
])

print("\nDetector-plane separation validation")
print(separation_df.to_string(index=False))

# The two tangent entrance points must map to the same detector point.
assert np.isclose(delta_x_contact_um, 0.0, atol=1e-5)
assert np.isclose(delta_y_contact_um, 0.0, atol=1e-5)


# ============================================================
# Contact-field pupil consistency
# ============================================================

upper_contact_rays = (
    rays_df[(rays_df["Image"] == "Upper") & (rays_df["Field"] == "270 deg")]
    .sort_values("Pupil")
)
lower_contact_rays = (
    rays_df[(rays_df["Image"] == "Lower") & (rays_df["Field"] == "90 deg")]
    .sort_values("Pupil")
)

contact_residual_um = 1000.0 * (
    upper_contact_rays[["X [mm]", "Y [mm]"]].to_numpy()
    - lower_contact_rays[["X [mm]", "Y [mm]"]].to_numpy()
)
print("\nMaximum contact-ray discrepancy [um]:", np.max(np.abs(contact_residual_um)))

# ============================================================
# Image-plane visualization
# ============================================================

fig, ax = plt.subplots(figsize=(7, 7))

upper_df = df[df["Image"] == "Upper"]
lower_df = df[df["Image"] == "Lower"]

for image_df, color in [
    (upper_df, "tab:blue"),
    (lower_df, "tab:orange"),
]:

    ax.scatter(
        image_df["X [mm]"],
        image_df["Y [mm]"],
        marker="x",
        linewidths=1.5,
        s=50,
        color=color,
        zorder=3,
    )

# ============================================================
# Hexapolar spot diagrams
# ============================================================

valid_rays = rays_df[rays_df["Valid"]]

for image_name, color in [
    ("Upper", "tab:blue"),
    ("Lower", "tab:orange"),
]:

    image_rays = valid_rays[
        valid_rays["Image"] == image_name
    ]

    ax.scatter(
        image_rays["X [mm]"],
        image_rays["Y [mm]"],
        s=3,
        alpha=0.35,
        color=color,
        marker=".",
        zorder=2,
    )
ax.set_aspect("equal", adjustable="box")

apply_scientific_style(
    ax,
    xlabel="X [mm]",
    ylabel="Y [mm]",
)



plt.tight_layout()
plt.show()

# ============================================================
# Polarimetric image convex hulls
# ============================================================

fig, ax = plt.subplots(figsize=(7, 7))

valid_rays = rays_df[rays_df["Valid"]]

for image_name, color in [
    ("Upper", "tab:blue"),
    ("Lower", "tab:orange"),
]:

    image_rays = valid_rays[
        valid_rays["Image"] == image_name
    ]

    points = image_rays[
        ["X [mm]", "Y [mm]"]
    ].to_numpy()

    # Compute convex hull
    hull = ConvexHull(points)

    # Plot ray intersections
    ax.scatter(
        points[:, 0],
        points[:, 1],
        s=3,
        alpha=0.25,
        color=color,
        marker=".",
    )

    # Draw convex hull boundary
    hull_vertices = points[hull.vertices]

    hull_vertices = np.vstack([
        hull_vertices,
        hull_vertices[0],
    ])

    ax.plot(
        hull_vertices[:, 0],
        hull_vertices[:, 1],
        color=color,
        linewidth=1.5,
    )

ax.set_aspect("equal", adjustable="box")

apply_scientific_style(
    ax,
    xlabel="X [mm]",
    ylabel="Y [mm]",
)


plt.tight_layout()
plt.show()

# ============================================================
# Polarimetric image overlap analysis
# ============================================================

polygons = {}

for image_name in ["Upper", "Lower"]:

    image_rays = rays_df[
        (rays_df["Image"] == image_name) &
        (rays_df["Valid"])
    ]

    points = image_rays[
        ["X [mm]", "Y [mm]"]
    ].to_numpy()

    hull = ConvexHull(points)

    polygons[image_name] = Polygon(
        points[hull.vertices]
    )


# Compute geometric intersection
upper_polygon = polygons["Upper"]
lower_polygon = polygons["Lower"]

intersection = upper_polygon.intersection(lower_polygon)

# Convert mm² to um²
conversion = 1e6

upper_area = upper_polygon.area * conversion
lower_area = lower_polygon.area * conversion
overlap_area = intersection.area * conversion

upper_overlap_fraction = (
    overlap_area / upper_area
)

lower_overlap_fraction = (
    overlap_area / lower_area
)


# ============================================================
# Overlap summary
# ============================================================

overlap_df = pd.DataFrame([
    {
        "Image": "Upper",
        "Area [um²]": upper_area,
        "Overlap [um²]": overlap_area,
        "Overlap [%]": 100.0 * upper_overlap_fraction,
    },
    {
        "Image": "Lower",
        "Area [um²]": lower_area,
        "Overlap [um²]": overlap_area,
        "Overlap [%]": 100.0 * lower_overlap_fraction,
    },
])

print("\nPolarimetric image overlap")
print(overlap_df.to_string(index=False))