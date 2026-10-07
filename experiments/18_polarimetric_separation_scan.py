"""Experiment 18: detector-plane overlap versus polarimetric image separation.

Two computational images originate from one physical 100-um circular fiber.
Each image is represented by its center and eight boundary fields. For each
field, 91 hexapolar pupil rays are traced in the selected Zemax configuration.
Convex-hull overlap is a geometric sampling metric, not an intensity leakage
or polarimetric cross-talk measurement.
"""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyzdde.zdde as pyz
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon

from utils import emar_utils
from utils.emar_plots import apply_scientific_style

# ============================================================
# Physical and numerical parameters
# ============================================================

fiber_diameter_um = 100.0
fiber_radius_um = fiber_diameter_um / 2.0
hx_max = 0.286480
hy_max = 0.285400

d_values_um = np.arange(110.0, 115.5, 0.5)

order = 60
base_config = ((order - 60) // 7) + 1
wave_num = 6
surf = 54

zemax_file = (
    emar_utils.ZEMAX_DIR
    / "WP - Con prismas diseñados - camara.zmx"
)

# Nine points: fiber center and eight points on its circumference.
angles_rad = np.deg2rad(np.arange(0.0, 360.0, 45.0))
x_offsets_um = np.concatenate(([0.0], fiber_radius_um * np.cos(angles_rad)))
y_offsets_um = np.concatenate(([0.0], fiber_radius_um * np.sin(angles_rad)))
field_labels = ["center"] + [f"{a} deg" for a in range(0, 360, 45)]

pupil_px, pupil_py = emar_utils.sample_hexapolar_pupil(
    n_rings=5, include_center=True
)
n_pupil = len(pupil_px)
assert n_pupil == 91, "Unexpected pupil sampling size."


def trace_image(ln, center_y_um):
    """Trace all nine field positions and return a hull and diagnostics."""
    hx_values = (x_offsets_um / fiber_radius_um) * hx_max
    hy_values = ((center_y_um + y_offsets_um) / fiber_radius_um) * hy_max

    intersections = []
    n_errors = 0
    n_vignetted = 0
    n_nonfinite = 0

    for field_label, hx, hy in zip(field_labels, hx_values, hy_values):
        for px, py in zip(pupil_px, pupil_py):
            ray = ln.zGetTrace(
                waveNum=wave_num,
                mode=0,
                surf=surf,
                hx=float(hx),
                hy=float(hy),
                px=float(px),
                py=float(py),
            )
            error, vignette = ray[0], ray[1]
            if error != 0:
                n_errors += 1
            if vignette != 0:
                n_vignetted += 1
            xy = np.array([ray[2], ray[3]], dtype=float)
            if not np.isfinite(xy).all():
                n_nonfinite += 1
            if error == 0 and vignette == 0 and np.isfinite(xy).all():
                intersections.append(xy)

    if len(intersections) < 3:
        raise RuntimeError("Insufficient valid rays to construct a convex hull.")

    points = np.asarray(intersections)
    hull = ConvexHull(points)
    polygon = Polygon(points[hull.vertices])
    diagnostics = {
        "Valid": len(points),
        "Errors": n_errors,
        "Vignetted": n_vignetted,
        "Nonfinite": n_nonfinite,
    }
    return polygon, diagnostics


# ============================================================
# Zemax separation scan
# ============================================================

if not zemax_file.is_file():
    raise FileNotFoundError(f"Zemax file not found: {zemax_file}")

ln = pyz.createLink()
if ln is None:
    raise RuntimeError("Could not establish a PyZDDE connection to Zemax.")

rows = []
try:
    ln.zLoadFile(str(zemax_file))
    ln.zSetConfig(base_config)
    ln.zGetUpdate()
    wavelength_um = ln.zGetWave(wave_num).wavelength

    print("\nPolarimetric separation scan")
    print(f"Order / configuration : {order} / {base_config}")
    print(f"Wavelength            : {wavelength_um:.9f} um")
    print(f"Fields / pupil rays   : {len(field_labels)} / {n_pupil}")
    print(f"Separation range      : {d_values_um[0]:.1f}–{d_values_um[-1]:.1f} um")

    for d_um in d_values_um:
        upper, upper_diag = trace_image(ln, +d_um / 2.0)
        lower, lower_diag = trace_image(ln, -d_um / 2.0)

        overlap_mm2 = upper.intersection(lower).area
        gap_um = upper.distance(lower) * 1000.0
        
        area_upper_mm2 = upper.area
        area_lower_mm2 = lower.area
        
        rows.append({
            "d [um]": float(d_um),
            "Upper area [um²]": area_upper_mm2 * 1e6,
            "Lower area [um²]": area_lower_mm2 * 1e6,
            "Overlap area [um²]": overlap_mm2 * 1e6,
            "Upper overlap [%]": 100.0 * overlap_mm2 / area_upper_mm2,
            "Lower overlap [%]": 100.0 * overlap_mm2 / area_lower_mm2,
            "Gap [um]": gap_um,
            "Valid rays": upper_diag["Valid"] + lower_diag["Valid"],
            "Errors": upper_diag["Errors"] + lower_diag["Errors"],
            "Vignetted": upper_diag["Vignetted"] + lower_diag["Vignetted"],
            "Nonfinite": upper_diag["Nonfinite"] + lower_diag["Nonfinite"],
        })
finally:
    pyz.closeLink()

scan_df = pd.DataFrame(rows)

# ============================================================
# Compact diagnostics and results
# ============================================================

print("\nRay-tracing diagnostics")
print(f"Traced     : {len(d_values_um) * 2 * len(field_labels) * n_pupil}")
print(f"Valid      : {scan_df['Valid rays'].sum()}")
print(f"Errors     : {scan_df['Errors'].sum()}")
print(f"Vignetted  : {scan_df['Vignetted'].sum()}")
print(f"Nonfinite  : {scan_df['Nonfinite'].sum()}")

print("\nGeometric overlap by separation")

columns = [
    "d [um]",
    "Overlap area [um²]",
    "Upper overlap [%]",
    "Lower overlap [%]",
    "Gap [um]",
]

print(
    scan_df[columns].to_string(
        index=False,
        float_format=lambda v: f"{v:.4f}"
    )
)

# Numerical tolerance is only for floating-point roundoff, not a physical
# allowance. A zero overlap area also includes point/edge contact.
area_tolerance_um2 = 1e-8
no_overlap = scan_df[scan_df["Overlap area [um²]"] <= area_tolerance_um2]
if no_overlap.empty:
    print("\nNo zero-area overlap found in the scanned range.")
else:
    first_zero_d = no_overlap.iloc[0]["d [um]"]
    print(f"\nFirst sampled separation with zero overlap area: {first_zero_d:.1f} um")
    print("This is a scan-grid result, not a refined critical separation.")

# ============================================================
# Overlap versus entrance separation
# ============================================================

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(scan_df["d [um]"], scan_df["Upper overlap [%]"], "o-", label="Upper")
ax.plot(scan_df["d [um]"], scan_df["Lower overlap [%]"], "s--", label="Lower")
ax.set_ylim(bottom=0)
apply_scientific_style(
    ax,
    xlabel=r"Center-to-center separation $d$ [$\mu$m]",
    ylabel="Convex-hull overlap [%]",
)
ax.legend(frameon=False)
plt.tight_layout()
plt.show()
