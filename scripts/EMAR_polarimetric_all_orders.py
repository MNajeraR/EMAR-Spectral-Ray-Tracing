"""
EMAR: Polarimetric Echelle Spectral Format Analysis.

Computes and visualizes the detector spectral format for echelle orders
60–144, considering two polarimetric images of a single 100-um fiber.

The analysis uses 100 wavelengths per order, nine field points per
polarimetric image, and 91 hexapolar pupil rays per field.

Detector envelopes are cached independently for each order. The script
generates three complementary spectral-format visualizations and evaluates
geometric overlap between adjacent orders.
"""
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import pyzdde.zdde as pyz
from tqdm.auto import tqdm

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))
from utils import emar_utils
from utils.emar_plots import apply_scientific_style

ORDERS = tuple(range(60, 145))
fiber_radius_um = 50.0
fiber_diameter_um = 2.0 * fiber_radius_um
edge_separation_um = 250.0
d_um = fiber_diameter_um + edge_separation_um
hx_max, hy_max = 0.286480, 0.285400
surf, row_order, row_wave1 = 54, 1, 2
n_wavelengths = 100
n_comparison_points = 1000
FORCE_RETRACE = False

zemax_file = emar_utils.ZEMAX_DIR / "WP - Con prismas diseñados - camara.zmx"
orders_file = emar_utils.RESULTS_DIR / "echelle_orders_60_144.csv"
output_dir = emar_utils.RESULTS_DIR / f"polarimetric_all_orders_d{d_um:g}um"
cache_dir = output_dir / "envelopes"
full_csv = output_dir / "polarimetric_envelopes_m60_144.csv"
overlap_csv = output_dir / "adjacent_order_clearance.csv"

angles = np.deg2rad(np.arange(0.0, 360.0, 45.0))
x_offsets_um = np.r_[0.0, fiber_radius_um * np.cos(angles)]
y_offsets_um = np.r_[0.0, fiber_radius_um * np.sin(angles)]
pupil_px, pupil_py = emar_utils.sample_hexapolar_pupil(
    n_rings=5, include_center=True
)
assert len(pupil_px) == 91

# Retained only for the detailed 60/61 view.
COLORS = {
    (60, "Upper"): "tab:blue",
    (60, "Lower"): "tab:cyan",
    (61, "Upper"): "tab:orange",
    (61, "Lower"): "tab:red",
}

def trace_order(ln, order, wavelengths):
    """Trace both polarimetric images and collect detector envelopes."""

    base_config = ((order - 60) // 7) + 1

    data = {
        label: {
            "x": [[] for _ in wavelengths],
            "y": [[] for _ in wavelengths],
        }
        for label in ("Upper", "Lower")
    }

    image_centers = (
        ("Upper", d_um),
        ("Lower", 0.0),
    )

    total_fields = 2 * len(x_offsets_um)

    wavelength_min = min(wavelengths)
    wavelength_max = max(wavelengths)

    progress = tqdm(
        total=total_fields,
        desc=f"m={order}",
        unit="field",
        leave=False,
        dynamic_ncols=False,
        ncols=105,
        mininterval=1.0,
        bar_format="{desc} |{bar:25}| {n_fmt}/{total_fmt} ({percentage:.0f}%)",
    )

    try:
        for label, center_y_um in image_centers:

            hx_values = x_offsets_um / fiber_radius_um * hx_max
            hy_values = (
                (center_y_um + y_offsets_um)
                / fiber_radius_um * hy_max
            )

            for field_index, (hx, hy) in enumerate(
                zip(hx_values, hy_values), start=1
            ):

                progress.set_description(
                    f"m={order} | {label:5s} | "
                    f"Field {field_index}/9 | "
                    f"lambda={wavelength_min:.4f}-{wavelength_max:.4f} um",
                    refresh=False,
                )

                # Trace all wavelengths and pupil rays for this field.
                with redirect_stdout(io.StringIO()):
                    result = emar_utils.trace_all_order(
                        ln=ln,
                        base_config=base_config,
                        target_order=order,
                        wavelengths=wavelengths,
                        surf=surf,
                        px=pupil_px,
                        py=pupil_py,
                        row_order=row_order,
                        row_wave1=row_wave1,
                        hx=float(hx),
                        hy=float(hy),
                    )

                spots = result["spots"]

                if len(spots) != len(wavelengths):
                    raise RuntimeError(
                        f"Unexpected wavelength count for m={order}"
                    )

                for i, wavelength in enumerate(wavelengths):

                    spot = spots[float(wavelength)]

                    x = np.asarray(spot["x"], dtype=float)
                    y = np.asarray(spot["y"], dtype=float)

                    if x.shape != y.shape:
                        raise RuntimeError(
                            "Mismatched detector coordinate shapes"
                        )

                    good = np.isfinite(x) & np.isfinite(y)

                    if good.any():
                        data[label]["x"][i].append(x[good])
                        data[label]["y"][i].append(y[good])

                # Update only after completing the entire field.
                progress.update(1)

    finally:
        progress.close()

    rows = []

    for label, samples in data.items():

        for i, wavelength in enumerate(wavelengths):

            if not samples["x"][i]:
                raise RuntimeError(
                    f"No valid rays: m={order}, "
                    f"{label}, lambda={wavelength}"
                )

            x = np.concatenate(samples["x"][i])
            y = np.concatenate(samples["y"][i])

            rows.append({
                "order": order,
                "image": label,
                "wavelength_um": wavelength,
                "x_centroid_mm": np.mean(x),
                "x_min_mm": np.min(x),
                "x_max_mm": np.max(x),
                "y_centroid_mm": np.mean(y),
                "y_min_mm": np.min(y),
                "y_max_mm": np.max(y),
                "n_valid_rays": len(x),
            })

    return pd.DataFrame(rows)


def interpolate_envelope(df, x_common):
    """Interpolate Y extrema using wavelength-wise spectral X centroids."""
    ordered = df.sort_values("x_centroid_mm")
    x = ordered["x_centroid_mm"].to_numpy()
    if np.any(np.diff(x) <= 0):
        raise RuntimeError("Nonmonotonic/duplicate spectral X centroids")
    return {
        key: np.interp(x_common, x, ordered[key].to_numpy())
        for key in ("y_min_mm", "y_max_mm", "y_centroid_mm")
    }


def style_axes(ax, xlabel, ylabel, *, equal=False):
    """Keep axis typography consistent with the project's scientific plots."""
    apply_scientific_style(ax)
    ax.set_title("")
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.tick_params(axis="both", which="major", labelsize=16,
                   direction="in", top=True, right=True,
                   length=6, width=1.3)
    ax.tick_params(axis="both", which="minor", direction="in",
                   top=True, right=True, length=3)
    ax.minorticks_on()
    ax.set_aspect("equal" if equal else "auto", adjustable="box")




def cache_path(order):
    return cache_dir / f"m{order}_d{d_um:g}um_{n_wavelengths}w_9fields.csv"


def load_or_trace_all():
    """Reuse complete per-order CSVs and trace only missing orders."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(orders_file).set_index("order")
    frames = []
    ln = None
    try:
        for order in ORDERS:
            path = cache_path(order)
            if path.is_file() and not FORCE_RETRACE:
                frame = pd.read_csv(path)
                expected = {(order, "Upper"), (order, "Lower")}
                actual = set(zip(frame["order"], frame["image"]))
                if len(frame) != 2 * n_wavelengths or actual != expected:
                    raise ValueError(f"Invalid cached order file: {path}")
                print(f"m={order}: loaded from cache")
            else:
                if ln is None:
                    if not zemax_file.is_file():
                        raise FileNotFoundError(zemax_file)
                    ln = pyz.createLink()
                    if ln is None:
                        raise RuntimeError("Could not connect to Zemax")
                    status = ln.zLoadFile(str(zemax_file))
                    if isinstance(status, (int, float)) and status != 0:
                        raise RuntimeError(f"Zemax load failed: {status}")
                row = table.loc[order]
                reference = np.asarray(
                    [row[f"wave_{i}"] for i in range(1, 12)], dtype=float
                )
                wavelengths = np.linspace(reference[0], reference[-1],
                                          n_wavelengths)
                frame = trace_order(ln, order, wavelengths)
                frame.to_csv(path, index=False)
                print(f"m={order}: saved")
            frames.append(frame)
    finally:
        if ln is not None:
            ln.close()
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(full_csv, index=False)
    return df


def plot_trace(ax, df, order, label, color, *, alpha=0.18,
               linewidth=0.8, linestyle="-"):
    sub = df[(df["order"] == order) & (df["image"] == label)]
    sub = sub.sort_values("x_centroid_mm")
    x = sub["x_centroid_mm"].to_numpy()
    ylo = sub["y_min_mm"].to_numpy()
    yhi = sub["y_max_mm"].to_numpy()
    yc = sub["y_centroid_mm"].to_numpy()
    ax.fill_between(x, ylo, yhi, color=color, alpha=alpha, linewidth=0)
    ax.plot(x, yc, color=color, linewidth=linewidth, linestyle=linestyle)


def order_colormap():
    # Red at m=60; blue at m=144.
    return plt.get_cmap("coolwarm_r"), mcolors.Normalize(60, 144)


def add_order_colorbar(fig, ax, cmap, norm):
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.025)
    cb.set_label("Echelle order", fontsize=16)
    cb.set_ticks([60, 80, 100, 120, 144])
    cb.ax.tick_params(labelsize=13)
    return cb


FIGURE_SIZE = (9, 9)


def draw_full_format(df):
    """Full spectral-format overview; deliberately enlarged vertical scale."""
    fig = plt.figure(figsize=FIGURE_SIZE)
    ax = fig.add_axes([0.14, 0.13, 0.70, 0.76])
    style_axes(ax, "Detector X [mm]", "Detector Y [mm]")
    cmap, norm = order_colormap()
    for order in ORDERS:
        color = cmap(norm(order))
        for label in ("Lower", "Upper"):
            plot_trace(ax, df, order, label, color,
                       alpha=0.15, linewidth=0.75,
                       linestyle="-" if label == "Upper" else "--")
    from matplotlib.cm import ScalarMappable
    colorbar_ax = fig.add_axes([0.87, 0.13, 0.025, 0.76])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=colorbar_ax)
    cb.set_label("Order", fontsize=16)
    cb.set_ticks([60, 80, 100, 120, 144])
    cb.ax.tick_params(labelsize=13)
    fig.savefig(output_dir / "01_full_spectral_format.png", dpi=300)
    return fig


def draw_adjacent_zoom(df):
    """Same two clean inset panels as experiment 19 for m60 and m61."""
    pair = df[df["order"].isin((60, 61))]
    lo = max(pair.loc[pair.order == m, "x_centroid_mm"].min()
             for m in (60, 61))
    hi = min(pair.loc[pair.order == m, "x_centroid_mm"].max()
             for m in (60, 61))
    if hi <= lo:
        raise RuntimeError("Orders 60/61 have no common X interval")
    x_mm = np.linspace(lo, hi, n_comparison_points)
    interpolated = {}
    for order in (60, 61):
        for label in ("Upper", "Lower"):
            sub = pair[(pair.order == order) & (pair.image == label)]
            interpolated[(order, label)] = interpolate_envelope(sub, x_mm)

    fig = plt.figure(figsize=FIGURE_SIZE)
    ax = fig.add_axes([0.14, 0.13, 0.82, 0.76])
    style_axes(ax, "Detector X [mm]", "Detector Y [mm]")
    for key, env in interpolated.items():
        color = COLORS[key]
        ax.fill_between(x_mm, env["y_min_mm"], env["y_max_mm"],
                        color=color, alpha=0.23, linewidth=0)
        ax.plot(x_mm, env["y_centroid_mm"], color=color, linewidth=1.4,
                label=f"m={key[0]} {key[1]}")
    ax.legend(frameon=False, fontsize=11, loc="upper center", ncol=4,
              bbox_to_anchor=(0.5, 1.01))

    for i, (f0, f1) in enumerate(((0.00, 0.20), (0.80, 1.00))):
        left = x_mm[0] + f0 * (x_mm[-1] - x_mm[0])
        right = x_mm[0] + f1 * (x_mm[-1] - x_mm[0])
        mask = (x_mm >= left) & (x_mm <= right)
        ylo = min(np.min(env["y_min_mm"][mask])
                  for env in interpolated.values())
        yhi = max(np.max(env["y_max_mm"][mask])
                  for env in interpolated.values())
        padding = max(0.01, 0.12 * (yhi - ylo))
        inset = inset_axes(
            ax, width="35%", height="35%",
            loc="lower left" if i == 0 else "lower right",
            borderpad=1.0
        )
        for key, env in interpolated.items():
            color = COLORS[key]
            inset.fill_between(x_mm, env["y_min_mm"], env["y_max_mm"],
                               color=color, alpha=0.23, linewidth=0)
            inset.plot(x_mm, env["y_centroid_mm"],
                       color=color, linewidth=1.2)
        inset.set_xlim(left, right)
        inset.set_ylim(ylo - padding, yhi + padding)
        inset.set_aspect("auto")
        inset.set_xticks([])
        inset.set_yticks([])
        inset.set_xlabel("")
        inset.set_ylabel("")
        for spine in inset.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.0)
        ax.add_patch(Rectangle(
            (left, ylo - padding), right - left,
            yhi - ylo + 2 * padding,
            fill=False, edgecolor="0.45",
            linewidth=0.9, linestyle="--", zorder=5
        ))
    fig.savefig(output_dir / "02_orders_60_61_zoom.png", dpi=300)
    return fig


def draw_detector_aspect(df):
    """Monochrome detector view with true physical X/Y scale."""
    fig = plt.figure(figsize=FIGURE_SIZE, facecolor="black")
    ax = fig.add_axes([0.13, 0.13, 0.82, 0.82])
    ax.set_facecolor("black")
    for order in ORDERS:
        for label in ("Lower", "Upper"):
            sub = df[(df["order"] == order) & (df["image"] == label)]
            sub = sub.sort_values("x_centroid_mm")
            ax.plot(sub["x_centroid_mm"], sub["y_centroid_mm"],
                    color="white", linewidth=0.55, alpha=0.9)

    ax.set_xlabel("X [mm]", fontsize=16, color="white")
    ax.set_ylabel("Y [mm]", fontsize=16, color="white")
    ax.set_title("")
    ax.tick_params(axis="both", which="both", colors="white",
                   direction="in", top=True, right=True)
    for spine in ax.spines.values():
        spine.set_color("white")
    # Equal X/Y units and a square plotting area: no rectangular detector.
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    span = max(x1 - x0, y1 - y0)
    xc, yc = (x0 + x1) / 2, (y0 + y1) / 2
    ax.set_xlim(xc - span / 2, xc + span / 2)
    ax.set_ylim(yc - span / 2, yc + span / 2)
    ax.set_aspect("equal", adjustable="box")
    fig.savefig(output_dir / "03_detector_true_aspect.png",
                dpi=300, facecolor=fig.get_facecolor())
    return fig


def evaluate_adjacent_clearance(df):
    """Check sampled Y-envelope clearance for every adjacent-order image pair.

    Negative clearance indicates overlap of the sampled Y envelopes at
    a shared detector X. It is not a full two-dimensional ray intersection.
    """
    rows = []
    labels = ("Upper", "Lower")
    for order_a, order_b in zip(ORDERS[:-1], ORDERS[1:]):
        for label_a in labels:
            a = df[(df["order"] == order_a) & (df["image"] == label_a)]
            a = a.sort_values("x_centroid_mm")
            xa = a["x_centroid_mm"].to_numpy()
            if np.any(np.diff(xa) <= 0):
                raise ValueError(f"Nonmonotonic X centroids: m={order_a}, {label_a}")
            for label_b in labels:
                b = df[(df["order"] == order_b) & (df["image"] == label_b)]
                b = b.sort_values("x_centroid_mm")
                xb = b["x_centroid_mm"].to_numpy()
                if np.any(np.diff(xb) <= 0):
                    raise ValueError(f"Nonmonotonic X centroids: m={order_b}, {label_b}")
                lo = max(xa[0], xb[0])
                hi = min(xa[-1], xb[-1])
                if hi <= lo:
                    rows.append({
                        "order_a": order_a, "image_a": label_a,
                        "order_b": order_b, "image_b": label_b,
                        "min_clearance_um": np.nan,
                        "x_at_min_mm": np.nan, "overlap": False,
                        "common_x": False,
                    })
                    continue
                x = np.linspace(lo, hi, n_comparison_points)
                ea = interpolate_envelope(a, x)
                eb = interpolate_envelope(b, x)
                # Interval gap: positive for either vertical ordering,
                # negative if the two Y intervals intersect.
                gap_a_above_b = ea["y_min_mm"] - eb["y_max_mm"]
                gap_b_above_a = eb["y_min_mm"] - ea["y_max_mm"]
                clearance_um = 1000.0 * np.maximum(gap_a_above_b,
                                                     gap_b_above_a)
                k = int(np.argmin(clearance_um))
                rows.append({
                    "order_a": order_a, "image_a": label_a,
                    "order_b": order_b, "image_b": label_b,
                    "min_clearance_um": float(clearance_um[k]),
                    "x_at_min_mm": float(x[k]),
                    "overlap": bool(np.any(clearance_um < 0)),
                    "common_x": True,
                })

    report = pd.DataFrame(rows)
    report.to_csv(overlap_csv, index=False)
    overlaps = report[report["overlap"]]
    print(f"\nAdjacent-order pairs checked: {len(report)}")
    if overlaps.empty:
        print("No sampled Y-envelope overlaps found between adjacent orders.")
    else:
        print(f"\nWARNING: {len(overlaps)} adjacent-order image pairs overlap:")
        for row in overlaps.itertuples(index=False):
            print(f"  m={row.order_a} {row.image_a} / "
                  f"m={row.order_b} {row.image_b}: "
                  f"min clearance={row.min_clearance_um:.3f} um "
                  f"at X={row.x_at_min_mm:.4f} mm")
    print(f"Clearance report: {overlap_csv}")
    return report


def main():
    df = load_or_trace_all()
    draw_full_format(df)
    draw_adjacent_zoom(df)
    draw_detector_aspect(df)
    evaluate_adjacent_clearance(df)
    plt.show()
    print(f"Saved figures and data to: {output_dir}")
    return df


if __name__ == "__main__":
    envelopes_df = main()
