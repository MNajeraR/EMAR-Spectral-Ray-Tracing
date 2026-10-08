"""Experiment 19: polarimetric clearance between echelle orders 60 and 61.

One physical 100-um fiber, two polarimetric images separated by d_um in the
input focal plane. Nine source fields per image, 91 pupil rays per field,
100 wavelengths per order. Geometric Y envelopes are compared at common X.

After the first Zemax run, detector envelopes are cached as CSV. Set
FORCE_RETRACE=True to recompute them.
"""

import sys
from pathlib import Path
from tqdm.auto import tqdm
from contextlib import redirect_stdout
import io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyzdde.zdde as pyz
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from utils import emar_utils
from utils.emar_plots import apply_scientific_style

# Numerical parameters
ORDERS = (60, 61)

# fiber_radius_um = 50.0
# fiber_diameter_um = 2.0 * fiber_radius_um

# edge_separation_um = 100.0
# d_um = fiber_diameter_um + edge_separation_um

fiber_radius_um = 50.0
fiber_diameter_um = 2.0 * fiber_radius_um

#edge_separation_um = 550.0
edge_separation_um = 250.0
d_um = fiber_diameter_um + edge_separation_um  # 500 um

hx_max = 0.286480
hy_max = 0.285400

surf = 54
row_order = 1
row_wave1 = 2
n_wavelengths = 100
n_comparison_points = 1000

# FORCE_RETRACE = True
FORCE_RETRACE = False

zemax_file = (emar_utils.ZEMAX_DIR /
              "WP - Con prismas diseñados - camara.zmx")
orders_file = emar_utils.RESULTS_DIR / "echelle_orders_60_144.csv"
cache_file = emar_utils.RESULTS_DIR / (
    f"polarimetric_envelopes_m60_61_d{d_um:g}um_"
    f"{n_wavelengths}w_9fields.csv"
)
clearance_file = emar_utils.RESULTS_DIR / (
    f"polarimetric_clearance_m60_61_d{d_um:g}um.csv"
)
format_figure = emar_utils.RESULTS_DIR / "polarimetric_format_m60_61_zoom.png"
clearance_figure = emar_utils.RESULTS_DIR / "polarimetric_clearance_m60_61.png"

angles = np.deg2rad(np.arange(0.0, 360.0, 45.0))
x_offsets_um = np.r_[0.0, fiber_radius_um * np.cos(angles)]
y_offsets_um = np.r_[0.0, fiber_radius_um * np.sin(angles)]
pupil_px, pupil_py = emar_utils.sample_hexapolar_pupil(
    n_rings=5, include_center=True
)
assert len(pupil_px) == 91

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
        leave=True,
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


def load_or_trace():
    """Use cached detector envelopes, unless retracing is requested."""
    if cache_file.is_file() and not FORCE_RETRACE:
        df = pd.read_csv(cache_file)
        expected = {(order, label) for order in ORDERS
                    for label in ("Upper", "Lower")}
        actual = set(zip(df["order"], df["image"]))
        if expected != actual or len(df) != len(expected) * n_wavelengths:
            raise RuntimeError("Cached CSV does not match the current sampling")
        print(f"Loaded cached envelopes: {cache_file}")
        return df

    if not zemax_file.is_file():
        raise FileNotFoundError(f"Zemax file not found: {zemax_file}")
    if not orders_file.is_file():
        raise FileNotFoundError(f"Order table not found: {orders_file}")

    table = pd.read_csv(orders_file).set_index("order")
    wavelengths_by_order = {}
    for order in ORDERS:
        row = table.loc[order]
        reference = np.asarray(
            [row[f"wave_{i}"] for i in range(1, 12)], dtype=float
        )
        wavelengths_by_order[order] = np.linspace(
            reference[0], reference[-1], n_wavelengths
        )

    ln = pyz.createLink()
    if ln is None:
        raise RuntimeError("Could not connect to Zemax through PyZDDE")
    try:
        status = ln.zLoadFile(str(zemax_file))
        if isinstance(status, (int, float)) and status != 0:
            raise RuntimeError(f"Zemax model load failed: {status}")
        frames = [trace_order(ln, order, wavelengths_by_order[order])
                  for order in ORDERS]
    finally:
        ln.close()

    df = pd.concat(frames, ignore_index=True)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache_file, index=False)
    print(f"Saved detector envelopes: {cache_file}")
    return df


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


def draw_spectral_format(interpolated, x_common):
    """Show full spectral format with two detector-X zoom regions."""
    x_mm = x_common
    fig, ax = plt.subplots(figsize=(12, 7))
    style_axes(ax, "Detector X [mm]", "Detector Y [mm]")

    for key, env in interpolated.items():
        color = COLORS[key]
        ax.fill_between(x_mm, env["y_min_mm"],
                        env["y_max_mm"],
                        color=color, alpha=0.23, linewidth=0)
        ax.plot(x_mm, env["y_centroid_mm"],
                color=color, linewidth=1.4,
                label=f"m={key[0]} {key[1]}")

    ax.legend(frameon=False, fontsize=11, loc="upper center", ncol=4,
              bbox_to_anchor=(0.5, 1.01))
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.13, top=0.92)

    # Zoom regions are located near the left and right ends of common X.
    # Their Y bounds are derived from the sampled envelopes, not hardcoded.
    zoom_positions = [(0.00, 0.20), (0.80, 1.00)]
    for i, (f0, f1) in enumerate(zoom_positions):
        left = x_mm[0] + f0 * (x_mm[-1] - x_mm[0])
        right = x_mm[0] + f1 * (x_mm[-1] - x_mm[0])
        mask = (x_mm >= left) & (x_mm <= right)
        y_values = []
        for env in interpolated.values():
            y_values.extend((env["y_min_mm"][mask],
                             env["y_max_mm"][mask]))
        ymin = min(np.min(a) for a in y_values)
        ymax = max(np.max(a) for a in y_values)
        padding = max(0.01, 0.12 * (ymax - ymin))

        inset = inset_axes(ax, width="35%", height="35%",
                           loc="lower left" if i == 0 else "lower right",
                           borderpad=1.0,
                           )
        for key, env in interpolated.items():
            color = COLORS[key]
            inset.fill_between(x_mm, env["y_min_mm"],
                               env["y_max_mm"],
                               color=color, alpha=0.23, linewidth=0)
            inset.plot(x_mm, env["y_centroid_mm"],
                       color=color, linewidth=1.2)
        inset.set_xlim(left, right)
        inset.set_ylim(ymin - padding, ymax + padding)
        inset.set_aspect("auto")
        inset.set_xticks([])
        inset.set_yticks([])
        
        inset.set_xlabel("")
        inset.set_ylabel("")
        
        for spine in inset.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.)
        ax.add_patch(Rectangle(
            (left, ymin - padding), right - left,
            ymax - ymin + 2 * padding,
            fill=False, edgecolor="0.45", linewidth=0.9,
            linestyle="--", zorder=5,
        ))

    fig.savefig(format_figure, dpi=300, bbox_inches="tight")
    return fig


def draw_clearance(x_common, profiles):
    fig, ax = plt.subplots(figsize=(11, 6))
    style_axes(ax, "Detector X [mm]",
               r"Inter-order clearance [$\mu$m]")
    colors = {
        ("Upper", "Upper"): "tab:blue",
        ("Upper", "Lower"): "tab:orange",
        ("Lower", "Upper"): "tab:green",
        ("Lower", "Lower"): "tab:red",
    }
    for labels, clearance in profiles.items():
        ax.plot(x_common, clearance, linewidth=1.6,
                color=colors[labels],
                label=f"m60 {labels[0]} / m61 {labels[1]}")
    ax.axhline(0, color="0.35", linewidth=0.9, linestyle="--")
    # Highlight regions where the clearance is negative.
    for labels, clearance_um in profiles.items():
    
        overlap_mask = clearance_um < 0.0
    
        if np.any(overlap_mask):
            ax.fill_between(
                x_common,
                clearance_um,
                0.0,
                where=overlap_mask,
                interpolate=True,
                color="tab:red",
                alpha=0.30,
            )
    
            min_index = np.argmin(clearance_um)
    
            ax.scatter(
                x_common[min_index],
                clearance_um[min_index],
                color="tab:red",
                s=45,
                zorder=5,
            )
    ax.legend(frameon=False, fontsize=12, loc="best")
    fig.tight_layout()
    fig.savefig(clearance_figure, dpi=300, bbox_inches="tight")
    return fig


def main():
    envelopes_df = load_or_trace()
    x_lo = max(envelopes_df.loc[envelopes_df.order == order,
                                "x_centroid_mm"].min() for order in ORDERS)
    x_hi = min(envelopes_df.loc[envelopes_df.order == order,
                                "x_centroid_mm"].max() for order in ORDERS)
    if x_hi <= x_lo:
        raise RuntimeError("Orders 60/61 have no common spectral X interval")
    x_common = np.linspace(x_lo, x_hi, n_comparison_points)

    interpolated = {}
    for order in ORDERS:
        for label in ("Upper", "Lower"):
            subset = envelopes_df[(envelopes_df.order == order) &
                                  (envelopes_df.image == label)]
            interpolated[(order, label)] = interpolate_envelope(
                subset, x_common
            )

    rows = []
    profiles = {}
    for label_60 in ("Upper", "Lower"):
        for label_61 in ("Upper", "Lower"):
            a = interpolated[(60, label_60)]
            b = interpolated[(61, label_61)]
            if np.median(a["y_centroid_mm"] - b["y_centroid_mm"]) >= 0:
                upper, lower = a, b
            else:
                upper, lower = b, a
            clearance_um = 1000.0 * (
                upper["y_min_mm"] - lower["y_max_mm"]
            )
            profiles[(label_60, label_61)] = clearance_um
            k = int(np.argmin(clearance_um))
            rows.append({
                "Image m60": label_60,
                "Image m61": label_61,
                "Min clearance [um]": clearance_um[k],
                "Median clearance [um]": np.median(clearance_um),
                "X at min [mm]": x_common[k],
                "Overlap": bool(np.any(clearance_um < 0)),
            })

    result_df = pd.DataFrame(rows)
    # Report detector regions with negative clearance.
    for (label_60, label_61), clearance_um in profiles.items():
    
        overlap_mask = clearance_um < 0.0
    
        if not np.any(overlap_mask):
            continue
    
        overlap_indices = np.flatnonzero(overlap_mask)
        segments = np.split(
            overlap_indices,
            np.where(np.diff(overlap_indices) > 1)[0] + 1
        )
    
        print(f"\nOverlap: m60 {label_60} / m61 {label_61}")
    
        for segment in segments:
            x_start = x_common[segment[0]]
            x_end = x_common[segment[-1]]
    
            local_min_index = segment[np.argmin(clearance_um[segment])]
    
            print(
                f"  Detector X interval : "
                f"{x_start:.4f} to {x_end:.4f} mm\n"
                f"  Maximum overlap     : "
                f"{-clearance_um[local_min_index]:.4f} um\n"
                f"  X at maximum overlap: "
                f"{x_common[local_min_index]:.4f} mm"
            )
    result_df.to_csv(clearance_file, index=False)
    print("\nAdjacent-order polarimetric clearance: m=60/61")
    print(f"Input image separation : {d_um:.1f} um")
    print(f"Wavelengths / order    : {n_wavelengths}")
    print(f"Fields / pupil rays    : {len(x_offsets_um)} / {len(pupil_px)}")
    print(f"Common detector X      : {x_lo:.4f} to {x_hi:.4f} mm")
    print(result_df.to_string(index=False,
                              float_format=lambda value: f"{value:.4f}"))
    print("Negative clearance indicates sampled Y-envelope overlap.")

    draw_spectral_format(interpolated, x_common)
    draw_clearance(x_common, profiles)
    plt.show()
    return result_df, envelopes_df


if __name__ == "__main__":
    result_df, envelopes_df = main()
