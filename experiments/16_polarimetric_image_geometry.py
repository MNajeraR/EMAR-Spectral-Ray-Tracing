import numpy as np
import matplotlib.pyplot as plt
from utils.emar_plots import apply_scientific_style


# ============================================================
# Physical parameters
# ============================================================

fiber_diameter_um = 100.0
fiber_radius_um = fiber_diameter_um / 2.0

# Center-to-center separation between the two polarimetric images
d_um = 95.0
#d_um = 100.0
#d_um = 120.0

# Image centers along Y
center_plus_y_um = +d_um / 2.0
center_minus_y_um = -d_um / 2.0

# Free edge-to-edge separation between the two images
free_separation_um = d_um - fiber_diameter_um

# ============================================================
# Upper polarimetric image
# ============================================================

image_plus_x_um = np.array([
    0.0,                  # center
    +fiber_radius_um,     # +X
    -fiber_radius_um,     # -X
    0.0,                  # +Y
    0.0,                  # -Y
])

image_plus_y_um = np.array([
    center_plus_y_um,                    # center
    center_plus_y_um,                    # +X
    center_plus_y_um,                    # -X
    center_plus_y_um + fiber_radius_um,  # +Y
    center_plus_y_um - fiber_radius_um,  # -Y
])

# ============================================================
# Lower polarimetric image
# ============================================================

image_minus_x_um = np.array([
    0.0,                  # center
    +fiber_radius_um,     # +X
    -fiber_radius_um,     # -X
    0.0,                  # +Y
    0.0,                  # -Y
])

image_minus_y_um = np.array([
    center_minus_y_um,                    # center
    center_minus_y_um,                    # +X
    center_minus_y_um,                    # -X
    center_minus_y_um + fiber_radius_um,  # +Y
    center_minus_y_um - fiber_radius_um,  # -Y
])

# ============================================================
# Check coordinates
# ============================================================

labels = ["center", "+X", "-X", "+Y", "-Y"]

print("\nUpper polarimetric image")
print("-------------------------")

for label, x, y in zip(labels, image_plus_x_um, image_plus_y_um):
    print(f"{label:>6s}: X = {x:7.2f} um, Y = {y:7.2f} um")
    
print("\nLower polarimetric image")
print("-------------------------")

for label, x, y in zip(labels, image_minus_x_um, image_minus_y_um):
    print(f"{label:>6s}: X = {x:7.2f} um, Y = {y:7.2f} um")    

# ============================================================
# Separation check
# ============================================================

print("\nPolarimetric image separation")
print("-----------------------------")
print(f"Fiber-image diameter : {fiber_diameter_um:.2f} um")
print(f"Center separation    : {d_um:.2f} um")
print(f"Free separation      : {free_separation_um:.2f} um")

if free_separation_um > 0:
    print("Status               : separated")
elif np.isclose(free_separation_um, 0.0):
    print("Status               : touching")
else:
    print("Status               : overlapping")
    
    
# ============================================================
# Visualization
# ============================================================

fig, ax = plt.subplots(
    figsize=(7, 7)
)

# Circular images
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

# Representative field positions
ax.scatter(
    image_plus_x_um,
    image_plus_y_um,
    marker="x",
    linewidths=1.5,
    s=50,
)

ax.scatter(
    image_minus_x_um,
    image_minus_y_um,
    marker="x",
    linewidths=1.5,
    s=50,
)

# Preserve physical geometry
ax.set_aspect(
    "equal",
    adjustable="box",
)

# EMAR scientific style
apply_scientific_style(
    ax,
    xlabel=r"X [$\mu$m]",
    ylabel=r"Y [$\mu$m]",
)

# ax.set_title(
#     "Polarimetric Image Geometry",
#     fontsize=18,
# )

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