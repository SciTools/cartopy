"""
Adaptive resampling in Cartopy
==============================

When a line that is straight in geographic coordinates is projected onto a map,
it can become highly curved.  Cartopy's resampling algorithm recursively
subdivides each segment, projecting the source-space midpoint and checking
whether it lands close enough to the straight chord between the projected
endpoints. If not, it bisects and recurses adding enough points until
the curve is approximated to within the specified threshold.

The acceptance threshold is ``dest_projection.threshold`` (in target projection
units, usually metres).  A smaller threshold means more subdivisions and a more
accurate curve; a larger threshold means fewer points and a coarser result.

This script demonstrates the effect by projecting the same 60°N latitude line
at three different thresholds.
"""

import matplotlib.pyplot as plt
import numpy as np
import shapely

import cartopy.crs as ccrs


src_crs = ccrs.PlateCarree()
source_line = shapely.LineString([(-100, 60), (100, 60)])


def make_crs(threshold_km):
    """Return a fresh LambertConformal CRS with the given threshold (km)."""
    crs = ccrs.LambertConformal(central_longitude=0, standard_parallels=(20, 60))
    crs.threshold = threshold_km * 1e3
    return crs


dest_crs_default = make_crs(100)

# Three thresholds to compare (coarse, default, fine)
thresholds_km = [1, 100, 1000]
labels = ["1 km  (fine)", "100 km  (default)", "1,000 km  (coarse)"]
colors = ["red", "green", "blue"]

# Project each threshold variant and collect (x, y) arrays of accepted pts
results = []
for thr_km in thresholds_km:
    dest = make_crs(thr_km)
    projected = dest.project_geometry(source_line, src_crs)
    results.append(np.array(projected.geoms[0].coords))

ax = plt.figure(figsize=(12, 7), layout="constrained").add_subplot(
    1, 1, 1, projection=dest_crs_default
)

ax.set_extent([-120, 120, 28, 82], crs=src_crs)
ax.coastlines(linewidth=0.5, color="0.65")
ax.gridlines(linewidth=0.3, color="0.82", linestyle="--")

for coords, thr_km, label, color in zip(results, thresholds_km, labels, colors):
    n = len(coords)
    ax.plot(
        coords[:, 0],
        coords[:, 1],
        color=color,
        lw=2,
        transform=dest_crs_default,
        label=f"threshold = {label}  ({n} pts)",
    )
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        color=color,
        s=40,
        transform=dest_crs_default,
    )

ax.set_title(
    "60°N latitude line: PlateCarree → LambertConformal\n"
    "Dots show accepted sample points at each threshold."
)
ax.legend(loc="lower center")

plt.show()
