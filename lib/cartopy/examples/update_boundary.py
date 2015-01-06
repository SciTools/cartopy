"""
Update Boundary
------------------

This example demonstrates how to modify the boundary of a given axes.
"""

__tags__ = ['Lines and polygons']
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
import matplotlib.path as mpath

def update_boundary(ax, path):
    """Modifies the boundary of a given axes."""
    background = mpatches.PathPatch(path, edgecolor='none', facecolor='white',
                                    zorder=-1, clip_on=False, transform=ax.transData)
    background.orig_path = path
    background.reclip = True
    ax.background_patch.remove()
    ax.background_patch = background
    ax.patch = background

    outline = mpatches.PathPatch(path, facecolor='none', edgecolor='k',
                                 zorder=2.5, clip_on=False,
                                 transform=ax.transData)
    outline.orig_path = path
    outline.reclip = True
    ax.outline_patch.remove()
    ax.outline_patch = outline

    with ax.hold_limits():
        ax.add_patch(outline)
        ax.add_patch(background)

def main():
    fig = plt.figure()
    ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
    ax2 = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())

    ax1.coastlines(); ax2.coastlines()

    star_path = mpath.Path.unit_regular_star(5, 0.5)
    star_path = mpath.Path(star_path.vertices.copy() * 80, star_path.codes.copy())
    update_boundary(ax2, star_path)

    plt.show()

if __name__ == '__main__':
    main()
