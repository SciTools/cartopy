"""
Modifying the boundary/neatline of a map in cartopy
---------------------------------------------------

This example demonstrates how to modify the boundary/neatline
of an axes. We construct a star with coordinates in a Plate Carree
coordinate system, and use the star as the outline of the map.

Notice how changing the projection of the map represents a *projected*
star shaped boundary.

"""
__tags__ = ['Miscellanea']
import matplotlib.path as mpath
import matplotlib.pyplot as plt

import cartopy.crs as ccrs


def main():
    ax = plt.axes([0, 0, 1, 1], projection=ccrs.PlateCarree())
    ax.coastlines()

    # Construct a star in longitudes and latitudes.
    star_path = mpath.Path.unit_regular_star(5, 0.5)
    star_path = mpath.Path(star_path.vertices.copy() * 80,
                           star_path.codes.copy())

    # Use the star as the boundary.
    ax.set_boundary(star_path, transform=ccrs.PlateCarree())

    plt.show()


if __name__ == '__main__':
    main()
