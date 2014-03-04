"""
Rotated pole boxes
------------------

This example demonstrates the way a box is warped when it is defined
in a rotated pole coordinate system.

Try changing the ``box_top`` to ``44``, ``46`` and ``75`` to see the effect
that including the pole in the polygon has.

"""

__tags__ = ['Lines and polygons']
import matplotlib.pyplot as plt

import cartopy.crs as ccrs


def main():
    rotated_pole = ccrs.RotatedPole(pole_latitude=45, pole_longitude=180)

    box_top = 45
    x, y = [-44, -44, 45, 45, -44], [-45, box_top, box_top, -45, -45]

    ax = plt.subplot(211, projection=rotated_pole)
    ax.stock_img()
    ax.coastlines()
    ax.plot(x, y, marker='o', transform=rotated_pole)
    ax.fill(x, y, color='coral', transform=rotated_pole, alpha=0.4)
    ax.gridlines()

    ax = plt.subplot(212, projection=ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines()
    ax.plot(x, y, marker='o', transform=rotated_pole)
    ax.fill(x, y, transform=rotated_pole, color='coral', alpha=0.4)
    ax.gridlines()

    plt.show()


if __name__ == '__main__':
    main()
