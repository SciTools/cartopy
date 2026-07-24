"""
Gridlines and tick labels
-------------------------

These examples demonstrate how to quickly add longitude
and latitude gridlines and tick labels on a non-rectangular projection.

As you can see on the first example,
longitude labels may be drawn on left and right sides,
and latitude labels may be drawn on bottom and top sides.
Thanks to the ``dms`` keyword, minutes are used when appropriate
to display fractions of degree.

In the second example, labels are still drawn at the map edges
despite its complexity, and some others are also drawn within the map
boundary.

In the third example, labels are drawn only on the left and bottom sides.
"""
import matplotlib.pyplot as plt

import cartopy.crs as ccrs


def main():



    plt.figure(figsize=(6.9228, 3))
    ax1 = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
    ax1.coastlines(resolution='110m')
    ax1.gridlines(draw_labels=True)
    ax1.zebra_frame(use_extent=True)

    plt.figure(figsize=(7, 3))
    ax2 = plt.axes(projection=ccrs.PlateCarree())
    ax2.coastlines(resolution='110m')
    gl = ax2.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    ax2.zebra_frame(colors=['blue', 'green'])
    plt.show()
    print('Done')


if __name__ == '__main__':
    main()
