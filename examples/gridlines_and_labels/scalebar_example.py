"""
Gridlines scalebar
-------------------------

These examples demonstrate how to quickly add longitude
and latitude gridlines and scalebar on a non-rectangular projection.


In this example, labels are drawn only on the left and bottom sides,\
and a scalebar is created at the bottom-center of the axes.
"""
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
from cartopy.mpl.scalebar import fancy_scalebar


def add_scalebar_to_ax(ax):

    fancy_scalebar(ax,
                   location=(0.5, 0.2),
                   length=5000,
                   unit_name='km',
                   dy=2,
                   max_stripes=3,
                   angle=0)


def main():

    plt.figure(figsize=(7, 3))
    ax2 = plt.axes(projection=ccrs.PlateCarree())
    ax2.coastlines(resolution='110m')
    add_scalebar_to_ax(ax2)
    gl = ax2.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    plt.show()


if __name__ == '__main__':
    main()
