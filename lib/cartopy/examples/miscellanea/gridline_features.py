"""
Gridline Features
-----------------

Show how to produce customized gridlines and tick labels with the
:class:`~cartopy.mpl.gridliner.Gridliner` class.

.. _examples-gridliner:
"""
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def main():
    ax = plt.axes(projection=ccrs.Mercator())
    ax.coastlines()

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')

    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xlines = False
    gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.xlabel_style = {'color': 'red', 'weight': 'bold'}

    plt.show()


if __name__ == '__main__':
    main()
