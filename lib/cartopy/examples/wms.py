__tags__ = ['Web services']
"""
Interactive WMS (Web Map Service)
---------------------------------

This example demonstrates the interactive pan and zoom capability
supported by an OGC web services Web Map Service (WMS) aware axes.

"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


def main():
    ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
    ax.coastlines()

    ax.add_wms(wms='http://vmap0.tiles.osgeo.org/wms/vmap0',
               layers=['basic'])

    plt.show()


if __name__ == '__main__':
    main()
