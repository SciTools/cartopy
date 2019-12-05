"""
WMS (Web Map Service): managing the boundary and token to the wms server
------------------------------------------------------------------------

This example shows how to request several layers from an OGC web services
Web Map Service (WMS) with different projection and how to manage the area.

It shows how to add a token to the request for server which requests it like
ecmwf.

It shows how to set a boundary to better display on the same x limits
(geoographically speaking) the same area with two projections

It shows how to set the boundary color or not draw it.

An http request for this kind of data:
https://apps.ecmwf.int/wms/?token=public&request=getmap&layers=t850_public,boundaries,foreground&srs=EPSG:3857&bbox=-10000000,1000000,15026376,20048966
"""
__tags__ = ['Web services']


import matplotlib.pyplot as plt
from matplotlib.path import Path
import cartopy.crs as ccrs


def full_geo_extent(axe, proj):
    """
    because a plate carree plot is wider for the same height
    the plot 2 may cover a wider window area
    to avoid that:
    get the geo extent that corresponds to the full original axes extent (ax1)
    plot 2 is plotted in the full window width with an extended geographic area
    so that when we cut it to the same xs as plot 1
    we get the same geographic extent as plot 1

    return the geographic extent corresponding to the full window extent
    """
    proj_geo = ccrs.Geodetic()

    ratioxao = (axe.get_position(original=True).width /
                axe.get_position(original=False).width)

    x01, x11 = axe.get_xlim()

    xo01 = x01 - (x11 - x01) * (ratioxao - 1) / 2.
    xo11 = x11 + (x11 - x01) * (ratioxao - 1) / 2.

    y01, y11 = axe.get_ylim()

    xog01, yog01 = proj_geo.transform_point(xo01, y01, proj)
    xog11, yog11 = proj_geo.transform_point(xo11, y11, proj)

    return xog01, xog11, yog01, yog11


def set_outline_color(axe, color='black'):
    """
    Set the color of the outline
    Depending on the cartopy.geoaxes version the outline is
    either a spines or a patch
    spines is from the current cartopy source version
    """
    if 'geo' in axe.spines:
        axe.spines['geo'].set_edgecolor(color)
    elif 'outline_patch' in axe.__dict__:
        axe.outline_patch.set_edgecolor(color)


def main():
    """
    Plot the current temperature et 850 hpa from ECMWF
    in Mercator and Plate Carree
    Adjust the axes to have the longitudes at the same x on both plots
    Use the boundary to clip the second plot to fit the same geographic extent
    If we use the same geo_extent then the second plot would not have the
    longitudes at the same x as the first one, due to the autoscale.
    """
    wms = 'https://apps.ecmwf.int/wms/?token=public'
    layers = ['t850_public', 'foreground', 'boundaries', 'grid']
    geo_extent = (-10, 40, 32, 48)  # x0, x1, y0, y1

    fig = plt.figure(figsize=(13, 9), dpi=72)
    plt.subplots_adjust(top=0.925, bottom=0.0, hspace=0.04)
    plt.suptitle(layers[0])
    proj1 = ccrs.epsg(3857)
    proj2 = ccrs.PlateCarree()

    ax1 = fig.add_subplot(2, 1, 1, projection=proj1)
    ax1.set_extent(geo_extent, crs=proj2)
    ax1.add_wms(wms=wms, layers=layers)
    ax1.set_title('Mercator epsg 3857')

    # the outline is accessed differently according to the version
    set_outline_color(ax1, color='red')

    fig.canvas.draw()  # necessary to get the active position bounding box
    geo_extento1 = full_geo_extent(ax1, proj1)

    ax2 = fig.add_subplot(2, 1, 2, projection=proj2)
    ax2.set_extent(geo_extento1, crs=proj2)
    # the ax2 extent is set so as to correspond to ax1

    ax2.set_title('Plate Carree')

    verts = [
        (geo_extent[0], geo_extent[2]),  # left, bottom
        (geo_extent[0], geo_extent[3]),  # left, top
        (geo_extent[1], geo_extent[3]),  # right, top
        (geo_extent[1], geo_extent[2]),  # right, bottom
        (geo_extent[0], geo_extent[2]),  # ignored
    ]
    codes = [Path.MOVETO] + [Path.LINETO] * 3 + [Path.CLOSEPOLY]
    path = Path(verts, codes)

    ax2.set_boundary(path, transform=ccrs.PlateCarree())
    ax2.add_wms(wms=wms, layers=layers)

    set_outline_color(ax2, color='none')

    plt.show()


if __name__ == '__main__':
    main()
