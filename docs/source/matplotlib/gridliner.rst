Cartopy map gridlines and tick labels
=====================================

The :class:`~cartopy.mpl.gridliner.Gridliner` instance, often created by calling the
:meth:`cartopy.mpl.geoaxes.GeoAxes.gridlines` method on a
:class:`cartopy.mpl.geoaxes.GeoAxes` instance, has a variety of attributes which can be
used to determine draw time behaviour of the gridlines and labels.


In this first example, gridines and tick labels are plotted in a
non-rectangular projection, with most default values and
no tuning of the gridliner attributes:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    rotated_crs = ccrs.RotatedPole(pole_longitude=120.0, pole_latitude=70.0)

    ax = plt.axes(projection=rotated_crs)
    ax.set_extent([-6, 3, 48, 58], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m')
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

    plt.show()


The following contrived example makes use of many of the features of the Gridliner
class to produce customized gridlines and tick labels:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import cartopy.crs as ccrs

    from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                    LatitudeLocator)


    ax = plt.axes(projection=ccrs.Mercator())
    ax.coastlines()

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.left_labels = False
    gl.xlines = False
    gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
    gl.ylocator = LatitudeLocator()
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.ylabel_style = {'size': 15, 'color': 'gray'}
    gl.xlabel_style = {'color': 'red', 'weight': 'bold'}

    plt.show()
