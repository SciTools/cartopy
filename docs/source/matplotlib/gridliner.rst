Cartopy map gridlines and tick labels 
=====================================

The :class:`~cartopy.mpl.gridliner.Gridliner` instance, often created by calling the 
:meth:`cartopy.mpl.geoaxes.GeoAxes.gridlines` method on a 
:class:`cartopy.mpl.geoaxes.GeoAxes` instance, has a variety of attributes which can be
used to determine draw time behaviour of the gridlines and labels.

.. important::

    The current :class:`~cartopy.mpl.gridliner.Gridliner` interface is likely to undergo
    a significant change in the versions following v0.6 in order to fix some of the underying
    limitations of the current implementation.
    

.. autoclass:: cartopy.mpl.gridliner.Gridliner
    :members:
    :undoc-members:

    
    
The following contrived example makes use of many of the features of the Gridliner
class to produce customized gridlines and tick labels:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import cartopy.crs as ccrs
    
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    
    
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
