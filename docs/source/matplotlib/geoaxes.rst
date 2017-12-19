Cartopy Matplotlib integration reference document
=================================================

The primary class for integrating cartopy into Matplotlib is the GeoAxes, which is a subclass of
a normal Matplotlib Axes. The GeoAxes class adds extra functionality to an axes which is specific
to drawing maps. The majority of the methods which have been specialised from the original Axes
are there to add improved -expected- behaviour, but some are to work around limitations that the
standard Matplotlib axes treats data in a Cartesian way (most of which either have, or will be,
submitted back to the Matplotlib project).


.. autoclass:: cartopy.mpl.geoaxes.GeoAxes
    :members:
    :exclude-members: cla, contour, contourf, draw, format_coord, imshow, pcolor, pcolormesh, scatter
    :undoc-members:


.. class:: cartopy.mpl.geoaxes.GeoAxesSubplot

    See :class:`cartopy.mpl.geoaxes.GeoAxes`.


.. autofunction:: cartopy.mpl.patch.path_to_geos

.. autofunction:: cartopy.mpl.patch.geos_to_path
