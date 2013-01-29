.. _cartopy_developer_interfaces:

Plotting, using the cartopy shapereader
=======================================

Plotting shapefiles becomes easy with use of the cartopy shapereader.

Shapefiles can be filtered if neccessary and various data is available, from
cultural (roads, cities, ...) to physical (rivers, coastlines, ...).

Loading shapefiles from various sources can be done in a few simple steps,
see :doc:`using_the_shapereader`.

In addition to the geometries (concept from
`shapely <http://pypi.python.org/pypi/Shapely>`_-the package underneath, used
for handling geometries) which is passed to the
:func:`~cartopy.mpl.geoaxes.GeoAxes.add_geometries` method to
add geometries to the current axes, there is the Feature class, passed to
:func:`~cartopy.mpl.geoaxes.GeoAxes.add_feature`, to add features to
the current axes.  A feature consists of geometrie(s) and any additional
information to plot these geometries (facecolor and crs for example).

Feature objects (collections of lines and polygons) act as geometry wrappers,
allowing ease of use.  Though these are not neccessary to make functionality
available to the user, they allow a simpler interface to be used and written.

Pre-defined features exist for the small-scale (1:110m)
`Natural Earth <http://www.naturalearthdata.com>`_ datasets detailed
below:

=======================================  ==================================================
Name                                     Description
=======================================  ==================================================
.. py:data:: cartopy.feature.BORDERS     Country boundaries.
.. py:data:: cartopy.feature.COASTLINE   Coastline, including major islands.
.. py:data:: cartopy.feature.LAKES       Natural and artificial lakes.
.. py:data:: cartopy.feature.LAND        Land polygons, including major islands.
.. py:data:: cartopy.feature.OCEAN       Ocean polygons.
.. py:data:: cartopy.feature.RIVERS      Single-line drainages, including lake centerlines.
=======================================  ==================================================

But any Natural Earth dataset can easily be used by creating an
instance of :class:`cartopy.feature.NaturalEarthFeature`.


.. currentmodule:: cartopy.feature

.. autoclass:: Feature
    :members:
    :undoc-members:
.. autoclass:: ShapelyFeature
    :members:
    :undoc-members:
.. autoclass:: NaturalEarthFeature
    :members:
    :undoc-members:
.. autoclass:: GSHHSFeature
    :members:
    :undoc-members:

Example1: Plotting countries, rivers and cities using geometries directly
-------------------------------------------------------------------------

.. literalinclude:: plot_shape1.py

.. plot:: tutorials/plot_shape1.py

.. note::
    When adding a geometry with
    :func:`~cartopy.mpl.geoaxes.GeoAxes.add_geometries` the coordinate
    reference system of the supplied geometry must be supplied.

Example2: Plotting countries, rivers and lakes using features
-------------------------------------------------------------

.. literalinclude:: plot_shape2.py

.. plot:: tutorials/plot_shape2.py

.. note::
    Features include all the necessary information for plotting when adding
    them with the :func:`~cartopy.mpl.geoaxes.GeoAxes.add_feature` method.

.. seealso::
   :doc:`../matplotlib/geoaxes`
