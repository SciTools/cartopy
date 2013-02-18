.. _cartopy_developer_interfaces:

The cartopy Feature interface with matplotlib
=============================================
Feature objects are collections of points, lines and polygons.

To add features to the current axes, the reader is refered to the
:func:`~cartopy.mpl.geoaxes.GeoAxes.add_feature` method of the
:class:`~cartopy.mpl.geoaxes.GeoAxes` class.

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

Pre-defined features exist for the small-scale (1:110m)
`Natural Earth <http://www.naturalearthdata.com>`_ datasets:

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

.. note:
    Any Natural Earth dataset can easily be used by creating an
    instance of :class:`cartopy.feature.NaturalEarthFeature`.


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
    Features can include matplotlib artist keywords for plotting, either
    associated with the feature itself, or when they are added with the
    :func:`~cartopy.mpl.geoaxes.GeoAxes.add_feature` method.

.. seealso::
   :doc:`../matplotlib/geoaxes`
