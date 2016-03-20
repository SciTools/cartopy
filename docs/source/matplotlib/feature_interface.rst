.. _cartopy_feature_interface:

The cartopy Feature interface
=============================

The :ref:`data copyright, license and attribution  <referencing_copyright>` can be blended on the map using `text annotations (mpl docs) <http://matplotlib.org/users/annotations_intro.html>`_ as shown in `feature_creation <../examples/feature_creation.html>`_. 

.. currentmodule:: cartopy.feature

.. autoclass:: Feature
    :members:
    :undoc-members:


---------

Specific Feature subclasses have been defined for common functionality, such as accessing
Natural Earth or GSHHS shapefiles.


.. autoclass:: ShapelyFeature
    
.. autoclass:: NaturalEarthFeature
    
.. autoclass:: GSHHSFeature

----------

To simplify some very common cases, some pre-defined Features exist as :mod:`cartopy.feature`
constants. The pre-defined Features are all small-scale (1:110m) 
`Natural Earth <http://www.naturalearthdata.com>`_ datasets, and can be added with methods
such as :func:`GeoAxes.add_feature <cartopy.mpl.geoaxes.GeoAxes.add_feature>`:

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

.. note::

    Any Natural Earth dataset can easily be used by creating an
    instance of :class:`cartopy.feature.NaturalEarthFeature`. For
    example::
    
        import cartopy.feature as cfeature
        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                                edgecolor='face',
                                                facecolor=cfeature.COLORS['land'])


A dictionary of some useful colors for drawing features also exists:

.. autodata:: COLORS

For a full list of names in this dictionary:

    >>> import cartopy.feature
    >>> sorted(cartopy.feature.COLORS.keys())
    ['land', 'land_alt1', 'water']
    

------------


Example of using the Feature class with the matplotlib interface
----------------------------------------------------------------

.. literalinclude:: /examples/feature_creation.py

.. plot:: examples/feature_creation.py

