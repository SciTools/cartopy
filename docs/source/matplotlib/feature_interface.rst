.. _cartopy_feature_interface:

The cartopy Feature interface
=============================

The :ref:`data copyright, license and attribution  <referencing_copyright>` can be blended on the map using `text annotations (mpl docs) <https://matplotlib.org/stable/users/explain/text/annotations.html>`_ as shown in `feature_creation <../gallery/lines_and_polygons/feature_creation.html>`_.

Specific Feature subclasses have been defined for common functionality, such as accessing
Natural Earth or GSHHS shapefiles. A list of these can be found in :ref:`the reference documentation <api.feature>`.

To simplify some very common cases, some pre-defined Features exist as :mod:`cartopy.feature`
constants. The pre-defined Features are all small-scale (1:110m)
`Natural Earth <https://www.naturalearthdata.com>`_ datasets, and can be added with methods
such as :func:`GeoAxes.add_feature <cartopy.mpl.geoaxes.GeoAxes.add_feature>`:

=======================  ================================================================
Name                     Description
=======================  ================================================================
.. py:data:: BORDERS     Country boundaries.

.. py:data:: COASTLINE   Coastline, including major islands.

.. py:data:: LAKES       Natural and artificial lakes.

.. py:data:: LAND        Land polygons, including major islands.

.. py:data:: OCEAN       Ocean polygons.

.. py:data:: RIVERS      Single-line drainages, including lake centerlines.

.. py:data:: STATES      Internal, first-order administrative boundaries (limited to the
                         United States at this scale).
                         Natural Earth have first-order admin boundaries for most
                         countries at the 1:10,000,000 scale; these may be
                         accessed with ``cartopy.feature.STATES.with_scale('10m')``
=======================  ================================================================

.. note::

    Any Natural Earth dataset can be used by creating an
    instance of :class:`cartopy.feature.NaturalEarthFeature`. For
    example::

        import cartopy.feature as cfeature
        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                                edgecolor='face',
                                                facecolor=cfeature.COLORS['land'])


A dictionary of some useful colors for drawing features also exists in :attr:`cartopy.feature.COLORS`.

For a full list of names in this dictionary:

    >>> import cartopy.feature
    >>> sorted(cartopy.feature.COLORS.keys())
    ['land', 'land_alt1', 'water']


------------


Example of using the Feature class with the Matplotlib interface
----------------------------------------------------------------

.. figure:: ../gallery/lines_and_polygons/images/sphx_glr_feature_creation_001.png
   :target: ../gallery/lines_and_polygons/feature_creation.html
   :align: center
   :scale: 50
