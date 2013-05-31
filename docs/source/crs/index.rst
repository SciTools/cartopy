Coordinate reference systems in Cartopy
---------------------------------------

The :class:`cartopy.crs.CRS` class is the very core of Cartopy, all coordinate reference systems
in Cartopy have :class:`~cartopy.crs.CRS` as a parent class, meaning all projections have
the interface described below.

.. autoclass:: cartopy.crs.CRS
   :members:


The most common :class:`~cartopy.crs.CRS` subclass is itself another abstract class;
the :class:`cartopy.crs.Projection` class represents a 2 dimensional coordinate system
which could be drawn directly as a map (i.e. on a piece of paper). :class:`~cartopy.crs.Projection` is the parent class of
all projections in :ref:`cartopy_projections`.

.. autoclass:: cartopy.crs.Projection
    :members:


There are a few non-:class:`~cartopy.crs.Projection` subclasses. These represent
coordinate reference systems which are 3 dimensional and could not be drawn directly on a piece of paper.

.. autoclass:: cartopy.crs.Geodetic

.. autoclass:: cartopy.crs.Geocentric

.. autoclass:: cartopy.crs.RotatedGeodetic
