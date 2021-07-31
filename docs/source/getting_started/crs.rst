
Coordinate reference systems
----------------------------

The :class:`cartopy.crs.CRS` class is the very core of cartopy, all coordinate reference systems
in cartopy have :class:`~cartopy.crs.CRS` as a parent class, meaning all projections have
the interface described below.

.. autoclass:: cartopy.crs.CRS
   :members:
   :noindex:

   .. data:: globe

      The :class:`~cartopy.crs.Globe` instance of this CRS.


The :class:`~cartopy.crs.Globe` class is used to encapsulate the underlying sphere or ellipsoid of any cartopy CRS.
All CRSs have an associated :class:`~cartopy.crs.Globe`, though often it is just the default :class:`~cartopy.crs.Globe`
which represents the reference ellipsoid (i.e. "wgs84").

.. autoclass:: cartopy.crs.Globe
   :members:
   :exclude-members: to_proj4_params
   :noindex:


The most common :class:`~cartopy.crs.CRS` subclass is itself another abstract class;
the :class:`cartopy.crs.Projection` class represents a 2 dimensional coordinate system
which could be drawn directly as a map (i.e. on a flat piece of paper). :class:`~cartopy.crs.Projection` is the parent class of
all projections in the :ref:`cartopy_projections`.

.. autoclass:: cartopy.crs.Projection
   :members:
   :noindex:


There are a few non-:class:`~cartopy.crs.Projection` subclasses. These represent
coordinate reference systems which are 3 dimensional and could not be drawn directly on a piece of paper.

.. autoclass:: cartopy.crs.Geodetic
   :noindex:

.. autoclass:: cartopy.crs.Geocentric
   :noindex:

.. autoclass:: cartopy.crs.RotatedGeodetic
   :noindex:

There is also a function for calling epsg.io with a specified code, returning the corresponding cartopy projection, see below.

.. autofunction:: cartopy.crs.epsg
   :noindex:
