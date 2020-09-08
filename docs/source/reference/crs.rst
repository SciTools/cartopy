.. _api.crs:

Coordinate reference systems (CRS)
----------------------------------

.. currentmodule:: cartopy.crs

The :class:`cartopy.crs.CRS` class is the very core of cartopy, all coordinate reference systems
in cartopy have :class:`~cartopy.crs.CRS` as a parent class.

Base CRS's
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   CRS
   Globe
   Projection
   Geodetic
   Geocentric
   RotatedGeodetic
   epsg

.. _api.geodesic:

Geodesic calculations
~~~~~~~~~~~~~~~~~~~~~

The :mod:`cartopy.geodesic` module defines the :class:`cartopy.crs.geodesic.Geodesic` class which can interface with the Proj
geodesic functions. See the `Proj geodesic page`_ for more background
information.

.. _Proj geodesic page: https://proj.org/geodesic.html

.. currentmodule:: cartopy

.. autosummary::
   :toctree: generated/

   geodesic.Geodesic


List of projections
~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   projections