.. _api.crs:

Coordinate reference systems (CRS)
----------------------------------

.. module:: cartopy.crs

The :class:`cartopy.crs.CRS` class is the very core of cartopy, all coordinate reference systems
in cartopy have :class:`~cartopy.crs.CRS` as a parent class.

Base CRS's
~~~~~~~~~~

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_inherited.rst

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

.. module:: cartopy.geodesic

The :mod:`cartopy.geodesic` module defines the :class:`cartopy.geodesic.Geodesic` class which can interface with the Proj
geodesic functions. See the `Proj geodesic page`_ for more background
information.

.. _Proj geodesic page: https://proj.org/geodesic.html

.. autosummary::
   :toctree: generated/

   Geodesic

List of projections
~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   projections
