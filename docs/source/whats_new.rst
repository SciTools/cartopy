What's new in Cartopy 0.5
*************************

:Release: 0.5.0
:Date: 7 Dec 2012

This document explains the new/changed features of Cartopy in version 0.5.

Release 0.5 of Cartopy continues the work to expand the feature-set of
Cartopy to encompass common operations, and provide performance
improvements.


Cartopy 0.5 features
====================

A summary of the main features added with version 0.5:

* An improved feature API to support future expansion and
  sophistication, and a wider range of pre-defined Natural Earth
  datasets.


Incompatible changes
--------------------
None

Deprecations
------------
* The method :meth:`Axes.natural_earth_shp()` has been replaced by the
  method :meth:`Axes.add_feature()` and the :mod:`cartopy.feature`
  module.


Feature API
===========

Features (i.e. collections of lines and polygons) are now added to a map
via the :meth:`~cartopy.mpl.geoaxes.GeoAxes.add_feature`
method.

Pre-defined features exist for the small-scale (1:110m)
`Natural Earth <http://www.naturalearthdata.com>`_ datasets detailed
below:

=========== ==================================================
Name        Description
=========== ==================================================
BORDERS     Country boundaries.
COASTLINE   Coastline, including major islands.
LAKES       Natural and artificial lakes.
LAND        Land polygons, including major islands.
OCEAN       Ocean polygons.
RIVERS      Single-line drainages, including lake centerlines.
=========== ==================================================

But any Natural Earth dataset can easily be used by creating an
instance of :class:`cartopy.feature.NaturalEarthFeature`.

For example, one can draw a map of Africa demonstrating all of the
pre-defined features with the following code:

.. literalinclude:: features.py

.. plot:: features.py
