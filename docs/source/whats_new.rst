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

A new features api is now available, see :doc:`tutorials/using_the_shapereader`.

.. literalinclude:: features.py

.. plot:: features.py
