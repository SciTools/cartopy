What's new in Cartopy 0.7
*************************

:Release: 0.7.0
:Date: 21 Mar 2013

* Carwyn Pelley added support for 2D arrays of points to :meth:`cartopy.crs.CRS.transform_points`. (:pull:`192`)
* Phil Elson added control for the gridlines and tick labels drawn with
  :meth:`cartopy.mpl.geoaxes.GeoAxes.gridlines`. (:pull:`238`)
* Various documentation enhancements have been added. (:pull:`247`, :pull:`244` :pull:`240` and :pull:`242`)

This is a quick release which targets two very specific requirements. The goals outlined in the development plan at
``v0.6`` still remain the primary target for ``v0.8`` and beyond.


What's new in Cartopy 0.6
*************************

:Release: 0.6.0
:Date: 19 Feb 2013

* Patrick Peglar added the ability to draw ticks for some limited projections
  when using the :py:func:`~cartopy.mpl.geoaxes.GeoAxes.gridlines` method on an Axes.

* Phil Elson and Carwyn Pelley extended the cartopy documentation to include
  new tutorials such as :ref:`using_the_shapereader`.

* Ian Edwards :doc:`added a new example <examples/favicon>` to create a favicon for cartopy.

* Phil Elson :doc:`added a new example <examples/hurricane_katrina>` to show polygon analysis
  and visualisation with Shapely and Cartopy.

* Edward Campbell added a new :py:class:`cartopy.crs.EuroPP` projection for UTM zone 32.

* Andrew Dawson added a ``central_longitude`` keyword for the Stereographic family of projections.

* Phil Elson added a :py:class:`~cartopy.io.Downloader` class which allows
  automatic downloading of shapefiles (currently from Natural Earth and GSHHS).
  The extension requires no user action and can be configured via the :py:data:`cartopy.config` dictionary.


Development plans for Cartopy 0.7 and beyond
============================================

* Improve the projection definitions to support better control over datum definitions
  and consider adding WKT support (:issue:`ticket <153>`).

* Begin work on vector field support (barbs, quiver, streamlines etc.).

* Continue identifying and implementing performance enhancements (particularly in contour drawing).

* Extend the number of projections for which it is possible to draw tick marks.


-----------


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

.. literalinclude:: /examples/features.py

.. plot:: examples/features.py
