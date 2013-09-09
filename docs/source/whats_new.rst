What's new in cartopy 0.9
*************************

:Release: 0.9.0
:Date: 12 September 2013

* We are very pleased to announce that Bill Little was added to the cartopy
  core development team. Bill has made some excellent contributions to cartopy,
  and `his presentation at EuroScipy'13 on
  "Iris & Cartopy" <https://www.euroscipy.org/schedule/presentation/35/>`_
  was voted best talk of the conference.
* Other talks and tutorials during this release cycle include Phil Elson's `talk at SciPy'13
  (with video) <http://conference.scipy.org/scipy2013/presentation_detail.php?id=132>`_,
  `Thomas Lecocq's tutorial at EuroSciPy <https://www.euroscipy.org/schedule/presentation/27/>`_
  and a forthcoming `talk at FOSS4G <http://2013.foss4g.org/conf/programme/presentations/29/>`_.
* Christoph Gohlke updated cartopy to support Windows 7.
* The Plate Carree projection was updated to fully handle arbitrary globe definitions.
* Peter Killick updated the Mercator class' default globe to WGS84. His refactor paved the way
  for some follow on work to fully implement the Google Spherical Mercator (EPSG:3857) projection.

    |image_eyja_volcano|_

    .. |image_eyja_volcano| image:: examples/eyja_volcano_01_00.thumb.png

    .. _image_eyja_volcano: examples/eyja_volcano.html

* The TransverseMercator class saw a tidy up to include several common arguments (:issue:`ticket <309>`)
* Bill Little added the Geostationary projection to allow geolocation of satellite imagery.
  
    |image_geostationary|_

    .. |image_geostationary| image:: examples/geostationary_01_00.thumb.png

    .. _image_geostationary: examples/geostationary.html

* Byron Blay added the :class:`Lambert conformal conic projection <cartopy.crs.LambertConformal>`.


What's new in cartopy 0.8
*************************

:Release: 0.8.0
:Date: 3 June 2013

* Bill Little added support for the OSNI projection and enhanced the image nest capability. (:pull:`263`)
* :class:`cartopy.io.img_nest.Img` has been extended to include a
  :func:`cartopy.io.img_nest.Img.from_world_file` static method for
  easier loading of georeferenced images.
* Phil Elson added a major performance improvement when plotting data from PlateCarree onto a
  PlateCarree map. (:pull:`260`)
* Byron Blay and Richard Hattersley added a :class:`cartopy.crs.Globe` class to encapsulate ellipsoid and optionally
  datum information for CRSs. Globe handling in many projections, including Stereographic, has been added.


What's new in cartopy 0.7
*************************

:Release: 0.7.0
:Date: 21 Mar 2013

* Carwyn Pelley added support for 2D arrays of points to :meth:`cartopy.crs.CRS.transform_points`. (:pull:`192`)
* Phil Elson added control for the gridlines and tick labels drawn with
  :meth:`cartopy.mpl.geoaxes.GeoAxes.gridlines`. (:pull:`238`)
* Various documentation enhancements have been added. (:pull:`247`, :pull:`244` :pull:`240` and :pull:`242`)

This is a quick release which targets two very specific requirements. The goals outlined in the development plan at
``v0.6`` still remain the primary target for ``v0.8`` and beyond.


What's new in cartopy 0.6
*************************

:Release: 0.6.0
:Date: 19 Feb 2013

* Patrick Peglar added the ability to draw ticks for some limited projections
  when using the :py:func:`~cartopy.mpl.geoaxes.GeoAxes.gridlines` method on an Axes.

* Phil Elson and Carwyn Pelley extended the cartopy documentation to include
  new tutorials such as :ref:`using_the_shapereader`.

* Ian Edwards :doc:`added a new example <examples/favicon>` to create a favicon for cartopy.

* Phil Elson :doc:`added a new example <examples/hurricane_katrina>` to show polygon analysis
  and visualisation with Shapely and cartopy.

* Edward Campbell added a new :py:class:`cartopy.crs.EuroPP` projection for UTM zone 32.

* Andrew Dawson added a ``central_longitude`` keyword for the Stereographic family of projections.

* Phil Elson added a :py:class:`~cartopy.io.Downloader` class which allows
  automatic downloading of shapefiles (currently from Natural Earth and GSHHS).
  The extension requires no user action and can be configured via the :py:data:`cartopy.config` dictionary.


Development plans for cartopy 0.7 and beyond
============================================

* Improve the projection definitions to support better control over datum definitions
  and consider adding WKT support (:issue:`ticket <153>`).

* Begin work on vector field support (barbs, quiver, streamlines etc.).

* Continue identifying and implementing performance enhancements (particularly in contour drawing).

* Extend the number of projections for which it is possible to draw tick marks.


-----------


What's new in cartopy 0.5
*************************

:Release: 0.5.0
:Date: 7 Dec 2012

This document explains the new/changed features of cartopy in version 0.5.

Release 0.5 of cartopy continues the work to expand the feature-set of
cartopy to encompass common operations, and provide performance
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
