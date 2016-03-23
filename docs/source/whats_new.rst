What's New in cartopy 0.14
==========================

:Release: 0.14.0
:Date: 24th March 2016

Features
--------

* Zachary Tessler and Raj Kesavan added the :class:`cartopy.crs.Sinusoidal` projection,
  allowing MODIS data to be visualised in its native projection. Additionally, a
  prepared :data:`cartopy.crs.Sinusoidal.MODIS` projection has been made available for
  convenience.

* Joseph Hogg and Daniel Atton Beckmann added the :class:`cartopy.geodesic.Geodesic`
  class which wraps the proj.4 geodesic library. This allows users to solve the direct and
  inverse geodesic problems (calculating distances between points etc). It also contains a
  convenience function that returns geodetic circles. This is used by
  :meth:`cartopy.mpl.geoaxes.GeoAxes.tissot` which draws Tissot's indicatrices on the axes.

   |tissot|_

     .. |tissot| image:: examples/tissot_00_00.png

     .. _tissot: examples/tissot.html

* The SRTM3 data source has been changed to the `LP DAAC Data Pool
  <https://lpdaac.usgs.gov/data_access/data_pool>`_. The Data Pool is more
  consistent, fixing several missing tiles, and the data is void-filled.
  Consequently, the :func:`cartopy.srtm.fill_gaps` function has been deprecated
  as it has no purpose within the STRM context. The
  :ref:`SRTM example<examples-srtm_shading>` has also
  been updated to skip the void-filling step. Additionally, this data source
  provides SRTM at a higher resolution of 1 arc-second, which may be accessed
  via :class:`cartopy.io.srtm.SRTM1Source`.

* All downloaders will use secure connections where available. Not
  every service supports this method, and so those will use non-secured HTTP connections
  instead. (See :pull:`736` for full details.)

* Cartopy now supports, and is tested against, matplotlib 1.3 and 1.5 as well as
  numpy 1.7, 1.8 and 1.10.

* Daniel Eriksson added a new example to the gallery:

  |image_aurora|_

    .. |image_aurora| image:: examples/aurora_forecast_00_00.png

    .. _image_aurora: examples/aurora_forecast.html

Incompatible changes
--------------------
* :meth:`cartopy.crs.CRS.transform_point` now issues NaNs when invalid transforms are identified.


Deprecations
------------
* :data:`cartopy.crs.GOOGLE_MERCATOR` has been moved to :data:`cartopy.crs.Mercator.GOOGLE`.


What's new in cartopy 0.13
==========================

:Release: 0.13.0
:Date: 30th June 2015

Features
--------

* Andrea Smith fixed the cartopy CRS class such that 3d transforms such as :class:`cartopy.crs.Geocentric`
  now correctly apply deg2rad and rad2deg. (:pull:`625`)

* Peter Killick fixed the cartopy.crs.Mercator projection for non-zero central longitudes. (:pull:`633`)

* Conversion between matplotlib :class:`matplotlib.path.Path` and
  :class:`shapely.geometry.Geometry <Shapely geometry>` using
  :func:`cartopy.mpl.patch.path_to_geos` and :func:`cartopy.mpl.patch.geos_to_path` now
  handles degenerate point paths.

* Update of tools/feature_download.py to allow mass download of feature data rather than
  on-demand downloading.

* A new example was added to the gallery:

  |image_eccentric_ellipse|_

    .. |image_eccentric_ellipse| image:: examples/eccentric_ellipse_00_00.png

    .. _image_eccentric_ellipse: examples/eccentric_ellipse.html



What's new in cartopy 0.12
==========================

:Release: 0.12.0
:Date: 14th April 2015

Features
--------

* We are very pleased to announce that Elliott Sales de Andrade was added to the cartopy
  core development team. Elliott has added several new projections in this release, as well
  as setting up cartopy's Python 3 testing on TravisCI and generally improving the cartopy
  codebase.

* Installing cartopy became much easier for conda users. A ``scitools`` channel has been
  added which makes getting cartopy and all of its dependencies on Linux, OSX and
  Windows possible with::

     conda install -c scitools cartopy

* Support for Python 3, specifically 3.3 and 3.4, has been added. Some features that depend
  on OWSLib will not be available as it does not support Python 3.

* Two new projections, :class:`~cartopy.crs.AzimuthalEquidistant` and
  :class:`~cartopy.crs.AlbersEqualArea` have been added. See the :ref:`cartopy_projections`
  for the full list of projections now available in cartopy.

* The Web Map Service (WMS) interface has been extended to support on-the-fly reprojection
  of imagery if the service does not support the projection of the map being drawn.
  The following example demonstrates the process by adding WMS imagery to an Interrupted
  Goode Homolosine map - unsurprisingly this WMS service does not provide IGH imagery, so
  cartopy has had to reproject them from a projection the WMS does support:

    .. plot:: examples/wms.py
       :width: 200pt

* Peter Killick added an interface for accessing MapBox tiles using the MapBox
  Developer API. A MapBox client can be created with,
  :class:`~cartopy.io.img_tiles.MapboxTiles` and as with the other imagery from a simple URL
  based imagery service, it can be added to a :class:`~cartopy.mpl.geoaxes.GeoAxes` with the
  :meth:`~cartopy.mpl.geoaxes.GeoAxes.add_image` method. The following example demonstrates the
  interface for another source of imagery:

    .. plot:: examples/image_tiles.py
       :width: 200pt

* Some improvements were made to the geometry transformation algorithm to improve
  the stability of geometry winding. Several cases of geometries being incorrectly
  inverted when transformed have now been resolved. (:pull:`545`)

* Mark Hedley added the ``central_rotated_longitude`` keyword to
  :class:`cartopy.crs.RotatedPole`, which is particularly useful for limited area
  rotated pole models in areas such as New Zealand:

    .. plot::
       :width: 200pt

        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs

        rpole = ccrs.RotatedPole(pole_longitude=171.77,
                                 pole_latitude=49.55,
                                 central_rotated_longitude=180)
        ax = plt.axes(projection=rpole)
        ax.set_global()
        ax.gridlines()
        ax.stock_img()
        ax.coastlines()
        plt.show()

* A new method has been added to the :class:`~cartopy.mpl.geoaxes.GeoAxes` to
  allow control of the neatline of a map drawn with the matplotlib interface.
  The method, :meth:`~cartopy.mpl.geoaxes.GeoAxes.set_boundary`, takes a
  :class:`matplotlib Path<matplotlib.path.Path>` object, which means that
  arbitrary shaped edges can be achieved:

    .. plot:: examples/star_shaped_boundary.py
       :width: 200pt

* A new SRTM3 RasterSource has been implemented allowing interactive pan/zoom
  of 3 arc-second elevation data from the Shuttle Radar Topography Mission.
  The :ref:`SRTM example<examples-srtm_shading>` has also been updated to use the
  new interface.

* New additions to the gallery:

  * |image_un_flag|_

    .. |image_un_flag| image:: examples/un_flag_00_00.thumb.png

    .. _image_un_flag: examples/un_flag.html

  * |image_always_circular_stereo|_

    .. |image_always_circular_stereo| image:: examples/always_circular_stereo_00_00.thumb.png

    .. _image_always_circular_stereo: examples/always_circular_stereo.html

  * |image_tube_stations|_

    .. |image_tube_stations| image:: examples/tube_stations_00_00.thumb.png

    .. _image_tube_stations: examples/tube_stations.html

  * |image_wms|_

    .. |image_wms| image:: examples/wms_00_00.thumb.png

    .. _image_wms: examples/wms.html

  * |image_image_tiles|_

    .. |image_image_tiles| image:: examples/image_tiles_00_00.thumb.png

    .. _image_image_tiles: examples/image_tiles.html



Deprecations
------------
* The SRTM module has been re-factored for simplicity and to take advantage
  of the new :ref:`raster source interface <raster-source-interface>`. Some
  methods have therefore been deprecated and will be removed in future
  releases. The function :func:`cartopy.io.srtm.srtm` has been replaced with
  the :meth:`cartopy.io.srtm.SRTM3Source.single_tile` method. Similarly,
  :func:`cartopy.io.srtm.srtm_composite` and
  :func:`cartopy.io.srtm.SRTM3_retrieve` have been replaced with the
  :meth:`cartopy.io.srtm.SRTM3Source.combined` and
  :meth:`cartopy.io.srtm.SRTM3Source.srtm_fname` methods respectively.

* The :class:`cartopy.io.RasterSource.fetch_raster` interface has been
  changed such that a sequence of :class:`cartopy.io.LocatedImage` must be
  returned, rather than a single image and its associated extent.

* The ``secant_latitudes`` keyword in :class:`cartopy.crs.LambertConformal` has
  been deprecated in favour of ``standard_parallels``.


-----------



What's new in cartopy 0.11
==========================

:Release: 0.11.0
:Date: 19 June 2014


* Richard Hattersley added :func:`~cartopy.crs.epsg` support for generating
  a Cartopy projection at run-time based on the EPSG code of a projected
  coordinate system. This mechanism utilises https://epsg.io/ as a coordinate
  system resource and employs EPSG request caching using
  `pyepsg <https://github.com/rhattersley/pyepsg>`_

* Phil Elson added :class:`~cartopy.io.ogc_clients.WMSRasterSource` which
  provides interactive pan and zoom OGC web services support for a Web Map
  Service (WMS) aware axes. This capability may be added to an axes via the
  :meth:`~cartopy.mpl.geoaxes.GeoAxes.add_wms` method. Generic interactive
  slippy map panning and zooming capability is managed through the new
  :class:`~cartopy.mpl.slippy_image_artist.SlippyImageArtist` and use of the
  :meth:`~cartopy.mpl.geoaxes.GeoAxes.add_raster` method.

* :class:`~cartopy.io.ogc_clients.WMTSRasterSource` was added by Richard
  Hattersley to provide interactive pan and zoom OGC web services support for
  a Web Map Tile Service (WMTS) aware axes, which is available through the
  :meth:`~cartopy.mpl.geoaxes.GeoAxes.add_wmts` method. This includes support
  for the Google Mercator projection and efficient WTMS tile caching. This new
  capability determines how to match up the available tiles projections
  with the target projection and chooses the zoom level to best match the pixel
  density in the rendered image.

.. plot:: examples/wmts.py
    :width: 300pt

* Thomas Lecocq added functionality to :mod:`cartopy.io.srtm` allowing
  intelligent filling of missing elevation data, as well as a function to
  compute elevation shading for relief style mapping. An example has been added
  which uses both of these functions to produce a grayscale shaded relief map:

.. plot:: examples/srtm_shading.py
   :width: 300pt

* Lion Krischer extended the capability of
  :class:`~cartopy.io.img_tiles.GoogleTiles` to allow support for **street**,
  **satellite**, **terrain** and **street_only** style Google Map tiles.

* Nat Wilson's contribution brought us a major step closer to Python 3 compatibility.

* Support for the :class:`~cartopy.crs.UTM` projection was added by Mark Hedley.

* Andrew Dawson has added a new convenience utility function
  :func:`~cartopy.util.add_cyclic_point` to add a cyclic point to an array and
  optionally to a corresponding 1D coordinate.

* Andrew Dawson added formatters for producing longitude/latitude tick labels for
  rectangular projections. The formatters are customizable and can be used to produce
  nice tick labels in a variety of styles:

.. plot:: examples/tick_labels.py
   :width: 300pt


-----------


What's new in cartopy 0.10
==========================

:Release: 0.10.0
:Date: 17 January 2014

We are very pleased to announce that Andrew Dawson was added to the cartopy
core development team. In this release Andrew has single-handedly
implemented comprehensive vector transformation and visualisation
capabilities, including: 

* The ability to transform vector fields between different coordinate
  reference systems via the :meth:`~cartopy.crs.CRS.transform_vectors`
  CRS method.

* :meth:`GeoAxes.quiver <cartopy.mpl.geoaxes.GeoAxes.quiver>` and
  :meth:`GeoAxes.barbs <cartopy.mpl.geoaxes.GeoAxes.barbs>` for arrow and
  barb plotting. More information is available at :ref:`vector_plotting`.

* A regridding function for "regularising" a vector field in the target
  coordinate system. See also
  :func:`cartopy.vector_transform.vector_scalar_to_grid`. Both
  :meth:`~cartopy.mpl.geoaxes.GeoAxes.quiver` and
  :meth:`~cartopy.mpl.geoaxes.GeoAxes.barbs` accept the ``regrid_shape``
  keyword to trigger this behaviour automatically. 
  
* :meth:`GeoAxes.streamplot <cartopy.mpl.geoaxes.GeoAxes.streamplot>` adds
  the ability to draw streamlines in any projection from a vector field in
  any other projection.

.. plot:: examples/barbs.py
    :width: 300pt


-----------


What's new in cartopy 0.9
=========================

:Release: 0.9.0
:Date: 12 September 2013

* We are very pleased to announce that Bill Little was added to the cartopy
  core development team. Bill has made some excellent contributions to cartopy,
  and `his presentation at EuroScipy'13 on
  "Iris & Cartopy" <https://www.euroscipy.org/2013/schedule/presentation/35/>`_
  was voted best talk of the conference.
* Other talks and tutorials during this release cycle include Phil Elson's `talk at SciPy'13
  (with video) <http://conference.scipy.org/scipy2013/presentation_detail.php?id=132>`_,
  `Thomas Lecocq's tutorial at EuroSciPy
  <https://www.euroscipy.org/2013/schedule/presentation/27/>`_
  and a forthcoming `talk at FOSS4G <http://2013.foss4g.org/conf/programme/presentations/29/>`_.
* Christoph Gohlke updated cartopy to support Windows 7.
* The Plate Carree projection was updated to fully handle arbitrary globe definitions.
* Peter Killick updated the Mercator class' default globe to WGS84. His refactor paved the way
  for some follow on work to fully implement the Google Spherical Mercator (EPSG:3857) projection.

    |image_eyja_volcano|_

    .. |image_eyja_volcano| image:: examples/eyja_volcano_00_00.thumb.png

    .. _image_eyja_volcano: examples/eyja_volcano.html

* The TransverseMercator class saw a tidy up to include several common arguments (:pull:`pull request <309>`)
* Bill Little added the Geostationary projection to allow geolocation of satellite imagery.
  
    |image_geostationary|_

    .. |image_geostationary| image:: examples/geostationary_00_00.thumb.png

    .. _image_geostationary: examples/geostationary.html

* Byron Blay added the :class:`Lambert conformal conic projection <cartopy.crs.LambertConformal>`.


-----------



What's new in cartopy 0.8
=========================

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


-----------



What's new in cartopy 0.7
=========================

:Release: 0.7.0
:Date: 21 Mar 2013

* Carwyn Pelley added support for 2D arrays of points to :meth:`cartopy.crs.CRS.transform_points`. (:pull:`192`)
* Phil Elson added control for the gridlines and tick labels drawn with
  :meth:`cartopy.mpl.geoaxes.GeoAxes.gridlines`. (:pull:`238`)
* Various documentation enhancements have been added. (:pull:`247`, :pull:`244` :pull:`240` and :pull:`242`)

This is a quick release which targets two very specific requirements. The goals outlined in the development plan at
``v0.6`` still remain the primary target for ``v0.8`` and beyond.



-----------


What's new in cartopy 0.6
=========================

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
--------------------------------------------

* Improve the projection definitions to support better control over datum definitions
  and consider adding WKT support (:issue:`ticket <153>`).

* Begin work on vector field support (barbs, quiver, streamlines etc.).

* Continue identifying and implementing performance enhancements (particularly in contour drawing).

* Extend the number of projections for which it is possible to draw tick marks.


-----------


What's new in cartopy 0.5
=========================

:Release: 0.5.0
:Date: 7 Dec 2012

This document explains the new/changed features of cartopy in version 0.5.

Release 0.5 of cartopy continues the work to expand the feature-set of
cartopy to encompass common operations, and provide performance
improvements.


Cartopy 0.5 features
--------------------

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
-----------

A new features api is now available, see :doc:`tutorials/using_the_shapereader`.

.. literalinclude:: /examples/features.py

.. plot:: examples/features.py
