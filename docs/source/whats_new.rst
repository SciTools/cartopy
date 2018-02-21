What's New in cartopy 0.16
==========================

:Release: 0.16.0
:Date: 21st Feb 2018

Features
--------

* We are very pleased to announce that Ryan May has been added to the cartopy
  core development team. Ryan (@dopplershift) brings a wealth of experience,
  and has already made significant contributions to the Matplotlib interface,
  extended projections, and helped modernise the development infrastructure.

* The :ref:`gnomonic-projection` was brought up-to-date to include the
  ``central_longitude`` argument. (:pull:`855`)

* Ryan May improved the formulation of the boundary ellipse for the
  :class:`~cartopy.crs.Geostationary` projection and added the
  ``sweep_angle_axis`` keyword argument. (:pull:`890`, :pull:`897`) 

* Elliott Sales de Andrade made a number of micro-optimisations to the
  Matplotlib interface, fixed a number of documentation issues with
  Python 3 and added Matplotlib 2.0 & 2.1 compatibility. (:pull:`886`,
  :pull:`901`, :pull:`780`, :pull:`773`, :pull:`977`)

* Tick padding was added to the gridliner.
  :data:`cartopy.mpl.gridliner.Gridliner.xpadding` and
  :data:`~cartopy.mpl.gridliner.Gridliner.ypadding` relate. (:pull:`783`)

* Ryan May added the :meth:`~cartopy.feature.NaturalEarthFeature.with_scale`
  method to the NaturalEarthFeature class.
  For example, it is now possible to access higher resolution land features
  with ``cartopy.feature.LAND.with_scale('50m')``. In addition to this,
  :data:`cartopy.feature.STATES` was added to easily access administrative
  area boundaries that were previously only accessible by manually
  constructing :class:`~cartopy.feature.NaturalEarthFeature` instances
  (as is done in the :ref:`sphx_glr_gallery_feature_creation.py` example).
  (:pull:`898`)

* Daryl Herzmann and Robert Redl improved cartopy's internal conversion
  between shapely objects and Matplotlib Paths. (:pull:`885` & :pull:`1021`)

* Åsmund Steen Skjæveland fixed :meth:`cartopy.mpl.geoaxes.GeoAxes.tissot`
  to use the documented units of kilometres, where before it had been using
  metres. (:pull:`904`)

* Andrew Dawson wrote a new tutorial for the user guide:
  :ref:`understanding_transform`. (:pull:`914`)

.. figure:: _images/understanding_transform-6.png
   :target: tutorials/understanding_transform.html
   :align: center

* Daniel Kirkham and Daryl Herzmann made significant improvements to the
  stability of polygon transformation. The changes reduce the frequency
  of messages such as
  ``Unidentified problem with geometry, linestring being re-added`` and
  ``Self-intersection at or near point <X> <Y>`` occurring.
  (:pull:`974` and :pull:`903`)

* Chris Holdgraf and Corinne Bosley worked collaboratively to bring
  `sphinx-gallery <https://github.com/sphinx-gallery/sphinx-gallery>`_ to the
  cartopy docs. (:pull:`969`)

* Ray Bell neatened up many of the examples to explicitly pass the coordinate
  system when calling :meth:`~cartopy.mpl.geoaxes.GeoAxes.set_extent`.
  (:pull:`975`)

* Ryan May changed the default zorder of LAND and OCEAN to -1, thus fixing
  an issue with LAND/OCEAN appearing above some data elements such as
  vectors. (:pull:`916`)

* Kevin Donkers added the 60 UTM projections example to the gallery
  in :pull:`954`:

.. figure:: gallery/images/sphx_glr_utm_all_zones_001.png
   :target: gallery/utm_all_zones.html
   :align: center

* Andrey Kiselev added support for reading shapes with a third (Z) dimension.
  (:pull:`958`)

* Corinne Bosley standardised the docstring format for improved readability
  and visual consistency. (:pull:`987`)

* Cartopy now no longer enables :func:`shapely.speedups` at cartopy import
  time. (:pull:`990`)

* Mahé Perrette and Ryan May collaborated to improve the
  :class:`~cartopy.crs.Stereographic` projection. (:pull:`929`)



What's New in cartopy 0.15
==========================

:Release: 0.15.0
:Date: 1st February 2017

Features
--------

* The :class:`cartopy.crs.Mercator` class now allows a ``latitude_true_scale``
  to be specified.

* A ``tiles`` url can now be passed directly to the
  :class:`cartopy.io.img_tiles.GoogleTiles` class.

* The :meth:`~cartopy.mpl.geoaxes.GeoAxes.background_img` method has been
  added. This allows users to add a background image to the map, from a
  selection of pre-prepared images held in a directory specified by the
  CARTOPY_USER_BACKGROUNDS environment variable.

* The Web Map Tile Service (WMTS) interface has been extended so that WMTS
  layers can be added to geoaxes in different projections.

* The :class:`~cartopy.crs.NearsidePerspective` projection has been added.

* Optional kwargs can now be supplied to the
  :meth:`~cartopy.mpl.geoaxes.GeoAxes.add_wmts` method, which will be passed to
  the OGC WMTS ``gettile`` method.

* New additions to the gallery:

.. figure:: gallery/images/sphx_glr_axes_grid_basic_001.png
   :target: gallery/axes_grid_basic.html
   :align: center
   :scale: 70

.. figure:: gallery/images/sphx_glr_reprojected_wmts_001.png
   :target: gallery/reprojected_wmts.html
   :align: center
   :scale: 70

.. figure:: gallery/images/sphx_glr_wmts_time_001.png
   :target: gallery/wmts_time.html
   :align: center
   :scale: 70

-----------


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

  .. figure:: gallery/images/sphx_glr_tissot_001.png
     :target: gallery/tissot.html
     :align: center
     :scale: 70

* The SRTM3 data source has been changed to the `LP DAAC Data Pool
  <https://lpdaac.usgs.gov/data_access/data_pool>`_. The Data Pool is more
  consistent, fixing several missing tiles, and the data is void-filled.
  Consequently, the :func:`cartopy.srtm.fill_gaps` function has been deprecated
  as it has no purpose within the STRM context. The
  SRTM example has also been updated to skip the void-filling step.
  Additionally, this data source provides SRTM at a higher resolution of
  1 arc-second, which may be accessed via :class:`cartopy.io.srtm.SRTM1Source`.

* All downloaders will use secure connections where available. Not
  every service supports this method, and so those will use non-secured
  HTTP connections instead. (See :pull:`736` for full details.)

* Cartopy now supports, and is tested against, Matplotlib 1.3 and 1.5 as well as
  numpy 1.7, 1.8 and 1.10.

* Daniel Eriksson added a new example to the gallery:

  .. figure:: gallery/images/sphx_glr_aurora_forecast_001.png
     :target: gallery/aurora_forecast.html
     :align: center
     :scale: 70


Incompatible changes
--------------------
* :meth:`cartopy.crs.CRS.transform_point` now issues NaNs when invalid transforms are identified.


Deprecations
------------
* :data:`cartopy.crs.GOOGLE_MERCATOR` has been moved to :data:`cartopy.crs.Mercator.GOOGLE`.


-----------



What's new in cartopy 0.13
==========================

:Release: 0.13.0
:Date: 30th June 2015

Features
--------

* Andrea Smith fixed the cartopy CRS class such that 3d transforms such as :class:`cartopy.crs.Geocentric`
  now correctly apply deg2rad and rad2deg. (:pull:`625`)

* Peter Killick fixed the cartopy.crs.Mercator projection for non-zero central longitudes. (:pull:`633`)

* Conversion between Matplotlib :class:`matplotlib.path.Path` and
  :class:`shapely.geometry.Geometry <Shapely geometry>` using
  :func:`cartopy.mpl.patch.path_to_geos` and :func:`cartopy.mpl.patch.geos_to_path` now
  handles degenerate point paths.

* Update of tools/feature_download.py to allow mass download of feature data rather than
  on-demand downloading.

* A new example was added to the gallery:

  .. figure:: gallery/images/sphx_glr_eccentric_ellipse_001.png
     :target: gallery/eccentric_ellipse.html
     :align: center
     :scale: 70


-----------



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

  .. figure:: gallery/images/sphx_glr_wms_001.png
     :target: gallery/wms.html
     :align: center
     :scale: 70

* Peter Killick added an interface for accessing MapBox tiles using the MapBox
  Developer API. A MapBox client can be created with,
  :class:`~cartopy.io.img_tiles.MapboxTiles` and as with the other imagery from a simple URL
  based imagery service, it can be added to a :class:`~cartopy.mpl.geoaxes.GeoAxes` with the
  :meth:`~cartopy.mpl.geoaxes.GeoAxes.add_image` method. The following example demonstrates the
  interface for another source of imagery:

  .. figure:: gallery/images/sphx_glr_image_tiles_001.png
     :target: gallery/image_tiles.html
     :align: center
     :scale: 70

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
        fig = plt.figure(figsize=(10, 5))
        ax = plt.axes(projection=rpole)
        ax.set_global()
        ax.gridlines()
        ax.stock_img()
        ax.coastlines()
        plt.show()

* A new method has been added to the :class:`~cartopy.mpl.geoaxes.GeoAxes` to
  allow control of the neatline of a map drawn with the Matplotlib interface.
  The method, :meth:`~cartopy.mpl.geoaxes.GeoAxes.set_boundary`, takes a
  :class:`matplotlib Path<matplotlib.path.Path>` object, which means that
  arbitrary shaped edges can be achieved:

  .. figure:: gallery/images/sphx_glr_star_shaped_boundary_001.png
     :target: gallery/star_shaped_boundary.html
     :align: center
     :scale: 70

* A new SRTM3 RasterSource has been implemented allowing interactive pan/zoom
  of 3 arc-second elevation data from the Shuttle Radar Topography Mission.
  The SRTM example has also been updated to use the new interface.

* New additions to the gallery:


  .. figure:: gallery/images/sphx_glr_un_flag_001.png
     :target: gallery/un_flag.html
     :align: center
     :scale: 70

  .. figure:: gallery/images/sphx_glr_always_circular_stereo_001.png
     :target: gallery/always_circular_stereo.html
     :align: center
     :scale: 70

  .. figure:: gallery/images/sphx_glr_tube_stations_001.png
     :target: gallery/tube_stations.html
     :align: center
     :scale: 70

  .. figure:: gallery/images/sphx_glr_wms_001.png
     :target: gallery/wms.html
     :align: center
     :scale: 70

  .. figure:: gallery/images/sphx_glr_image_tiles_001.png
     :target: gallery/image_tiles.html
     :align: center
     :scale: 70


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

  .. figure:: gallery/images/sphx_glr_wmts_001.png
     :target: gallery/wmts.html
     :align: center
     :scale: 70

* Thomas Lecocq added functionality to :mod:`cartopy.io.srtm` allowing
  intelligent filling of missing elevation data, as well as a function to
  compute elevation shading for relief style mapping. An example has been added
  which uses both of these functions to produce a grayscale shaded relief map

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

  .. figure:: gallery/images/sphx_glr_tick_labels_001.png
     :target: gallery/tick_labels.html
     :align: center
     :scale: 70


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

  .. figure:: gallery/images/sphx_glr_barbs_001.png
     :target: gallery/barbs.html
     :align: center
     :scale: 70

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
  (with video) <https://conference.scipy.org/scipy2013/presentation_detail.php?id=132>`_,
  `Thomas Lecocq's tutorial at EuroSciPy
  <https://www.euroscipy.org/2013/schedule/presentation/27/>`_
  and a forthcoming `talk at FOSS4G <http://2013.foss4g.org/conf/programme/presentations/29/>`_.
* Christoph Gohlke updated cartopy to support Windows 7.
* The Plate Carree projection was updated to fully handle arbitrary globe definitions.
* Peter Killick updated the Mercator class' default globe to WGS84. His refactor paved the way
  for some follow on work to fully implement the Google Spherical Mercator (EPSG:3857) projection.


    .. figure:: gallery/images/sphx_glr_eyja_volcano_001.png
       :target: gallery/eyja_volcano.html
       :align: center
       :scale: 70

* The TransverseMercator class saw a tidy up to include several common arguments (:pull:`pull request <309>`)
* Bill Little added the Geostationary projection to allow geolocation of satellite imagery.

  .. figure:: gallery/images/sphx_glr_geostationary_001.png
     :target: gallery/geostationary.html
     :align: center
     :scale: 70

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

* Ian Edwards :doc:`added a new example <gallery/favicon>` to create a favicon for cartopy.

* Phil Elson :doc:`added a new example <gallery/hurricane_katrina>` to show polygon analysis
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

.. figure:: gallery/images/sphx_glr_features_001.png
   :target: gallery/features.html
   :align: center
   :scale: 70
