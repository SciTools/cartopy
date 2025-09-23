.. _cartopy_developer_interfaces:

Cartopy developer interfaces
============================

Cartopy exposes several interfaces to help make it easier to add new
functionality quickly and easily.

..
    Feature API
    -----------
    TODO

    Image API
    ---------
    TODO


Data/Download API
-----------------

In order to keep the size of a cartopy release low, the majority of data is
not included as standard. This means that, when developing new features, it is
often necessary to provide interfaces which can acquire data from external
sources (normally via HTTP). The :class:`~cartopy.io.Downloader` class
has been designed to make this process as easy as possible for developers to
extend, whilst still making it possible for users to configure in their own
way.


An example of specialising this class can be found in
:mod:`cartopy.io.shapereader.NEShpDownloader` which enables the downloading of
zipped shapefiles from the `<https://www.naturalearthdata.com>`_ website. All
known subclasses of :class:`~cartopy.io.Downloader` are listed below for
reference:

   * :class:`cartopy.io.shapereader.NEShpDownloader`
   * :class:`cartopy.io.srtm.SRTMDownloader`


Raster images
-------------

The abstraction between retrieval and visualisation of raster data means
that the :class:`cartopy.io.RasterSource` class exists to retrieve an image
(given sufficient context of projection, extent, resolution etc.) while in the
matplotlib interface the
:class:`cartopy.mpl.slippy_image_artist.SlippyImageArtist`
class feeds the appropriate information to the
:class:`~cartopy.io.RasterSource` and visualises it on a map.
The orchestration in Matplotlib is made more convenient
to the user of a :class:`~cartopy.mpl.geoaxes.GeoAxes` through the
:class:`~cartopy.mpl.geoaxes.GeoAxes.add_raster` method. Anything which exposes
the ``validate_projection`` and ``fetch_raster`` methods in the form described
in :class:`~cartopy.io.RasterSource` can be used as a slippy maps source in
this way.


.. currentmodule:: cartopy.mpl.slippy_image_artist

The :class:`SlippyImageArtist` class
provides panning and zooming of image sources which are able to
re-retrieve data (such as that from a web service) for efficient and
interactive visualisation. Generally the SlippyImageArtist is a developer's
interface, with users often creating a
:class:`SlippyImageArtist` instance through
the GeoAxes' :meth:`~cartopy.mpl.geoaxes.GeoAxes.add_raster` method.
