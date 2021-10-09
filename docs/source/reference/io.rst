.. _api.io:

.. currentmodule:: cartopy.io

Input/output capabilities (cartopy.io)
--------------------------------------

Cartopy has many built-in image and map acquisition capabilities. These
capabilities allow the maps to be loaded, saved, and retrieved in various
data formats.

.. _api.io.shapereader:

Shapefiles
~~~~~~~~~~

Cartopy provides a basic interface for accessing shapefiles.

.. autosummary::
    :toctree: generated/

    shapereader.Reader
    shapereader.BasicReader
    shapereader.Record
    shapereader.natural_earth
    shapereader.NEShpDownloader
    shapereader.gshhs
    shapereader.GSHHSShpDownloader


Image collections
~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    img_nest.Img
    img_nest.ImageCollection
    img_nest.NestedImageCollection


Image tiles
~~~~~~~~~~~

These classes provide an interface to the respective tile resources to
automatically load the proper tile and resolution depending on the desired domain.

.. autosummary::
    :toctree: generated/

    img_tiles.OSM
    img_tiles.GoogleTiles
    img_tiles.GoogleWTS
    img_tiles.MapQuestOSM
    img_tiles.MapQuestOpenAerial
    img_tiles.MapboxStyleTiles
    img_tiles.MapboxTiles
    img_tiles.OrdnanceSurvey
    img_tiles.QuadtreeTiles
    img_tiles.Stamen


Open Geospatial Consortium (OGC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several classes to enable interfacing with OGC clients.

.. autosummary::
    :toctree: generated/

    ogc_clients.WFSGeometrySource
    ogc_clients.WMSRasterSource
    ogc_clients.WMTSRasterSource


Shuttle Radar Topography Mission (SRTM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SRTM data can be accessed through the :mod:`cartopy.io.srtm` module
using classes and functions defined below.

.. autosummary::
    :toctree: generated/
    :recursive:

    srtm.SRTM1Source
    srtm.SRTM3Source
    srtm.SRTMDownloader
    srtm.read_SRTM
    srtm.read_SRTM1
    srtm.read_SRTM3
    srtm.add_shading


Base classes and functions
~~~~~~~~~~~~~~~~~~~~~~~~~~

These are the base classes that new resources can leverage
to implement a new reader or tile client.

.. autosummary::
    :toctree: generated/

    Downloader
    DownloadWarning
    LocatedImage
    RasterSource
    RasterSourceContainer
    PostprocessedRasterSource
    fh_getter
