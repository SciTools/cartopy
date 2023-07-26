.. _api.io:

Input/output capabilities (cartopy.io)
--------------------------------------

Cartopy has many built-in image and map acquisition capabilities. These
capabilities allow the maps to be loaded, saved, and retrieved in various
data formats.

.. _api.io.shapereader:

Shapefiles
~~~~~~~~~~

.. module:: cartopy.io.shapereader

:mod:`cartopy.io.shapereader` provides a basic interface for accessing shapefiles.

.. autosummary::
    :toctree: generated/

    Reader
    BasicReader
    FionaReader
    Record
    FionaRecord
    natural_earth
    NEShpDownloader
    gshhs
    GSHHSShpDownloader

Image collections
~~~~~~~~~~~~~~~~~

.. module:: cartopy.io.img_nest

:mod:`cartopy.io.img_nest` provides an interface for representing images.

.. autosummary::
    :toctree: generated/

    Img
    ImageCollection
    NestedImageCollection

Image tiles
~~~~~~~~~~~

.. module:: cartopy.io.img_tiles

Classes in :mod:`cartopy.io.img_tiles` provide an interface to the respective tile resources to
automatically load the proper tile and resolution depending on the desired domain.

.. autosummary::
    :toctree: generated/

    OSM
    GoogleTiles
    GoogleWTS
    MapQuestOSM
    MapQuestOpenAerial
    MapboxStyleTiles
    MapboxTiles
    OrdnanceSurvey
    QuadtreeTiles
    Stamen

Open Geospatial Consortium (OGC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: cartopy.io.ogc_clients

:mod:`cartopy.io.ogc_clients` contains several classes to enable interfacing with OGC clients.

.. autosummary::
    :toctree: generated/

    WFSGeometrySource
    WMSRasterSource
    WMTSRasterSource

Shuttle Radar Topography Mission (SRTM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: cartopy.io.srtm

The SRTM data can be accessed through the :mod:`cartopy.io.srtm` module
using classes and functions defined below.

.. autosummary::
    :toctree: generated/
    :recursive:

    SRTM1Source
    SRTM3Source
    SRTMDownloader
    read_SRTM
    read_SRTM1
    read_SRTM3
    add_shading

Base classes and functions
~~~~~~~~~~~~~~~~~~~~~~~~~~

These are the base classes in :mod:`cartopy.io` that new resources can leverage
to implement a new reader or tile client.

.. module:: cartopy.io

.. autosummary::
    :toctree: generated/

    Downloader
    DownloadWarning
    LocatedImage
    RasterSource
    RasterSourceContainer
    PostprocessedRasterSource
    fh_getter
