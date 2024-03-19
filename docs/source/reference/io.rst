.. _api.io:

Input/output capabilities (cartopy.io)
--------------------------------------

Cartopy has many built-in image and map acquisition capabilities. These
capabilities allow the maps to be loaded, saved, and retrieved in various
data formats.

.. _api.io.shapereader:

Shapereader
~~~~~~~~~~~

.. automodule:: cartopy.io.shapereader

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

.. automodule:: cartopy.io.img_nest

.. autosummary::
    :toctree: generated/

    Img
    ImageCollection
    NestedImageCollection

Image tiles
~~~~~~~~~~~

.. automodule:: cartopy.io.img_tiles

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
    StadiaMapsTiles
    Stamen

Open Geospatial Consortium (OGC) Clients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.io.ogc_clients

.. autosummary::
    :toctree: generated/

    WFSGeometrySource
    WMSRasterSource
    WMTSRasterSource

Shuttle Radar Topography Mission (SRTM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.io.srtm

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

.. automodule:: cartopy.io

.. autosummary::
    :toctree: generated/

    Downloader
    DownloadWarning
    LocatedImage
    RasterSource
    RasterSourceContainer
    PostprocessedRasterSource
    fh_getter
