.. _api.io:

Input/output capabilities
-------------------------

Cartopy has many built-in image and map acquisition capabilities. These
capabilities allow the maps to be loaded, saved, and retrieved in various
data formats.

.. _api.io.shapereader:

.. currentmodule:: cartopy.io.shapereader

Shapefiles
~~~~~~~~~~

Cartopy provides a basic interface for accessing shapefiles.

.. autosummary::
    :toctree: generated/

    Reader
    BasicReader
    Record
    natural_earth
    NEShpDownloader
    gshhs
    GSHHSShpDownloader



.. currentmodule:: cartopy.io.img_nest

Image collections
~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: generated/

    Img
    ImageCollection
    NestedImageCollection


.. currentmodule:: cartopy.io.img_tiles

Image tiles
~~~~~~~~~~~

These classes provide an interface to the respective tile resources to
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
    StamenTerrain


.. currentmodule:: cartopy.io.ogc_clients

Open Geospatial Consortium (OGC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several classes to enable interfacing with OGC clients.

.. autosummary::
    :toctree: generated/

    WFSGeometrySource
    WMSRasterSource
    WMTSRasterSource


.. currentmodule:: cartopy.io.srtm

Shuttle Radar Topography Mission (SRTM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SRTM data can be accessed through the :mod:`cartopy.io.srtm` module
using classes and functions defined below.

.. autosummary::
    :toctree: generated/

    SRTM1Source
    SRTM3Source
    SRTMDownloader
    SRTM3_retrieve
    read_SRTM
    read_SRTM1
    read_SRTM3
    srtm
    srtm_composite
    add_shading
    fill_gaps


.. currentmodule:: cartopy.io

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
