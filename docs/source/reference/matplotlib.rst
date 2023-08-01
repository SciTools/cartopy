.. _api.mpl:

Matplotlib interface (cartopy.mpl)
----------------------------------

Cartopy extends some Matplotlib capabilities to handle geographic
projections, such as non-rectangular axes and spines.

.. module:: cartopy.mpl

Geoaxes
~~~~~~~

.. module:: cartopy.mpl.geoaxes

The most primitive extension is the :class:`cartopy.mpl.geoaxes.GeoAxes` class, which
extends a Matplotlib Axes and adds a `transform` keyword
argument to many plotting methods to enable geographic projections and boundary wrapping
to occur on the axes.

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    GeoAxes
    GeoAxesSubplot
    GeoSpine
    InterProjectionTransform


Gridlines and ticks
~~~~~~~~~~~~~~~~~~~

.. module:: cartopy.mpl.gridliner

Cartopy can produce gridlines and ticks in any projection and add
them to the current geoaxes projection, providing a way to add detailed
location information to the plots.

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    Gridliner

.. module:: cartopy.mpl.ticker

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    LongitudeFormatter
    LatitudeFormatter
    LongitudeLocator
    LatitudeLocator

Artist extensions
~~~~~~~~~~~~~~~~~

.. module:: cartopy.mpl.feature_artist

Features and images can be added to a :class:`cartopy.mpl.geoaxes.GeoAxes` through
an extension of the Matplotlib Artist interfaces.

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    FeatureArtist

.. module:: cartopy.mpl.slippy_image_artist

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    SlippyImageArtist


Additional extensions
~~~~~~~~~~~~~~~~~~~~~

.. module:: cartopy.mpl.patch

Extra functionality that is primarily intended for developers. They describe
some of the capabilities for transforming between Shapely, and Matplotlib paths.

.. autosummary::
    :toctree: generated/

    geos_to_path
    path_segments
    path_to_geos
