.. _api.mpl:

Matplotlib interface (cartopy.mpl)
----------------------------------

Cartopy extends some Matplotlib capabilities to handle geographic
projections, such as non-rectangular axes and spines.

.. currentmodule:: cartopy.mpl

Geoaxes
~~~~~~~

The most primitive extension is the :class:`cartopy.mpl.geoaxes.GeoAxes` class, which
extends a Matplotlib Axes and adds a `transform` keyword
argument to many plotting methods to enable geographic projections and boundary wrapping
to occur on the axes.

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    geoaxes.GeoAxes
    geoaxes.GeoAxesSubplot
    geoaxes.GeoSpine
    geoaxes.InterProjectionTransform


.. currentmodule:: cartopy.mpl

Gridlines and ticks
~~~~~~~~~~~~~~~~~~~

Cartopy can produce gridlines and ticks in any projection and add
them to the current geoaxes projection, providing a way to add detailed
location information to the plots.

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst
    
    gridliner.Gridliner
    ticker.LongitudeFormatter
    ticker.LatitudeFormatter
    ticker.LongitudeLocator
    ticker.LatitudeLocator

Artist extensions
~~~~~~~~~~~~~~~~~

Features and images can be added to a :class:`cartopy.mpl.geoaxes.GeoAxes` through
an extension of the Matplotlib Artist interfaces.

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    feature_artist.FeatureArtist
    slippy_image_artist.SlippyImageArtist

.. currentmodule:: cartopy.mpl

Additional extensions
~~~~~~~~~~~~~~~~~~~~~

Extra functionality that is primarily intended for developers. They describe
some of the capabilities for transforming
between GEOS, Shapely, and Matplotlib paths.

.. autosummary::
    :toctree: generated/

    patch.geos_to_path
    patch.path_segments
    patch.path_to_geos