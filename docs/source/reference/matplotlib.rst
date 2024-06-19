.. _api.mpl:

Matplotlib interface (cartopy.mpl)
----------------------------------

Cartopy extends some Matplotlib capabilities to handle geographic
projections, such as non-rectangular axes and spines.

.. automodule:: cartopy.mpl

Geoaxes
~~~~~~~

.. automodule:: cartopy.mpl.geoaxes

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    GeoAxes
    GeoAxesSubplot
    GeoSpine
    InterProjectionTransform

Gridlines and ticks
~~~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.mpl.gridliner

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    Gridliner

.. automodule:: cartopy.mpl.ticker

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    LongitudeFormatter
    LatitudeFormatter
    LongitudeLocator
    LatitudeLocator

Artist extensions
~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.mpl.feature_artist

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    FeatureArtist

.. automodule:: cartopy.mpl.slippy_image_artist

.. autosummary::
    :toctree: generated/
    :template: autosummary/class_without_inherited.rst

    SlippyImageArtist

Patch
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.mpl.patch

.. autosummary::
    :toctree: generated/

    geos_to_path
    path_segments
    path_to_geos
