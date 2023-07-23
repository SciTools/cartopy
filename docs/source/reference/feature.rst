.. _api.feature:

Feature interface (cartopy.feature)
-----------------------------------

.. module:: cartopy.feature

The feature interface can be used and extended to add various "features"
to geoaxes, such as Shapely objects and Natural Earth Imagery. The default
zorder for Cartopy features is 1.5, which puts them above images and patches,
but below lines and text.

Feature attributes
~~~~~~~~~~~~~~~~~~

Some commonly used features are stored as attributes in the
feature module for ease of use.

.. autosummary::
   :toctree: generated/

    COLORS
    auto_scaler
    BORDERS
    STATES
    COASTLINE
    LAKES
    LAND
    OCEAN
    RIVERS

Feature classes
~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

    AdaptiveScaler
    Feature
    GSHHSFeature
    NaturalEarthFeature
    Scaler
    ShapelyFeature
    WFSFeature
    nightshade.Nightshade
