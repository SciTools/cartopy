.. _api.transformations:

Image and vector transformations
--------------------------------

There are several useful utility functions in cartopy to help transform
and reshape data when going from one projection to another.

Image transformations
~~~~~~~~~~~~~~~~~~~~~

.. module:: cartopy.img_transform

:mod:`cartopy.img_transform` contains:

.. autosummary::
    :toctree: generated/

    mesh_projection
    regrid
    warp_array
    warp_img

Vector transformations
~~~~~~~~~~~~~~~~~~~~~~

.. module:: cartopy.vector_transform

:mod:`cartopy.vector_transform` contains:

.. autosummary::
    :toctree: generated/

    vector_scalar_to_grid

Longitude wrapping
~~~~~~~~~~~~~~~~~~

.. module:: cartopy.util

:mod:`cartopy.util` contains:

.. autosummary::
    :toctree: generated/

    add_cyclic_point
    add_cyclic


LinearRing/LineString projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: cartopy.trace

:mod:`cartopy.trace` contains:

.. autosummary::
    :toctree: generated/

    project_linear
    Interpolator
    CartesianInterpolator
    SphericalInterpolator
    LineAccumulator
