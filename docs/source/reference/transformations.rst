.. _api.transformations:

Image and vector transformations
--------------------------------

There are several useful utility functions in cartopy to help transform
and reshape data when going from one projection to another.

Image transformations
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.img_transform

.. autosummary::
    :toctree: generated/

    mesh_projection
    regrid
    warp_array
    warp_img

Vector transformations
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.vector_transform

.. autosummary::
    :toctree: generated/

    vector_scalar_to_grid

Longitude wrapping
~~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.util

.. autosummary::
    :toctree: generated/

    add_cyclic_point
    add_cyclic


LinearRing/LineString projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: cartopy.trace

.. autosummary::
    :toctree: generated/

    project_linear
    Interpolator
    CartesianInterpolator
    SphericalInterpolator
    LineAccumulator
