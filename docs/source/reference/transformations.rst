.. _api.transformations:

Image and vector transformations
--------------------------------

There are several useful utility functions in cartopy to help transform
and reshape data when going from one projection to another.

.. currentmodule:: cartopy.img_transform

Image transformations
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

    mesh_projection
    regrid
    warp_array
    warp_img


.. currentmodule:: cartopy.vector_transform

Vector transformations
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

    vector_scalar_to_grid


.. currentmodule:: cartopy.util

Longitude wrapping
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

    add_cyclic_point


.. currentmodule:: cartopy.trace

LinearRing/LineString projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   project_linear
   Interpolator
   CartesianInterpolator
   SphericalInterpolator
   LineAccumulator