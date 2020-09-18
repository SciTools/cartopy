.. _api.transformations:

.. currentmodule:: cartopy

Image and vector transformations
--------------------------------

There are several useful utility functions in cartopy to help transform
and reshape data when going from one projection to another.

Image transformations
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

    img_transform.mesh_projection
    img_transform.regrid
    img_transform.warp_array
    img_transform.warp_img

Vector transformations
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

    vector_transform.vector_scalar_to_grid

Longitude wrapping
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

    util.add_cyclic_point

LinearRing/LineString projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   trace.project_linear
   trace.Interpolator
   trace.CartesianInterpolator
   trace.SphericalInterpolator
   trace.LineAccumulator