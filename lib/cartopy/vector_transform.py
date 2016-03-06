# (C) British Crown Copyright 2013 - 2016, Met Office
#
# This file is part of cartopy.
#
# cartopy is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cartopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cartopy.  If not, see <https://www.gnu.org/licenses/>.
"""
This module contains generic functionality to support Cartopy vector
transforms.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from scipy.interpolate import griddata


def _interpolate_to_grid(nx, ny, x, y, *scalars, **kwargs):
    """
    Interpolates two vector components and zero or more scalar fields,
    which can be irregular, to a regular grid.

    Kwargs:

    * target_extent:
        The extent in the target CRS that the grid should occupy, in the
        form ``(x-lower, x-upper, y-lower, y-upper)``. Defaults to cover
        the full extent of the vector field.

    """
    target_extent = kwargs.get('target_extent', None)
    if target_extent is None:
        target_extent = (x.min(), x.max(), y.min(), y.max())
    points = np.array([x.ravel(), y.ravel()]).T
    x0, x1, y0, y1 = target_extent
    x_grid, y_grid = np.meshgrid(np.linspace(x0, x1, nx),
                                 np.linspace(y0, y1, ny))
    s_grid_tuple = tuple()
    for s in scalars:
        s_grid_tuple += (griddata(points, s.ravel(), (x_grid, y_grid),
                                  method='linear'),)
    return (x_grid, y_grid) + s_grid_tuple


def vector_scalar_to_grid(src_crs, target_proj, regrid_shape, x, y, u, v,
                          *scalars, **kwargs):
    """
    Transform and interpolate a vector field to a regular grid in the
    target projection.

    Args:

    * src_crs:
        The :class:`~cartopy.crs.CRS` that represents the coordinate
        system the vectors are defined in.

    * target_proj:
        The :class:`~cartopy.crs.Projection` that represents the
        projection the vectors are to be transformed to.

    * regrid_shape:
        The regular grid dimensions. If a single integer then the grid
        will have that number of points in the x and y directions. A
        2-tuple of integers specify the size of the regular grid in the
        x and y directions respectively.

    * x, y:
        The x and y coordinates, in the source CRS coordinates,
        where the vector components are located.

    * u, v:
        The grid eastward and grid northward components of the
        vector field respectively. Their shapes must match.

    * scalars:
        Zero or more scalar fields to regrid along with the vector
        components. Each scalar field must have the same shape as the
        vector components.

    Kwargs:

    * target_extent:
        The extent in the target CRS that the grid should occupy, in the
        form ``(x-lower, x-upper, y-lower, y-upper)``. Defaults to cover
        the full extent of the vector field.

    Returns:

    * x_grid, y_grid:
        The x and y coordinates of the regular grid points as
        2-dimensional arrays.

    * u_grid, v_grid:
        The eastward and northward components of the vector field on
        the regular grid.

    * scalars_grid:
        The scalar fields on the regular grid. The number of returned
        scalar fields is the same as the number that were passed in.

    """
    if u.shape != v.shape:
        raise ValueError('u and v must be the same shape')
    if x.shape != u.shape:
        x, y = np.meshgrid(x, y)
        if not (x.shape == y.shape == u.shape):
            raise ValueError('x and y coordinates are not compatible '
                             'with the shape of the vector components')
    if scalars:
        for s in scalars:
            if s.shape != u.shape:
                raise ValueError('scalar fields must have the same '
                                 'shape as the vector components')
    try:
        nx, ny = regrid_shape
    except TypeError:
        nx = ny = regrid_shape
    if target_proj != src_crs:
        # Transform the vectors to the target CRS.
        u, v = target_proj.transform_vectors(src_crs, x, y, u, v)
        # Convert Coordinates to the target CRS.
        proj_xyz = target_proj.transform_points(src_crs, x, y)
        x, y = proj_xyz[..., 0], proj_xyz[..., 1]
    # Now interpolate to a regular grid in projection space, treating each
    # component as a scalar field.
    return _interpolate_to_grid(nx, ny, x, y, u, v, *scalars, **kwargs)
