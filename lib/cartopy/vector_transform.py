# (C) British Crown Copyright 2013, Met Office
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
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.
"""
This module contains generic functionality to support Cartopy vector
transforms.

"""
import numpy as np
from scipy.interpolate import griddata


def _interpolate_to_grid(x, y, s, nx, ny, target_extent=None):
    """Linear interpolation of possibly irregular data to a grid."""
    if target_extent is None:
        target_extent = (x.min(), x.max(), y.min(), y.max())
    points = np.array([x.ravel(), y.ravel()]).T
    values = s.ravel()
    x0, x1, y0, y1 = target_extent
    x_grid, y_grid = np.meshgrid(np.linspace(x0, x1, nx),
                                 np.linspace(y0, y1, ny))
    s_grid = griddata(points, values, (x_grid, y_grid), method='linear')
    return x_grid, y_grid, s_grid


def scalar_to_grid(src_crs, target_crs, x, y, s, regrid_shape,
                   target_extent=None):
    """
    Transform and interpolate a scalar field to a regular grid in the
    target projection.

    Args:

    * src_crs:
        The :class:`~cartopy.crs.CRS` that represents the coordinate
        system the scalar field is defined in.

    * target_crs:
        The :class:`~cartopy.crs.CRS` that represents the coordinate
        system the scalar field is to be transformed to.

    * x, y:
        The x and y coordinates, in the source CRS coordinates,
        where the scalar field points are located.

    * s:
        The scalar field.

    * nx, ny:
        The dimensions of the desired grid in the x and y
        directions.

    Kwargs:

    * target_extent:
        The extent in the target CRS that the grid should occupy, in the
        form ``(x-lower, x-upper, y-lower, y-upper)``. Defaults to cover
        the full extent of the scalar field.

    Returns:

    * x_grid, y_grid:
        The x and y coordinates of the regular grid points as
        2-dimensional arrays.

    * s_grid:
        The scalar field on the regular grid.

    """
    if x.shape != s.shape and y.shape != s.shape:
        x, y = np.meshgrid(x, y)
        if not (x.shape == y.shape == s.shape):
            raise ValueError('x and y coordinates are not compatible '
                             'with the shape of the scalar field')
    try:
        nx, ny = regrid_shape
    except TypeError:
        nx = ny = regrid_shape
    if target_crs != src_crs:
        # Convert coordinates to the target CRS.
        proj_xyz = target_crs.transform_points(src_crs, x, y)
        x, y = proj_xyz[..., 0], proj_xyz[..., 1]
    # Now interpolate to a regular grid in the target projection space.
    return _interpolate_to_grid(x, y, s, nx, ny, target_extent=target_extent)


def vector_to_grid(src_crs, target_crs, x, y, u, v, regrid_shape,
                   target_extent=None):
    """
    Transform and interpolate a vector field to a regular grid in the
    target projection.

    Args:

    * src_crs:
        The :class:`~cartopy.crs.CRS` that represents the coordinate
        system the vectors are defined in.

    * target_crs:
        The :class:`~cartopy.crs.CRS` that represents the coordinate
        system the vectors are to be transformed to.

    * x, y:
        The x and y coordinates, in the source CRS coordinates,
        where the vector components are located.

    * u, v:
        The grid eastward and grid northward components of the
        vector field respectively. Their shapes must match.

    * nx, ny:
        The dimensions of the desired grid in the x and y
        directions.

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

    """
    if u.shape != v.shape:
        raise ValueError('u and v must be the same shape')
    if x.shape != u.shape:
        x, y = np.meshgrid(x, y)
        if not (x.shape == y.shape == u.shape):
            raise ValueError('x and y coordinates are not compatible '
                             'with the shape of the vector components')
    try:
        nx, ny = regrid_shape
    except TypeError:
        nx = ny = regrid_shape
    if target_crs != src_crs:
        # Transform the vectors to the target CRS.
        u, v = target_crs.transform_vectors(src_crs, x, y, u, v)
        # Convert Coordinates to the target CRS.
        proj_xyz = target_crs.transform_points(src_crs, x, y)
        x, y = proj_xyz[..., 0], proj_xyz[..., 1]
    # Now interpolate to a regular grid in projection space, treating each
    # component as a scalar field.
    x_grid, y_grid, u_grid = _interpolate_to_grid(x, y, u, nx, ny,
                                                  target_extent=target_extent)
    x_grid, y_grid, v_grid = _interpolate_to_grid(x, y, v, nx, ny,
                                                  target_extent=target_extent)
    return x_grid, y_grid, u_grid, v_grid
