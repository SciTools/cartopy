# (C) British Crown Copyright 2011 - 2012, Met Office
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
This module contains generic functionality to support Cartopy image transformations.

"""

import matplotlib.image
import numpy
import scipy.spatial

import cartopy.crs as ccrs


def stack(*one_dim_arrays):
    """
    XXX Not used external to this module. Only used within the module 
    XXX by functions that are not used themselves. DELETE ME?

    Stack the given N :class:`numpy.ndarray`'s of similar shape M into a
    single :class:`numpy.ndarray` of shape (M, N,).

    Also see, :func:`~cartopy.img_transform.unpack`.

    Args:

    * one_dim_arrays:
        One or more :class:`~numpy.ndarray` instance of the same shape.

    Returns:
        A :class:`~numpy.ndarray` instance of shape (M, N).
    
    """

    # Ensure all arrays are of the same shape.
    first = one_dim_arrays[0]
    first_shape = first.shape
    assert all([arr.shape == first_shape for arr in one_dim_arrays])

    # Concatenate all the arrays together into one stacked array.
    arrs = list(one_dim_arrays)
    for arr in arrs:
        arr.shape = first_shape + (1,)
    r = numpy.concatenate(arrs, axis= -1)

    # Undo the shape change to the original arrays.
    for arr in arrs:
        arr.shape = first_shape
    return r


def unpack(data):
    """
    XXX Not used! DELETE ME?

    Unpack the :class:`~numpy.ndarray` of shape (..., M) into M separate :class:`~numpy.ndarray`'s.

    Also see, :func:`~cartopy.img_transform.stack`.

    Args:
    
    * data:
        The :class:`~numpy.ndarray` to be unpacked.
    
    Returns:
        A list of M :class:`~numpy.ndarray`'s.

    """

    return [data[..., i] for i in xrange(data.shape[-1])]


def ll_to_cart(lonslats):
    """
    XXX Not used! DELETE ME!

    Converts longitude and latitude coordinate data into Cartesian coordinate data.

    Args:

    * lonslats:
        A :class:`~numpy.ndarray` instance of shape (M, 2) containing M coordinate 
        data values in longitude, latitude coordinate order.

    Returns:
        A single :class:`~numpy.ndarray` of shape (M, 3) containing M 
        Cartesian coordinate data values in x, y, z order.

    """

    # Unpack into longitude and latitude values.
    lons, lats = unpack(lonslats)

    # Convert to Cartesian coordinates.
    # XXX Should be done using crs.Geocentric
    x = numpy.sin(numpy.deg2rad(90 - lats)) * numpy.cos(numpy.deg2rad(lons))
    y = numpy.sin(numpy.deg2rad(90 - lats)) * numpy.sin(numpy.deg2rad(lons))
    z = numpy.cos(numpy.deg2rad(90 - lats))
    return stack(x, y, z)


def mesh_projection(projection, nx, ny, x_extents=[None, None], y_extents=[None, None]):
    """
    Returns sample points in the given projection which span the entire projection range evenly.

    The range of the x-direction and y-direction sample points will be within the bounds
    of the projection or specified extents.

    Args:
    
    * projection:
        A :class:`~cartopy.crs.Projection` instance.

    * nx:
        The number of sample points in the projection x-direction. 

    * ny:
        The number of sample points in the projection y-direction.

    Kwargs:

    * x_extents:
        The (lower, upper) x-direction extent of the projection. 
        Defaults to the :attribute:`~cartopy.crs.Projection.x_limits`.

    * y_extents:
        The (lower, upper) y-direction extent of the projection.
        Defaults to the :attribute:`~cartopy.crs.Projection.y_limits`.
    
    Returns:
        A tuple of three items. The x-direction sample points :class:`numpy.ndarray` 
        of shape (nx, ny), y-direction sample points :class:`numpy.ndarray` of shape (nx, ny),
        and the extent of the projection range as (x-lower, x-upper, y-lower, y-upper).

    """

    # Establish the x-direction and y-direction extents.
    x_lower = x_extents[0] or projection.x_limits[0]
    x_upper = x_extents[1] or projection.x_limits[1]
    y_lower = y_extents[0] or projection.y_limits[0]
    y_upper = y_extents[1] or projection.y_limits[1]

    # Calculate evenly spaced sample points spanning the extent - excluding endpoint.
    x, xstep = numpy.linspace(x_lower, x_upper, nx, retstep=True, endpoint=False)
    y, ystep = numpy.linspace(y_lower, y_upper, ny, retstep=True, endpoint=False)

    # Offset the sample points to be within the extent range.
    x += 0.5 * xstep
    y += 0.5 * ystep

    # Generate the x-direction and y-direction meshgrids.
    x, y = numpy.meshgrid(x, y)
    return x, y, [x_lower, x_upper, y_lower, y_upper]


def projection_coords(projection, nx, ny):
    """
    XXX Not used and wrong! DELETE ME!

    Returns coords in the projection which span the entire projection range evenly.
    
    The return value is (natives, latslons, xyzs) for convenience.
    
    """

    x = numpy.linspace(projection.x_limits[0], projection.x_limits[1], nx).reshape(1, -1)
    y = numpy.linspace(projection.y_limits[0], projection.y_limits[1], ny).reshape(-1, 1)
    x, y = numpy.meshgrid(x, y)
    x, y = x.flat, y.flat
    native = stack(x, y)
    lons, lats = projection.unproject_points(x, y)
    ll = stack(lons, lats)
    xyz = ll_to_cart(ll)
    return native, ll, xyz


def get_img_coords_and_nxy(fname, projection):
    """
    XXX Not used! DELETE ME!

    """

    img = matplotlib.image.imread(fname)
    ny, nx = img.shape[:-1]
    _, _, xyz = projection_coords(projection, nx, ny)

    return img, xyz, nx, ny


def warp_img(fname, target_proj, source_proj=None, target_res=(400, 200)):
    """
    Regrid the image file from the source projection to the target projection.

    Args:

    * fname:
        Image filename to be loaded and warped.

    * target_proj:
        The target :class:`~cartopy.crs.Projection` instance for the image.

    Kwargs:

    * source_proj:
        The source :class:`~cartopy.crs.Projection` instance of the image.
        Defaults to a :class:`~cartopy.crs.PlateCarree` projection.

    * target_res:
        The (nx, ny) resolution of the target projection. Where nx defaults to
        400 sample points, and ny defaults to 200 sample points.

    """

    if source_proj is None:
        source_proj = ccrs.PlateCarree()

    raise NotImplementedError('Not yet implemented.')


def warp_array(array, target_proj, source_proj=None, target_res=(400, 200), source_extent=None, target_extent=None):
    """
    Regrid the data array from the source projection to the target projection.

    Also see, :function:`~cartopy.img_transform.regrid`.

    Args:

    * array:
        The :class:`numpy.ndarray` of data to be regridded to the target projection.

    * target_proj:
        The target :class:`~cartopy.crs.Projection` instance for the data.

    Kwargs:

    * source_proj:
        The source :class:`~cartopy.crs.Projection' instance of the data.
        Defaults to a :class:`~cartopy.crs.PlateCarree` projection.

    * target_res:
        The (nx, ny) resolution of the target projection. Where nx defaults to
        400 sample points, and ny defaults to 200 sample points.

    * source_extent:
        The (x-lower, x-upper, y-lower, y-upper) extent in native
        source projection coordinates.

    * target_extent:
        The (x-lower, x-upper, y-lower, y-upper) extent in native
        target projection coordinates.

    Returns:
        A tuple of the regridded :class:`numpy.ndarray` in the target projection and
        the (x-lower, x-upper, y-lower, y-upper) target projection extent.

    """

    # source_extent is in source coordinates.
    if source_extent is None:
        source_extent = [None] * 4
    # target_extent is in target coordinates.
    if target_extent is None:
        target_extent = [None] * 4

    source_x_extents = source_extent[:2]
    source_y_extents = source_extent[2:]

    target_x_extents = target_extent[:2]
    target_y_extents = target_extent[2:]

    if source_proj is None:
        source_proj = ccrs.PlateCarree()

    ny, nx = array.shape[:2]
    source_native_xy = mesh_projection(source_proj, nx, ny,
                                       x_extents=source_x_extents, y_extents=source_y_extents)

    # XXX Take into account the extents of the original to determine target_extents? 
    target_native_x, target_native_y, extent = mesh_projection(target_proj, target_res[0], target_res[1],
                                       x_extents=target_x_extents, y_extents=target_y_extents)

    array = regrid(array, source_native_xy[0], source_native_xy[1],
                           source_proj, target_proj,
                           target_native_x, target_native_y)
    return array, extent


def regrid(array, source_x_coords, source_y_coords, source_cs, target_proj,
           target_x_points, target_y_points):
    """
    Regrid the data array from the source projection to the target projection.

    Args:
    
    * array:
        The :class:`numpy.ndarray` of data to be regridded to the target projection.

    * source_x_coords:
        A 2-dimensional source projection :class:`numpy.ndarray` of x-direction sample points. 

    * source_y_coords:
        A 2-dimensional source projection :class:`numpy.ndarray` of y-direction sample points.

    * source_cs:
        The source :class:`~cartopy.crs.Projection` instance.

    * target_cs:
        The target :class:`~cartopy.crs.Projection` instance.

    * target_x_points:
        A 2-dimensional target projection :class:`numpy.ndarray` of x-direction sample points.

    * target_y_points:
        A 2-dimensional target projection :class:`numpy.ndarray` of y-direction sample points.

    Returns:
        The data array regridded in the target projection.

    """

    # n.b. source_cs is actually a projection (the coord system of the
    # source coordinates), but not necessarily the native projection of
    # the source array (i.e. you can provide a warped image with lat lon
    # coordinates).

    #XXX NB. target_x and target_y must currently be rectangular (i.e. be a 2d np array)
    xyz = source_cs.as_geocentric().transform_points(source_cs,
                                                     source_x_coords.flatten(),
                                                     source_y_coords.flatten())
    target_xyz = source_cs.as_geocentric().transform_points(target_proj,
                                                            target_x_points.flatten(),
                                                            target_y_points.flatten())

    kdtree = scipy.spatial.cKDTree(xyz)
    distances, indices = kdtree.query(target_xyz, k=1)
    mask = numpy.isinf(distances)

    desired_ny, desired_nx = target_x_points.shape
    if array.ndim == 2:
        # Handle missing neighours using a masked array
        if numpy.any(mask):
            indices = numpy.where(numpy.logical_not(mask), indices, 0)
            array_1d = numpy.ma.array(array.reshape(-1)[indices], mask=mask)
        else:
            array_1d = array.reshape(-1)[indices]

        new_array = array_1d.reshape(desired_ny, desired_nx)
    elif array.ndim == 3:
        # Handle missing neighours using a masked array
        if numpy.any(mask):
            indices = numpy.where(numpy.logical_not(mask), indices, 0)
            array_2d = array.reshape(-1, array.shape[-1])[indices]
            mask, array_2d = numpy.broadcast_arrays(mask.reshape(-1, 1), array_2d)
            array_2d = numpy.ma.array(array_2d, mask=mask)
        else:
            array_2d = array.reshape(-1, array.shape[-1])[indices]

        new_array = array_2d.reshape(desired_ny, desired_nx, array.shape[-1])
    else:
        raise ValueError('Expected array.ndim to be 2 or 3, got {}'.format(array.ndim))

    # Do double transform to clip points that do not map back and forth
    # to the same point to within a fixed fractional offset.
    # XXX THIS ONLY NEEDS TO BE DONE FOR (PSEUDO-)CYLINDRICAL PROJECTIONS (OR ANY OTHERS
    # WHICH HAVE THE CONCEPT OF WRAPPING)
    source_desired_xyz = source_cs.transform_points(target_proj,
                                                    target_x_points.flatten(),
                                                    target_y_points.flatten())
    back_to_target_xyz = target_proj.transform_points(source_cs,
                                                      source_desired_xyz[:, 0],
                                                      source_desired_xyz[:, 1])
    back_to_target_x = back_to_target_xyz[:, 0].reshape(desired_ny, desired_nx)
    back_to_target_y = back_to_target_xyz[:, 1].reshape(desired_ny, desired_nx)
    FRACTIONAL_OFFSET_THRESHOLD = 0.1 # data has moved by 10% of the map
    
    x_extent = numpy.abs(target_proj.x_limits[1] - target_proj.x_limits[0])
    y_extent = numpy.abs(target_proj.y_limits[1] - target_proj.y_limits[0])

    non_self_inverse_points = (numpy.abs(target_x_points - back_to_target_x) /
                               x_extent) > FRACTIONAL_OFFSET_THRESHOLD
    if numpy.any(non_self_inverse_points):
        if numpy.ma.isMaskedArray(new_array):
            new_array.mask[non_self_inverse_points] = True
        else:
            new_array = numpy.ma.array(new_array, mask=True)
            new_array.mask[...] = False
            if new_array.ndim == 3:
                for i in range(new_array.shape[2]):
                    new_array.mask[:, :, i] = non_self_inverse_points
            else:
                new_array.mask[...] = non_self_inverse_points
    non_self_inverse_points = (numpy.abs(target_y_points - back_to_target_y) /
                               y_extent) > FRACTIONAL_OFFSET_THRESHOLD
    if numpy.any(non_self_inverse_points):
        if numpy.ma.isMaskedArray(new_array):
            new_array.mask[non_self_inverse_points] = True
        else:
            new_array = numpy.ma.array(new_array, mask=non_self_inverse_points)
    return new_array
