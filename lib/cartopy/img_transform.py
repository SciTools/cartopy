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


import matplotlib.pyplot as plt
import scipy.spatial
import numpy

import cartopy.crs as ccrs


def stack(*one_dim_arrays):
    """Stack the given arrays into an array of (m, n) where m is the len of each array and n is the number of arrays given."""
    first = one_dim_arrays[0]
    first_shape = first.shape
    assert all([arr.shape == first_shape for arr in one_dim_arrays])
    
    arrs = list(one_dim_arrays)
    for arr in arrs:
        arr.shape = first_shape + (1, )
    r =  numpy.concatenate(arrs, axis=-1)
    for arr in arrs:
        arr.shape = first_shape    
    return r


def unpack(vals):
    """opposite of stack."""
    return [vals[..., i] for i in xrange(vals.shape[-1])]


def ll_to_cart(lonslats):
    lons, lats = unpack(lonslats)
    # XXX Should be done using crs.Geocentric
    x = numpy.sin(numpy.deg2rad(90 - lats)) * numpy.cos(numpy.deg2rad(lons))
    y = numpy.sin(numpy.deg2rad(90 - lats)) * numpy.sin(numpy.deg2rad(lons))
    z = numpy.cos(numpy.deg2rad(90 - lats))
    return stack(x, y, z)


def mesh_projection(projection, nx, ny, x_extents=[None, None], y_extents=[None, None]):
    """
    Returns coords in the projection which span the entire projection range evenly.
    
    Returns: 
    
        xs, ys, extent
        
    The return value is native coordinate system.
    
    """
    x_lower = x_extents[0] or projection.x_limits[0]
    x_upper = x_extents[1] or projection.x_limits[1]
    y_lower = y_extents[0] or projection.y_limits[0]
    y_upper = y_extents[1] or projection.y_limits[1]
    
    x, xstep = numpy.linspace(x_lower, x_upper, nx, retstep=True, endpoint=False)
    y, ystep = numpy.linspace(y_lower, y_upper, ny, retstep=True, endpoint=False)
    
    x += 0.5 * xstep
    y += 0.5 * ystep
    
    x, y = numpy.meshgrid(x, y)
    return x, y, [x_lower, x_upper, y_lower, y_upper]


def projection_coords(projection, nx, ny):
    """
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
    img = plt.imread(fname)
    
    ny, nx = img.shape[:-1]
    
    _, _, xyz = projection_coords(projection, nx, ny)
    
    return img, xyz, nx, ny


def warp_img(fname, target_proj, source_proj=None, target_res=(400, 200)):
    if source_proj is None:
        source_proj = ccrs.PlateCarree()
    
    raise NotImplementedError('Not yet implemented.')
    
        
def warp_array(array, target_proj, source_proj=None, target_res=(400, 200), source_extent=None, target_extent=None):
    # source_extent is in source coordinates
    if source_extent is None:
        source_extent = [None] * 4
    # target_extent is in target coordinates
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
    
    # xxx take into account the extents of the original to determine target_extents? 
    target_native_x, target_native_y, extent = mesh_projection(target_proj, target_res[0], target_res[1], 
                                       x_extents=target_x_extents, y_extents=target_y_extents)
    
    array = regrid(array, source_native_xy[0], source_native_xy[1], 
                           source_proj, target_proj, 
                           target_native_x, target_native_y)
    return array, extent
    
                
def regrid(array, source_x_coords, source_y_coords, source_cs, target_proj, target_x_points, target_y_points):
    # n.b. source_cs is actually a projection (the coord system of the source coordinates), 
    # but not necessarily the native projection of the source array (i.e. you can provide a warped image with lat lon coordinates).
    # XXX NB. target_x and target_y must currently be rectangular (i.e. be a 2d np array)
    xyz = source_cs.as_geocentric().transform_points(source_cs, source_x_coords.flatten(), source_y_coords.flatten())
    
#    lons, lats = source_cs.unproject_points(source_x_coords.flatten(), source_y_coords.flatten())
#    ll = stack(lons, lats)
#    xyz = ll_to_cart(ll)
    
#    native_xy = stack(target_x_points.flatten(), target_y_points.flatten())
#    target_geocent = target_proj.as_geocentric().transform_points(target_proj, target_x_points.flatten(), target_y_points.flatten())
#    target_ll = stack(target_geocent[:, 0], target_geocent[:, 1])
#    ll = stack(lons, lats)
#    target_xyz = ll_to_cart(ll)
    target_xyz = source_cs.as_geocentric().transform_points(target_proj, target_x_points.flatten(), target_y_points.flatten())
    
    kdtree = scipy.spatial.cKDTree(xyz)
    
    distances, indices = kdtree.query(target_xyz, k=1)
    
    desired_ny, desired_nx = target_x_points.shape
    if array.ndim == 3:
        # XXX Handle other dimensions
        new_array = array.reshape(-1, array.shape[-1])[indices].reshape(desired_ny, desired_nx, array.shape[-1])
        # XXX Handle alpha
        missing = (0, 0, 0)
    else:
        new_array = array.reshape(-1)[indices].reshape(desired_ny, desired_nx)
        missing = 0#(0, 0, 0)

    # clip any bad points (i.e. those that lie outside of the extent of the original image)
#    source_desired_x, source_desired_y = source_cs.project_points(target_ll[:, 0], target_ll[:, 1])
#    source_desired_x = source_desired_x.reshape([desired_ny, desired_nx]) 
#    source_desired_y = source_desired_y.reshape([desired_ny, desired_nx])
#    outof_extent_points = (
#                        (source_desired_x < source_cs.x_limits[0]) | 
#                        (source_desired_x > source_cs.x_limits[1]) | 
#                        (source_desired_y < source_cs.y_limits[0]) | 
#                        (source_desired_y > source_cs.y_limits[1])
#                        )
    # XXX Assumes the image is rectilinear in its native grid. This is not necessarily true...
#    outof_extent_points = (
#                        (source_desired_x < source_x_coords.min()) | 
#                        (source_desired_x > source_x_coords.max()) | 
#                        (source_desired_y < source_y_coords.min()) | 
#                        (source_desired_y > source_y_coords.max())
#                        )
#    # XXX handle different shape mask (because of other dimensions)
##    new_array = numpy.ma.array(new_array, mask=outof_extent_points)
#    new_array[outof_extent_points] = missing
#    
#    # clip any points which do not map back to their own point (i.e. 365 might map back to 5 in degrees and so would be clipped)
 #   recalculated_coords, _ = target_proj.transform_points(source_cs.as_geocentric(), target_ll[:, 0], target_ll[:, 1])
 #   recalculated_x = recalculated_x.reshape(desired_ny, desired_nx)
 #   native_x = native_xy[:, 0].reshape(desired_ny, desired_nx)
    
#    non_self_inverse_points = numpy.abs(recalculated_x - native_x) > 1e-3
#    # XXX useful for repeated values such as those in interrupted Goodes.
##    new_array.mask[non_self_inverse_points] = True
#    new_array[non_self_inverse_points] = missing

    return new_array

    
if __name__ == '__main__':
    import cartopy
    def bluemarble():
        source_proj = ccrs.PlateCarree()
        fname = '/data/local/dataZoo/cartography/raster/blue_marble_2000_1000.jpg'
        img_origin = 'lower'
        img = plt.imread(fname)
        return img, img_origin, source_proj
    
    def bluemarble_low():
        source_proj = ccrs.PlateCarree()
        fname = '/data/local/dataZoo/cartography/raster/blue_marble_720_360.png'
        img_origin = 'lower'
        img = plt.imread(fname)
        img = img[::-1]    
        return img, img_origin, source_proj
    
    def polar_bluemarble():
        source_proj = ccrs.SouthPolarStereo()
        source_proj._max = 5e6
        fname = '/home/h02/itpe/foo.png'
        img_origin = 'lower'
        img = plt.imread(fname)
        return img, img_origin, source_proj
    
    def orca1():
        from netCDF4 import Dataset
        fname = '/data/local/dataZoo/NetCDF/ORCA1/CICE_ORCA1_CF.nc'
        rootgrp = Dataset(fname, 'r', format='NETCDF4')
        lons = rootgrp.variables['TLON'][:]
        lats = rootgrp.variables['TLAT'][:]
        sst = rootgrp.variables['sst'][:]
        sst = sst.reshape(sst.shape[1:])
        return sst, lons, lats
     
    def orca2():
        from netCDF4 import Dataset
        fname = '/data/local/dataZoo/NetCDF/ORCA2/ORCA2_1d_00010101_00010101_grid_T_0000.nc'
        rootgrp = Dataset(fname, 'r', format='NETCDF4')
        lons = rootgrp.variables['nav_lon'][:]
        lats = rootgrp.variables['nav_lat'][:]
        sst = rootgrp.variables['sosstsst'][:][0, ...]
        lats = lats.astype(numpy.float64)
        lons = lons.astype(numpy.float64)
        return sst, lons, lats
     
    def do_image_warp(target_proj):
        ######### SOURCE #########################
        arr, img_origin, source_proj = bluemarble()
#        arr, img_origin, source_proj = bluemarble_low()
        arr, img_origin, source_proj = polar_bluemarble()
        
        ######### DO IT ########################
        image, extent = warp_array(arr, target_proj, source_proj, (750, 375)
                           )
        
        return image, extent, img_origin
    
    def do_data_regrid(target_proj):
        ############ SOURCE ##########################
#        data, lons, lats = orca1()
        data, lons, lats = orca2()
        # lons and lats are equivalent to a PlatCarree CS
        source_cs = ccrs.Geodetic()
        
        ############ TARGET ##########################
        target_x, target_y, extent = mesh_projection(target_proj, 500, 500)
        
        ######### DO IT ########################
        image = regrid(data, lons, lats, source_cs, target_proj, target_x, target_y)
        return image, extent, 'upward'
    
    
    ############################## MAIN STUFF ###############################
    
    
    ######### TARGET #########################
#    target_proj = ccrs.InterruptedGoodeHomolosine()
    target_proj = ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5)
#    target_proj = ccrs.Robinson()
#    target_proj = ccrs.Orthographic(central_longitude=-90, central_latitude=45)
#    target_proj = ccrs.LambertCylindrical()
    target_proj = ccrs.PlateCarree()

    
    if False:
        ############ IMAGE WARP ###############
        image, extent, image_origin = do_image_warp(target_proj)
    else:   
        ############ DATA REGRIDDING ###############
        image, extent, image_origin = do_data_regrid(target_proj)
        
    plt.axes(projection=ccrs.PlateCarree())
    plt.gca().coastlines()
    plt.imshow(image, origin=image_origin, extent=extent)
    plt.show()
