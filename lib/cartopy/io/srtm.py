"""
The Shuttle Radar Topography Mission (SRTM) is an international research effort that
obtained digital elevation models on a near-global scale from 56S to 60N, to
generate the most complete high-resolution digital topographic database of Earth prior
to the release of the ASTER GDEM in 2009.

   - Wikipedia (August 2012)

"""
import os

import numpy

import cartopy.crs as ccrs
from cartopy.io import fh_getter


def srtm(lon, lat):
    """
    Return (elevation, crs, extent) for the given longitude latitude.
    """
    fname = SRTM3_retrieve(lon, lat)
    return read_SRTM3(fname)

def srtm_composite(lon_min, lat_min, nx, ny):

    # XXX nx and ny have got confused in the code (numpy array ordering?). However, the interface works well.

    bottom_left_ll = (lon_min, lat_min)
    shape = numpy.array([1201, 1201])
    img = numpy.empty(shape * (nx, ny))

    for i in range(nx):
        for j in range(ny):
            x_img_slice = slice(i * shape[0], (i + 1) * shape[0])
            y_img_slice = slice(j * shape[1], (j + 1) * shape[1])

            tile_img, crs, extent = srtm(bottom_left_ll[0] + j, bottom_left_ll[1] + i)

#            print x_img_slice, y_img_slice
#            print img[y_img_slice, x_img_slice].shape
#            print extent
#            print

            img[x_img_slice, y_img_slice] = tile_img

    extent = (bottom_left_ll[0], bottom_left_ll[0] + ny, bottom_left_ll[1], bottom_left_ll[1] + nx)

    return img, crs, extent


def read_SRTM3(fh):
    fh, fname = fh_getter(fh, needs_filename=True)
    if fname.endswith('.zip'):
        from zipfile import ZipFile
        zfh = ZipFile(fh, 'r')
        fh = zfh.open(os.path.basename(fname[:-4]), 'r')

    elev = numpy.fromfile(fh, dtype=numpy.dtype('>i2'))
    elev.shape = (1201, 1201)

    fname = os.path.basename(fname)
    y_dir, y, x_dir, x = fname[0], int(fname[1:3]), fname[3], int(fname[4:7])


    if y_dir == 'S':
        y *= -1

    if x_dir == 'W':
        x *= -1

    return elev[::-1, ...], ccrs.PlateCarree(), [x, x + 1, y, y + 1]


def SRTM3_retrieve(lon, lat, data_dir=None):
    if data_dir is None:
        dname = os.path.dirname
        # be more clever in the data directory so that users can define a setting.
        data_dir = os.path.join(dname(dname(__file__)), 'data', 'SRTM3')

    x = '%s%03d' % ('E' if lon > 0 else 'W', abs(int(lon)))
    y = '%s%02d' % ('N' if lat > 0 else 'S', abs(int(lat)))
    # XXX figure out how to get the correct SRTM continent...
    SRTM_server = 'http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Eurasia/'

    filename = '{y}{x}.hgt'.format(base_url=SRTM_server, x=x, y=y)

    filepath = os.path.join(data_dir, filename)

    if not os.path.exists(filepath):
        # download the zip file
        import urllib2
        import cStringIO as StringIO
        from zipfile import ZipFile

        if not os.path.exists(data_dir):
            os.makedirs(data_dir)

        file_url = SRTM_server + filename + '.zip'

        srtm_online = urllib2.urlopen(file_url)
        zfh = ZipFile(StringIO.StringIO(srtm_online.read()), 'r')
        zfh.extract(filename, data_dir)

    return filepath


if __name__ == '__main__':
    fname = SRTM3_retrieve(-4, 52)
    img, crs, extent = read_SRTM3(fname)
