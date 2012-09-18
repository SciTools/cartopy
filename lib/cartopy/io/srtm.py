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


import json
_JSON_SRTM3_LOOKUP = os.path.join(os.path.dirname(__file__), os.path.splitext(os.path.basename(__file__))[0] + '.json')
_SRTM3_FILE_LOOKUP = json.load(open(_JSON_SRTM3_LOOKUP))


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

    # xxx extent may need to be wider by half a pixel
    return elev[::-1, ...], ccrs.PlateCarree(), [x, x + 1, y, y + 1]


def SRTM3_retrieve(lon, lat, data_dir=None):
    if data_dir is None:
        dname = os.path.dirname
        # be more clever in the data directory so that users can define a setting.
        data_dir = os.path.join(dname(dname(__file__)), 'data', 'SRTM3')

    x = '%s%03d' % ('E' if lon > 0 else 'W', abs(int(lon)))
    y = '%s%02d' % ('N' if lat > 0 else 'S', abs(int(lat)))

    filename = '{y}{x}.hgt'.format(x=x, y=y)

    filepath = os.path.join(data_dir, filename)

    if not os.path.exists(filepath):
        # download the zip file
        file_url = _SRTM3_FILE_LOOKUP.get(u'{y}{x}'.format(x=x, y=y), None)
        # no file exists in the SRTM dataset for these coordinates
        if file_url is None:
            return None

        import urllib2
        import cStringIO as StringIO
        from zipfile import ZipFile

        if not os.path.exists(data_dir):
            os.makedirs(data_dir)

        srtm_online = urllib2.urlopen(file_url)
        zfh = ZipFile(StringIO.StringIO(srtm_online.read()), 'r')
        zfh.extract(filename, data_dir)

    return filepath


def _create_srtm3_dict():
    """
    Returns a dictionary mapping srtm filename to the URL of the file.

    This is slow as it must query the SRTM server to identify the continent from
    which the tile comes. Hence a json file with this content exists in ```_JSON_SRTM3_LOOKUP```.

    The json file was created with::

        $> python -c "import cartopy.io.srtm as srtm; import json; json.dump(srtm._create_srtm3_dict(), open(srtm._JSON_SRTM3_LOOKUP, 'w'));"

    """
    # lazy imports. In most situations, these are not dependencies of cartopy.
    import urllib
    from BeautifulSoup import BeautifulSoup

    files = {}

    for continent in ['Australia', 'Africa', 'Eurasia', 'Islands', 'North_America', 'South_America']:

        url = "http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/%s" % continent
        f = urllib.urlopen(url)
        html = f.read()
        soup = BeautifulSoup(html)

        for link in soup('li'):
            name = str(link.text)
            if name != ' Parent Directory':
                # remove the '.hgt.zip'
                files[name[:-8]] = url + '/' + name
    return files