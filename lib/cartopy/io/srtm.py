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
The Shuttle Radar Topography Mission (SRTM) is an international research
effort that obtained digital elevation models on a near-global scale from
56S to 60N, to generate the most complete high-resolution digital topographic
database of Earth prior to the release of the ASTER GDEM in 2009.

   - Wikipedia (August 2012)

"""
import json
import os

import numpy

from cartopy import config
import cartopy.crs as ccrs
from cartopy.io import fh_getter, Downloader


def srtm(lon, lat):
    """
    Return (elevation, crs, extent) for the given longitude latitude.
    """
    fname = SRTM3_retrieve(lon, lat)
    if fname is None:
        raise ValueError('No srtm tile found for those coordinates.')
    return read_SRTM3(fname)


def srtm_composite(lon_min, lat_min, nx, ny):

    # XXX nx and ny have got confused in the code (numpy array ordering?).
    # However, the interface works well.

    bottom_left_ll = (lon_min, lat_min)
    shape = numpy.array([1201, 1201])
    img = numpy.empty(shape * (nx, ny))

    for i in range(nx):
        for j in range(ny):
            x_img_slice = slice(i * shape[0], (i + 1) * shape[0])
            y_img_slice = slice(j * shape[1], (j + 1) * shape[1])

            tile_img, crs, extent = srtm(bottom_left_ll[0] + j,
                                         bottom_left_ll[1] + i)
            img[x_img_slice, y_img_slice] = tile_img

    extent = (bottom_left_ll[0], bottom_left_ll[0] + ny,
              bottom_left_ll[1], bottom_left_ll[1] + nx)

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


def SRTM3_retrieve(lon, lat):
    x = '%s%03d' % ('E' if lon > 0 else 'W', abs(int(lon)))
    y = '%s%02d' % ('N' if lat > 0 else 'S', abs(int(lat)))

    srtm_downloader = Downloader.from_config(('SRTM', 'SRTM3'))
    return srtm_downloader.path({'config': config, 'x': x, 'y': y})


class SRTM3Downloader(Downloader):
    """
    Provides a SRTM3 download mechanism.

    """
    FORMAT_KEYS = ('config', 'x', 'y')

    _JSON_SRTM3_LOOKUP = os.path.join(os.path.dirname(__file__),
                                      'srtm.json')
    _SRTM3_LOOKUP_URL = json.load(open(_JSON_SRTM3_LOOKUP, 'r'))
    """
    The SRTM3 url lookup dictionary maps keys such as 'N43E043' to the url
    of the file to download.

    """

    def __init__(self,
                 target_path_template,
                 pre_downloaded_path_template='',
                 ):
        # adds some SRTM3 defaults to the __init__ of a Downloader
        # namely, the URl is determined on the fly using the
        # ``SRTM3Downloader._SRTM3_LOOKUP_URL`` dictionary
        Downloader.__init__(self, None,
                            target_path_template,
                            pre_downloaded_path_template)

    def url(self, format_dict):
        # override the url method, looking up the url from the
        # ``SRTM3Downloader._SRTM3_LOOKUP_URL`` dictionary
        key = u'{y}{x}'.format(**format_dict)
        url = SRTM3Downloader._SRTM3_LOOKUP_URL.get(key, None)
        return url

    def acquire_resource(self, target_path, format_dict):
        import cStringIO as StringIO
        from zipfile import ZipFile

        target_dir = os.path.dirname(target_path)
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)

        url = self.url(format_dict)

        srtm_online = self._urlopen(url)
        zfh = ZipFile(StringIO.StringIO(srtm_online.read()), 'r')

        zip_member_path = u'{y}{x}.hgt'.format(**format_dict)
        member = zfh.getinfo(zip_member_path)
        with open(target_path, 'wb') as fh:
            fh.write(zfh.open(member).read())

        srtm_online.close()
        zfh.close()

        return target_path

    @staticmethod
    def _create_srtm3_dict():
        """
        Returns a dictionary mapping srtm filename to the URL of the file.

        This is slow as it must query the SRTM server to identify the
        continent from which the tile comes. Hence a json file with this
        content exists in ``SRTM3Downloader._JSON_SRTM3_LOOKUP``.

        The json file was created with::

            import cartopy.io.srtm as srtm
            import json
            fh = open(srtm.SRTM3Downloader._JSON_SRTM3_LOOKUP, 'w')
            json.dump(srtm.SRTM3Downloader._create_srtm3_dict(), fh)

        """
        # lazy imports. In most situations, these are not
        # dependencies of cartopy.
        import urllib
        from BeautifulSoup import BeautifulSoup

        files = {}

        for continent in ['Australia', 'Africa', 'Eurasia', 'Islands',
                          'North_America', 'South_America']:

            url = "http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/%s" % continent
            f = urllib.urlopen(url)
            html = f.read()
            soup = BeautifulSoup(html)

            for link in soup('li'):
                name = str(link.text)
                if name != ' Parent Directory':
                    # remove the '.hgt.zip'
                    files[name[:-8]] = url + '/' + name
            f.close()
        return files

    @classmethod
    def default_downloader(cls):
        """
        Returns a typical downloader for this class. In general, this static
        method is used to create the default configuration in cartopy.config
        """
        target_path_template = os.path.join('{config[data_dir]}', 'SRTM',
                                            'SRTM3', '{y}{x}.hgt')
        return cls(target_path_template=target_path_template)


# add a generic SRTM downloader to the config 'downloaders' section.
config['downloaders'].setdefault(('SRTM', 'SRTM3'),
                                 SRTM3Downloader.default_downloader())
