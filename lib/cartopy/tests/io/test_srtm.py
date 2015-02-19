# (C) British Crown Copyright 2011 - 2014, Met Office
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

from __future__ import (absolute_import, division, print_function)

import unittest
import warnings

import numpy as np
from numpy.testing import assert_array_equal

import cartopy.crs as ccrs
import cartopy.io.srtm
from cartopy.tests.io.test_downloaders import download_to_temp


def test_srtm3_retrieve():
    # test that the download mechanism for srtm3 works
    with download_to_temp() as tmp_dir:
        with warnings.catch_warnings(record=True) as w:
            r = cartopy.io.srtm.SRTM3Source().srtm_fname(-4, 50)
            assert len(w) == 1
            assert issubclass(w[0].category, cartopy.io.DownloadWarning)

        assert r.startswith(tmp_dir), 'File not downloaded to tmp dir'

        img, _, _ = cartopy.io.srtm.read_SRTM3(r)

        # check that the data is fairly sensible
        msg = 'srtm data has changed. arbitrary value testing failed.'
        assert img.max() == 602, msg
        assert img.min() == -32768, msg
        assert img[-10, 12] == 78, msg + 'Got {}'.format(img[-10, 12])


def test_srtm3_out_of_range():
    # Somewhere over the pacific the elevation should be 0.
    img, _, _ = cartopy.io.srtm.SRTM3Source().combined(120, 2, 2, 2)
    assert_array_equal(img, np.zeros(np.array((1201, 1201)) * 2))


class TestSRTM3Source__single_tile(unittest.TestCase):
    def test_out_of_range(self):
        source = cartopy.io.srtm.SRTM3Source()
        msg = 'No srtm tile found for those coordinates.'
        with self.assertRaisesRegexp(ValueError, msg):
            source.single_tile(-25, 50)

    def test_in_range(self):
        source = cartopy.io.srtm.SRTM3Source()
        img, crs, extent = source.single_tile(-1, 50)
        self.assertIsInstance(img, np.ndarray)
        self.assertEqual(img.shape, (1201, 1201))
        self.assertEqual(img.dtype, np.dtype('>i2'))
        self.assertEqual(crs, ccrs.PlateCarree())
        self.assertEqual(extent, (-1, 0, 50, 51))

    def test_zeros(self):
        source = cartopy.io.srtm.SRTM3Source()
        _, _, extent = source.single_tile(0, 50)
        self.assertEqual(extent, (0, 1, 50, 51))


class TestSRTM3Source__combined(unittest.TestCase):
    def test_trivial(self):
        source = cartopy.io.srtm.SRTM3Source()
        e_img, e_crs, e_extent = source.single_tile(-3, 50)
        r_img, r_crs, r_extent = source.combined(-3, 50, 1, 1)
        assert_array_equal(e_img, r_img)
        self.assertEqual(e_crs, r_crs)
        self.assertEqual(e_extent, r_extent)

    def test_2by2(self):
        source = cartopy.io.srtm.SRTM3Source()
        e_img, _, e_extent = source.combined(-1, 50, 2, 1)
        self.assertEqual(e_extent, (-1, 1, 50, 51))
        imgs = [source.single_tile(-1, 50)[0],
                source.single_tile(0, 50)[0]]
        assert_array_equal(np.hstack(imgs), e_img)


class TestSRTM3Source_fetch_raster(unittest.TestCase):
    def test_as_combined(self):
        source = cartopy.io.srtm.SRTM3Source()
        e_img, e_crs, e_extent = source.combined(-1, 50, 2, 1)
        imgs = source.fetch_raster(ccrs.PlateCarree(),
                                   (-0.9, 0.1, 50.1, 50.999),
                                   None)
        self.assertEqual(len(imgs), 1)
        r_img, r_extent = imgs[0]
        self.assertEqual(e_extent, r_extent)
        assert_array_equal(e_img[::-1, :], r_img)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
