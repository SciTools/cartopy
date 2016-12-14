# (C) British Crown Copyright 2011 - 2016, Met Office
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

from __future__ import (absolute_import, division, print_function)

import unittest
import warnings

import numpy as np
from numpy.testing import assert_array_equal

import cartopy.crs as ccrs
import cartopy.io.srtm
from cartopy.tests.io.test_downloaders import download_to_temp


@unittest.skip('SRTM login not supported')
class TestRetrieve(unittest.TestCase):
    def _test_srtm_retrieve(self, Source, read_SRTM, max_, min_, pt):
        # test that the download mechanism for SRTM works
        with download_to_temp() as tmp_dir:
            with warnings.catch_warnings(record=True) as w:
                r = Source().srtm_fname(-4, 50)
                self.assertEqual(len(w), 1)
                self.assertTrue(issubclass(w[0].category,
                                           cartopy.io.DownloadWarning))

            self.assertTrue(r.startswith(tmp_dir),
                            'File not downloaded to tmp dir')

            img, _, _ = read_SRTM(r)

            # check that the data is fairly sensible
            msg = ('SRTM data has changed. Arbitrary value testing failed.'
                   ' Got {}.')
            self.assertEqual(img.max(), max_, msg=msg.format(img.max()))
            self.assertEqual(img.min(), min_, msg=msg.format(img.min()))
            self.assertEqual(img[-10, 12], pt, msg=msg.format(img[-10, 12]))

    def test_srtm3_retrieve(self):
        self._test_srtm_retrieve(cartopy.io.srtm.SRTM3Source,
                                 cartopy.io.srtm.read_SRTM3,
                                 602, -34, 78)

    def test_srtm1_retrieve(self):
        self._test_srtm_retrieve(cartopy.io.srtm.SRTM1Source,
                                 cartopy.io.srtm.read_SRTM1,
                                 602, -37, 50)

    def _test_srtm_out_of_range(self, Source, shape):
        # Somewhere over the pacific the elevation should be 0.
        img, _, _ = Source().combined(120, 2, 2, 2)
        assert_array_equal(img, np.zeros(np.array(shape) * 2))

    def test_srtm3_out_of_range(self):
        _test_srtm_out_of_range(self,
                                cartopy.io.srtm.SRTM3Source, (1201, 1201))

    def test_srtm1_out_of_range(self):
        _test_srtm_out_of_range(self,
                                cartopy.io.srtm.SRTM1Source, (3601, 3601))


@unittest.skip('SRTM login not supported')
class TestSRTMSource__single_tile(unittest.TestCase):
    def _out_of_range(self, source):
        msg = 'No srtm tile found for those coordinates.'
        with self.assertRaisesRegexp(ValueError, msg):
            source.single_tile(-25, 50)

    def test_out_of_range3(self):
        self._out_of_range(cartopy.io.srtm.SRTM3Source())

    def test_out_of_range1(self):
        self._out_of_range(cartopy.io.srtm.SRTM1Source())

    def _in_range(self, source, shape):
        img, crs, extent = source.single_tile(-1, 50)
        self.assertIsInstance(img, np.ndarray)
        self.assertEqual(img.shape, shape)
        self.assertEqual(img.dtype, np.dtype('>i2'))
        self.assertEqual(crs, ccrs.PlateCarree())
        self.assertEqual(extent, (-1, 0, 50, 51))

    def test_in_range3(self):
        self._in_range(cartopy.io.srtm.SRTM3Source(), (1201, 1201))

    def test_in_range1(self):
        self._in_range(cartopy.io.srtm.SRTM1Source(), (3601, 3601))

    def _zeros(self, source):
        _, _, extent = source.single_tile(0, 50)
        self.assertEqual(extent, (0, 1, 50, 51))

    def test_zeros3(self):
        self._zeros(cartopy.io.srtm.SRTM3Source())

    def test_zeros1(self):
        self._zeros(cartopy.io.srtm.SRTM1Source())


@unittest.skip('SRTM login not supported')
class TestSRTMSource__combined(unittest.TestCase):
    def _trivial(self, source):
        e_img, e_crs, e_extent = source.single_tile(-3, 50)
        r_img, r_crs, r_extent = source.combined(-3, 50, 1, 1)
        assert_array_equal(e_img, r_img)
        self.assertEqual(e_crs, r_crs)
        self.assertEqual(e_extent, r_extent)

    def test_trivial3(self):
        self._trivial(cartopy.io.srtm.SRTM3Source())

    def test_trivial1(self):
        self._trivial(cartopy.io.srtm.SRTM1Source())

    def _2by2(self, source):
        e_img, _, e_extent = source.combined(-1, 50, 2, 1)
        self.assertEqual(e_extent, (-1, 1, 50, 51))
        imgs = [source.single_tile(-1, 50)[0],
                source.single_tile(0, 50)[0]]
        assert_array_equal(np.hstack(imgs), e_img)

    def test_2by2_3(self):
        self._2by2(cartopy.io.srtm.SRTM3Source())

    def test_2by2_1(self):
        self._2by2(cartopy.io.srtm.SRTM1Source())


@unittest.skip('SRTM login not supported')
class TestSRTM3Source_fetch_raster(unittest.TestCase):
    def _as_combined(self, source):
        e_img, e_crs, e_extent = source.combined(-1, 50, 2, 1)
        imgs = source.fetch_raster(ccrs.PlateCarree(),
                                   (-0.9, 0.1, 50.1, 50.999),
                                   None)
        self.assertEqual(len(imgs), 1)
        r_img, r_extent = imgs[0]
        self.assertEqual(e_extent, r_extent)
        assert_array_equal(e_img[::-1, :], r_img)

    def test_as_combined3(self):
        self._as_combined(cartopy.io.srtm.SRTM3Source())

    def test_as_combined1(self):
        self._as_combined(cartopy.io.srtm.SRTM1Source())


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
