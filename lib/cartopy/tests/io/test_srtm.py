# (C) British Crown Copyright 2011 - 2018, Met Office
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

import warnings

import numpy as np
from numpy.testing import assert_array_equal
import pytest

import cartopy.crs as ccrs
import cartopy.io.srtm

from .test_downloaders import download_to_temp  # noqa: F401 (used as fixture)


pytestmark = [pytest.mark.network,
              pytest.mark.filterwarnings('ignore:SRTM requires an account'),
              pytest.mark.usefixtures('srtm_login_or_skip')]


@pytest.fixture
def srtm_login_or_skip(monkeypatch):
    import os
    try:
        srtm_username = os.environ['SRTM_USERNAME']
    except KeyError:
        pytest.skip('SRTM_USERNAME environment variable is unset.')
    try:
        srtm_password = os.environ['SRTM_PASSWORD']
    except KeyError:
        pytest.skip('SRTM_PASSWORD environment variable is unset.')

    from six.moves.urllib.request import (HTTPBasicAuthHandler,
                                          HTTPCookieProcessor,
                                          HTTPPasswordMgrWithDefaultRealm,
                                          build_opener)
    from six.moves.http_cookiejar import CookieJar

    password_manager = HTTPPasswordMgrWithDefaultRealm()
    password_manager.add_password(
        None,
        "https://urs.earthdata.nasa.gov",
        srtm_username,
        srtm_password)
    cookie_jar = CookieJar()
    opener = build_opener(HTTPBasicAuthHandler(password_manager),
                          HTTPCookieProcessor(cookie_jar))

    monkeypatch.setattr(cartopy.io, 'urlopen', opener.open)


class TestRetrieve(object):
    @pytest.mark.parametrize('Source, read_SRTM, max_, min_, pt', [
        (cartopy.io.srtm.SRTM3Source, cartopy.io.srtm.read_SRTM3,
         602, -34, 78),
        (cartopy.io.srtm.SRTM1Source, cartopy.io.srtm.read_SRTM1,
         602, -37, 50),
    ], ids=[
        'srtm3',
        'srtm1',
    ])
    def test_srtm_retrieve(self, Source, read_SRTM, max_, min_, pt,
                           download_to_temp):
        # test that the download mechanism for SRTM works
        with warnings.catch_warnings(record=True) as w:
            r = Source().srtm_fname(-4, 50)
            assert len(w) == 1
            assert issubclass(w[0].category, cartopy.io.DownloadWarning)

        assert r.startswith(str(download_to_temp)), \
            'File not downloaded to tmp dir'

        img, _, _ = read_SRTM(r)

        # check that the data is fairly sensible
        assert img.max() == max_
        assert img.min() == min_
        assert img[-10, 12] == pt

    @pytest.mark.parametrize('Source, shape', [
        (cartopy.io.srtm.SRTM3Source, (1201, 1201)),
        (cartopy.io.srtm.SRTM1Source, (3601, 3601)),
    ], ids=[
        'srtm3',
        'srtm1',
    ])
    def test_srtm_out_of_range(self, Source, shape):
        # Somewhere over the pacific the elevation should be 0.
        img, _, _ = Source().combined(120, 2, 2, 2)
        assert_array_equal(img, np.zeros(np.array(shape) * 2))


@pytest.mark.parametrize('Source', [
    cartopy.io.srtm.SRTM3Source,
    cartopy.io.srtm.SRTM1Source,
], ids=[
    'srtm3',
    'srtm1',
])
class TestSRTMSource__single_tile(object):
    def test_out_of_range(self, Source):
        source = Source()
        msg = 'No srtm tile found for those coordinates.'
        with pytest.raises(ValueError, message=msg):
            source.single_tile(-25, 50)

    def test_in_range(self, Source):
        if Source == cartopy.io.srtm.SRTM3Source:
            shape = (1201, 1201)
        elif Source == cartopy.io.srtm.SRTM1Source:
            shape = (3601, 3601)
        else:
            raise ValueError('Source is of unexpected type.')
        source = Source()
        img, crs, extent = source.single_tile(-1, 50)
        assert isinstance(img, np.ndarray)
        assert img.shape == shape
        assert img.dtype == np.dtype('>i2')
        assert crs == ccrs.PlateCarree()
        assert extent == (-1, 0, 50, 51)

    def test_zeros(self, Source):
        source = Source()
        _, _, extent = source.single_tile(0, 50)
        assert extent == (0, 1, 50, 51)


@pytest.mark.parametrize('Source', [
    cartopy.io.srtm.SRTM3Source,
    cartopy.io.srtm.SRTM1Source,
], ids=[
    'srtm3',
    'srtm1',
])
class TestSRTMSource__combined(object):
    def test_trivial(self, Source):
        source = Source()

        e_img, e_crs, e_extent = source.single_tile(-3, 50)
        r_img, r_crs, r_extent = source.combined(-3, 50, 1, 1)
        assert_array_equal(e_img, r_img)
        assert e_crs == r_crs
        assert e_extent == r_extent

    def test_2by2(self, Source):
        source = Source()

        e_img, _, e_extent = source.combined(-1, 50, 2, 1)
        assert e_extent == (-1, 1, 50, 51)
        imgs = [source.single_tile(-1, 50)[0],
                source.single_tile(0, 50)[0]]
        assert_array_equal(np.hstack(imgs), e_img)


@pytest.mark.parametrize('Source', [
    cartopy.io.srtm.SRTM3Source,
    cartopy.io.srtm.SRTM1Source,
], ids=[
    'srtm3',
    'srtm1',
])
def test_fetch_raster_ascombined(Source):
    source = Source()

    e_img, e_crs, e_extent = source.combined(-1, 50, 2, 1)
    imgs = source.fetch_raster(ccrs.PlateCarree(),
                               (-0.9, 0.1, 50.1, 50.999),
                               None)
    assert len(imgs) == 1
    r_img, r_extent = imgs[0]
    assert e_extent == r_extent
    assert_array_equal(e_img[::-1, :], r_img)
