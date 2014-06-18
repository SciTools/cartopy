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
from __future__ import absolute_import

import warnings

import numpy as np
from numpy.testing import assert_array_equal

import cartopy.io.srtm
from cartopy.tests.io.test_downloaders import download_to_temp


def test_srtm3_retrieve():
    # test that the download mechanism for srtm3 works
    with download_to_temp() as tmp_dir:
        with warnings.catch_warnings(record=True) as w:
            r = cartopy.io.srtm.SRTM3_retrieve(-4, 50)
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
    img, _, _ = cartopy.io.srtm.srtm_composite(120, 2, 2, 2)
    assert_array_equal(img, np.zeros(np.array((1201, 1201)) * 2))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-sv', '--with-doctest'], exit=False)
