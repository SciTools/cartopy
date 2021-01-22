# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import numpy as np
import numpy.ma as ma
from numpy.testing import assert_array_equal
import pytest

from cartopy.util import add_cyclic_point


class Test_add_cyclic_point:

    @classmethod
    def setup_class(cls):
        cls.lons = np.arange(0, 360, 60)
        cls.data2d = np.ones([3, 6]) * np.arange(6)
        cls.data4d = np.ones([4, 6, 2, 3]) * \
            np.arange(6)[..., np.newaxis, np.newaxis]

    def test_data_only(self):
        c_data = add_cyclic_point(self.data2d)
        r_data = np.concatenate((self.data2d, self.data2d[:, :1]), axis=1)
        assert_array_equal(c_data, r_data)

    def test_data_and_coord(self):
        c_data, c_lons = add_cyclic_point(self.data2d, coord=self.lons)
        r_data = np.concatenate((self.data2d, self.data2d[:, :1]), axis=1)
        r_lons = np.concatenate((self.lons, np.array([360.])))
        assert_array_equal(c_data, r_data)
        assert_array_equal(c_lons, r_lons)

    def test_data_only_with_axis(self):
        c_data = add_cyclic_point(self.data4d, axis=1)
        r_data = np.concatenate((self.data4d, self.data4d[:, :1]), axis=1)
        assert_array_equal(c_data, r_data)

    def test_data_and_coord_with_axis(self):
        c_data, c_lons = add_cyclic_point(self.data4d, coord=self.lons, axis=1)
        r_data = np.concatenate((self.data4d, self.data4d[:, :1]), axis=1)
        r_lons = np.concatenate((self.lons, np.array([360.])))
        assert_array_equal(c_data, r_data)
        assert_array_equal(c_lons, r_lons)

    def test_masked_data(self):
        new_data = ma.masked_less(self.data2d, 3)
        c_data = add_cyclic_point(new_data)
        r_data = ma.concatenate((self.data2d, self.data2d[:, :1]), axis=1)
        assert_array_equal(c_data, r_data)

    def test_invalid_coord_dimensionality(self):
        lons2d = np.repeat(self.lons[np.newaxis], 3, axis=0)
        with pytest.raises(ValueError):
            c_data, c_lons = add_cyclic_point(self.data2d, coord=lons2d)

    def test_invalid_coord_size(self):
        with pytest.raises(ValueError):
            c_data, c_lons = add_cyclic_point(self.data2d,
                                              coord=self.lons[:-1])

    def test_invalid_axis(self):
        with pytest.raises(ValueError):
            add_cyclic_point(self.data2d, axis=-3)
