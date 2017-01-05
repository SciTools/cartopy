# (C) British Crown Copyright 2017, Met Office
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

import matplotlib.pyplot as plt
from matplotlib.testing.decorators import cleanup
import numpy as np
from numpy.testing import assert_array_almost_equal

import cartopy.crs as ccrs



# TODO: Decide how I want to test this:
# Graphics test? (Probably not).
# Test for nans or infs in transformed arrays? (Maybe too contrived)
# Some kind of test of the geoaxes (like in test_contour.py) - Check mpl for options
def test_transform_to_orthographic():




def test_transform_to_transverse_mercator():




def test_transform_to_gnomonic():



def test_transform_to_geostationary():




if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
