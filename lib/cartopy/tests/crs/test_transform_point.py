# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
import warnings

import cartopy.crs as ccrs


def test_transform_point():
    # see https://github.com/SciTools/cartopy/pull/2194
    p = ccrs.PlateCarree()
    p2 = ccrs.Mercator()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        p2.transform_point(1, 2, p)
