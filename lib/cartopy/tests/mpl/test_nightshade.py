# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from datetime import datetime

import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
from cartopy.tests.mpl import ImageTesting


@pytest.mark.natural_earth
@ImageTesting(['nightshade_platecarree'])
def test_nightshade_image():
    # Test the actual creation of the image
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    dt = datetime(2018, 11, 10, 0, 0)
    ax.set_global()
    ax.add_feature(Nightshade(dt, alpha=0.75))
