# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt
from matplotlib.testing.decorators import cleanup
import pytest

import cartopy.crs as ccrs
from cartopy.tests.mpl import ImageTesting


@pytest.mark.natural_earth
@ImageTesting(['lambert_conformal_south'])
def test_lambert_south():
    # Reference image: https://www.icsm.gov.au/mapping/map_projections.html
    crs = ccrs.LambertConformal(central_longitude=140, cutoff=65,
                                standard_parallels=(-30, -60))
    ax = plt.axes(projection=crs)
    ax.coastlines()
    ax.gridlines()


@pytest.mark.natural_earth
@ImageTesting(['mercator_squashed'])
def test_mercator_squashed():
    globe = ccrs.Globe(semimajor_axis=10000, semiminor_axis=9000,
                       ellipse=None)
    crs = ccrs.Mercator(globe=globe, min_latitude=-40, max_latitude=40)
    ax = plt.axes(projection=crs)
    ax.coastlines()
    ax.gridlines()


@pytest.mark.natural_earth
@cleanup
def test_repr_html():
    pc = ccrs.PlateCarree()
    html = pc._repr_html_()

    assert html is not None
    assert '<svg ' in html
    assert '<pre>&lt;cartopy.crs.PlateCarree object at ' in html
