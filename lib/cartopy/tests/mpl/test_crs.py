# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs
from cartopy.tests.mpl import MPL_VERSION


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename="igh_land.png",
                               tolerance=(3.6
                                          if MPL_VERSION.release[:2] < (3, 5)
                                          else 0.5))
def test_igh_land():
    crs = ccrs.InterruptedGoodeHomolosine(emphasis="land")
    ax = plt.axes(projection=crs)
    ax.coastlines()
    ax.gridlines()
    return ax.figure


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename="igh_ocean.png",
                               tolerance=(4.5
                                          if MPL_VERSION.release[:2] < (3, 5)
                                          else 0.5))
def test_igh_ocean():
    crs = ccrs.InterruptedGoodeHomolosine(
        central_longitude=-160, emphasis="ocean"
    )
    ax = plt.axes(projection=crs)
    ax.coastlines()
    ax.gridlines()
    return ax.figure


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='lambert_conformal_south.png')
def test_lambert_south():
    # Reference image: https://www.icsm.gov.au/mapping/map_projections.html
    crs = ccrs.LambertConformal(central_longitude=140, cutoff=65,
                                standard_parallels=(-30, -60))
    ax = plt.axes(projection=crs)
    ax.coastlines()
    ax.gridlines()
    return ax.figure


@pytest.mark.natural_earth
def test_repr_html():
    pc = ccrs.PlateCarree()
    html = pc._repr_html_()

    assert html is not None
    assert '<svg ' in html
    assert '<pre>&lt;cartopy.crs.PlateCarree object at ' in html
