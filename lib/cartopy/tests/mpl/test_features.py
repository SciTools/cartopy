# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.pyplot as plt
import pytest

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.ogc_clients import _OWSLIB_AVAILABLE

from cartopy.tests.mpl import ImageTesting


@pytest.mark.filterwarnings("ignore:Downloading")
@pytest.mark.natural_earth
@ImageTesting(['natural_earth'], tolerance=0.97)
def test_natural_earth():
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, alpha=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.set_xlim((-20, 60))
    ax.set_ylim((-40, 40))


@pytest.mark.filterwarnings("ignore:Downloading")
@pytest.mark.natural_earth
@ImageTesting(['natural_earth_custom'])
def test_natural_earth_custom():
    ax = plt.axes(projection=ccrs.PlateCarree())
    feature = cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                           edgecolor='black',
                                           facecolor='none')
    ax.add_feature(feature)
    ax.set_xlim((-26, -12))
    ax.set_ylim((58, 72))


@ImageTesting(['gshhs_coastlines'], tolerance=0.95)
def test_gshhs():
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.set_extent([138, 142, 32, 42], ccrs.Geodetic())

    ax.stock_img()
    # Draw coastlines.
    ax.add_feature(cfeature.GSHHSFeature('coarse', edgecolor='red'))
    # Draw higher resolution lakes (and test overriding of kwargs)
    ax.add_feature(cfeature.GSHHSFeature('low', levels=[2],
                                         facecolor='green'), facecolor='blue')


@pytest.mark.network
@pytest.mark.skipif(not _OWSLIB_AVAILABLE, reason='OWSLib is unavailable.')
@ImageTesting(['wfs'])
def test_wfs():
    ax = plt.axes(projection=ccrs.OSGB(approx=True))
    url = 'https://nsidc.org/cgi-bin/atlas_south?service=WFS'
    typename = 'land_excluding_antarctica'
    feature = cfeature.WFSFeature(url, typename,
                                  edgecolor='red')
    ax.add_feature(feature)
