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
from shapely.geometry import Polygon


@pytest.mark.filterwarnings("ignore:Downloading")
@pytest.mark.natural_earth
@ImageTesting(['natural_earth'])
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


@pytest.mark.natural_earth
@ImageTesting(['UTM_all_zones'])
def test_utm_all():
    # this test fails because of padding (?) differences
    # probably solved with an rcParam
    zones = range(1, 61)
    fig = plt.figure(figsize=(18, 6))
    for zone in zones:
        ax = fig.add_subplot(1, len(zones), zone,
                             projection=ccrs.UTM(zone=zone,
                             southern_hemisphere=True))
        ax.add_feature(cfeature.LAND, facecolor="tab:blue")


@pytest.mark.natural_earth
@ImageTesting(['nonplatecarree_with_projected_hole'])
def test_nonpc_hole():
    lamaz_box = Polygon(((8.5e5, 1.1e6), (-8.5e5, 1.1e6),
                         (-5.6e5, -1.2e6), (5.6e5, -1.2e6),
                         (8.5e5, 1.1e6)))

    fig = plt.figure(figsize=(20, 7))
    proj1 = ccrs.LambertAzimuthalEqualArea(central_longitude=62.,
                                           central_latitude=-50.)
    ax1 = fig.add_subplot(1, 2, 1, projection=proj1)
    ax1.add_geometries([lamaz_box], crs=proj1, color="tab:blue")

    proj2 = ccrs.LambertAzimuthalEqualArea(central_longitude=-118.,
                                           central_latitude=50.)
    ax2 = fig.add_subplot(1, 2, 2, projection=proj2)
    ax2.add_geometries([lamaz_box], crs=proj1, color="tab:blue")
