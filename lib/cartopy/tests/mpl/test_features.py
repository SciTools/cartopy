# (C) British Crown Copyright 2011 - 2016, Met Office
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

import unittest

import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.ogc_clients import _OWSLIB_AVAILABLE

from cartopy.tests.mpl import ImageTesting


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


@ImageTesting(['natural_earth_custom'])
def test_natural_earth_custom():
    ax = plt.axes(projection=ccrs.PlateCarree())
    feature = cfeature.NaturalEarthFeature('physical', 'coastline', '50m',
                                           edgecolor='black',
                                           facecolor='none')
    ax.add_feature(feature)
    ax.set_xlim((-26, -12))
    ax.set_ylim((58, 72))


@ImageTesting(['gshhs_coastlines'])
def test_gshhs():
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.set_extent([138, 142, 32, 42], ccrs.Geodetic())

    ax.stock_img()
    # Draw coastlines.
    ax.add_feature(cfeature.GSHHSFeature('coarse', edgecolor='red'))
    # Draw higher resolution lakes (and test overriding of kwargs)
    ax.add_feature(cfeature.GSHHSFeature('low', levels=[2],
                                         facecolor='green'), facecolor='blue')


@unittest.skipIf(not _OWSLIB_AVAILABLE, 'OWSLib is unavailable.')
@ImageTesting(['wfs'])
def test_wfs():
    ax = plt.axes(projection=ccrs.OSGB())
    url = 'https://nsidc.org/cgi-bin/atlas_south?service=WFS'
    typename = 'land_excluding_antarctica'
    feature = cfeature.WFSFeature(url, typename,
                                  edgecolor='red')
    ax.add_feature(feature)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
