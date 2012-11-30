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
"""
This module defines :class:`Feature` instances, for use with
ax.add_feature().

"""
import numpy as np

import cartopy.crs


_COLOURS = {'land': np.array((240, 240, 220)) / 256.,
            'water': np.array((152, 183, 226)) / 256.}


_NATURAL_EARTH_GEOM_CACHE = {}
"""
Caches a mapping between (name, category, scale) and a tuple of the
resulting geometries.

Provides a significant performance benefit (when combined with object id
caching in GeoAxes.add_geometries) when producing multiple maps of the
same projection.

"""
from abc import ABCMeta, abstractmethod


class Feature(object):
    """
    The abstract base class for features.

    """

    __metaclass__ = ABCMeta

    def __init__(self, crs, kwargs):
        self._crs = crs
        self._kwargs = dict(kwargs)

    @property
    def crs(self):
        """The cartopy CRS for the geometries in this feature."""
        return self._crs

    @property
    def kwargs(self):
        """
        A dictionary of keyword arguments to be used when creating
        the matplotlib artists for this feature.

        """
        return dict(self._kwargs)

    @abstractmethod
    def geometries(self):
        """
        Must be overriden to return the shapely geometries for this
        dataset.

        """
        raise NotImplementedError()


class NaturalEarthFeature(Feature):
    """
    A simple interface to Natural Earth shapefiles.

    See http://www.naturalearthdata.com/

    """
    def __init__(self, category, name, scale, kwargs):
        """
        Args:

        * category:
            The category of the dataset, i.e. either 'cultural' or 'physical'.
        * name:
            The name of the dataset, e.g. 'admin-0-boundary-lines'.
        * scale:
            The dataset scale, i.e. one of '10m', '50m', or '110m'.
            Corresponding to 1:10,000,000, 1:50,000,000, and 1:110,000,000
            respectively.
        * kwargs:
            A dictionary of keyword arguments to be used when creating
            the matplotlib artists.

        """
        super(NaturalEarthFeature, self).__init__(cartopy.crs.PlateCarree(),
                                                  kwargs)
        self.category = category
        self.name = name
        self.scale = scale

    def geometries(self):
        """
        Returns the shapely geometries defined by this Natural
        Earth dataset.

        """
        import cartopy.io.shapereader as shapereader
        key = (self.name, self.category, self.scale)
        if key not in _NATURAL_EARTH_GEOM_CACHE:
            path = shapereader.natural_earth(resolution=self.scale,
                                             category=self.category,
                                             name=self.name)
            geometries = tuple(shapereader.Reader(path).geometries())
            _NATURAL_EARTH_GEOM_CACHE[key] = geometries
        else:
            geometries = _NATURAL_EARTH_GEOM_CACHE[key]
        return geometries


BORDERS = NaturalEarthFeature('cultural',
                              'admin_0_boundary_lines_land',
                              '110m',
                              {'edgecolor': 'black', 'facecolor': 'none'})
"""Small scale (1:110m) country boundaries."""


COASTLINE = NaturalEarthFeature('physical', 'coastline', '110m',
                                {'edgecolor': 'black', 'facecolor': 'none'})
"""Small scale (1:110m) coastline, including major islands."""


LAKES = NaturalEarthFeature('physical', 'lakes', '110m',
                            {'edgecolor': 'face',
                             'facecolor': _COLOURS['water']})
"""Small scale (1:110m) natural and artificial lakes."""


LAND = NaturalEarthFeature('physical', 'land', '110m',
                           {'edgecolor': 'face',
                            'facecolor': _COLOURS['land']})
"""Small scale (1:110m) land polygons, including major islands."""


OCEAN = NaturalEarthFeature('physical', 'ocean', '110m',
                            {'edgecolor': 'face',
                             'facecolor': _COLOURS['water']})
"""Small scale (1:110m) ocean polygons."""


RIVERS = NaturalEarthFeature('physical', 'rivers_lake_centerlines', '110m',
                             {'edgecolor': _COLOURS['water'],
                              'facecolor': 'none'})
"""Small scale (1:110m) single-line drainages, including lake centerlines."""
