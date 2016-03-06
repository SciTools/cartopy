# (C) British Crown Copyright 2014 - 2016, Met Office
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
"""
Provides support for converting EPSG codes to Projection instances.

"""

from __future__ import (absolute_import, division, print_function)

import cartopy.crs as ccrs
import numpy as np
import pyepsg
import shapely.geometry as sgeom


_GLOBE_PARAMS = {'datum': 'datum',
                 'ellps': 'ellipse',
                 'a': 'semimajor_axis',
                 'b': 'semiminor_axis',
                 'f': 'flattening',
                 'rf': 'inverse_flattening',
                 'towgs84': 'towgs84',
                 'nadgrids': 'nadgrids'}


class _EPSGProjection(ccrs.Projection):
    def __init__(self, code):
        projection = pyepsg.get(code)
        if not isinstance(projection, pyepsg.ProjectedCRS):
            raise ValueError('EPSG code does not define a projection')

        self.epsg_code = code

        proj4_str = projection.as_proj4().strip()
        terms = [term.strip('+').split('=') for term in proj4_str.split(' ')]
        globe_terms = filter(lambda term: term[0] in _GLOBE_PARAMS, terms)
        globe = ccrs.Globe(**{_GLOBE_PARAMS[name]: value for name, value in
                              globe_terms})
        other_terms = []
        for term in terms:
            if term[0] not in _GLOBE_PARAMS:
                if len(term) == 1:
                    other_terms.append([term[0], None])
                else:
                    other_terms.append(term)
        super(_EPSGProjection, self).__init__(other_terms, globe)

        # Convert lat/lon bounds to projected bounds.
        # GML defines gmd:EX_GeographicBoundingBox as:
        #   Geographic area of the entire dataset referenced to WGS 84
        # NB. We can't use a polygon transform at this stage because
        # that relies on the existence of the map boundary... the very
        # thing we're trying to work out! ;-)
        x0, x1, y0, y1 = projection.domain_of_validity()
        geodetic = ccrs.Geodetic()
        lons = np.array([x0, x0, x1, x1])
        lats = np.array([y0, y1, y1, y0])
        points = self.transform_points(geodetic, lons, lats)
        x = points[:, 0]
        y = points[:, 1]
        self.bounds = (x.min(), x.max(), y.min(), y.max())

    def __repr__(self):
        return '_EPSGProjection({})'.format(self.epsg_code)

    @property
    def boundary(self):
        x0, x1, y0, y1 = self.bounds
        return sgeom.LineString([(x0, y0), (x0, y1), (x1, y1), (x1, y0),
                                 (x0, y0)])

    @property
    def x_limits(self):
        x0, x1, y0, y1 = self.bounds
        return (x0, x1)

    @property
    def y_limits(self):
        x0, x1, y0, y1 = self.bounds
        return (y0, y1)

    @property
    def threshold(self):
        x0, x1, y0, y1 = self.bounds
        return min(x1 - x0, y1 - y0) / 100.
