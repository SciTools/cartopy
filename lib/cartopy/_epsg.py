# (C) British Crown Copyright 2014, Met Office
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
Provides support for converting EPSG codes to Projection instances.

"""
import cartopy.crs as ccrs
import pyepsg


class _EPSGProjection(ccrs._Proj4Projection):
    def __init__(self, code):
        projection = pyepsg.get(code)
        if not isinstance(projection, pyepsg.ProjectedCRS):
            raise ValueError('EPSG code does not define a projection')

        self.epsg_code = code

        proj4_str = projection.as_proj4().strip()
        x0, x1, y0, y1 = projection.domain_of_validity()

        super(_EPSGProjection, self).__init__(proj4_str, x0, x1, y0, y1)


    def __repr__(self):
        return '_EPSGProjection({})'.format(self.epsg_code)

