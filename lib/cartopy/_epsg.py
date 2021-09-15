# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
"""
Provide support for converting EPSG codes to Projection instances.

"""
import cartopy.crs as ccrs
from pyproj.crs import CRS as _CRS


class _EPSGProjection(ccrs.Projection):
    def __init__(self, code):
        crs = _CRS.from_epsg(code)
        if not crs.is_projected:
            raise ValueError('EPSG code does not define a projection')
        if not crs.area_of_use:
            raise ValueError("Area of use not defined.")

        self.epsg_code = code
        super().__init__(crs.to_wkt())

    def __repr__(self):
        return f'_EPSGProjection({self.epsg_code})'
