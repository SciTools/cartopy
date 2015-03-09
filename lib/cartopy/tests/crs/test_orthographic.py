# (C) British Crown Copyright 2015, Met Office
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

from __future__ import (absolute_import, division, print_function)

import cartopy.crs as ccrs


def test_ortho_geometry_transform():
    # Traditionally projecting geometries to Orthographic
    # as been the main source of transform issues. Let's test
    # several cases here to ensure that at least for orthographic,
    # changes to the transform functionality don't break known cases.
    import cartopy.io.shapereader as shpreader
    import shapely.geometry as sgeom

    shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='physical',
                                          name='ocean')
    oceans = list(shpreader.Reader(shpfilename).geometries())

    # Cases are ortho keywords vs the % of the visible globe which is ocean.
    cases = ([dict(central_longitude=-33, central_latitude=-80), 81],
             [dict(central_longitude=-89.719, central_latitude=10.51), 75],
             [dict(central_longitude=1.530, central_latitude=2.034), 61.8],
             [dict(central_longitude=-27.394, central_latitude=-5.054), 66.8],
             [dict(central_longitude=0, central_latitude=-75), 80.2],
             )

    problems = []

    for case, expected_area in cases:
        crs = ccrs.Orthographic(**case)
        # We add together the areas of each of the oceans projected (in the
        # same way that the Feature functionality does), and will compute the
        # fraction of the visible area which is ocean.
        total_area = 0

        geoms = []
        for ocean in oceans:
            projected = crs.project_geometry(ocean, ccrs.PlateCarree())
            if not projected.is_empty:
                geoms.append(projected)
            total_area += projected.area

        full_disk_area = sgeom.Polygon(crs.boundary).area
        actual_area = (total_area / full_disk_area) * 100.
        if not 0 < abs(actual_area - expected_area) < 0.5:
            problems.append('For case {}, expected to find an area ~{} but '
                            'actually got {}.'.format(case, expected_area,
                                                      actual_area))
    if problems:
        raise ValueError('Problems with transformations: \n  {}'
                         ''.format('\n  '.join(problems)))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
