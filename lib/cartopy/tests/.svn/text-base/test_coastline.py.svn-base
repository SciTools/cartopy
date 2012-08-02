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


import itertools
import os.path
import unittest

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pyppin
import shapely.geometry as sgeom

import cartopy


COASTLINE_PATH = os.path.join(os.path.dirname(cartopy.__file__), 'data', '110m_coastline')


class TestCoastline(unittest.TestCase):
    def test_robust(self):
        # Make sure all the coastlines can be projected without raising any
        # exceptions.
        projection = cartopy.crs.TransverseMercator(central_longitude=-90)
        reader = pyppin.Reader(COASTLINE_PATH)
        all_geometries = list(reader.geometries())
        geometries = []
        geometries += all_geometries
        #geometries += all_geometries[48:52] # Aus & Taz
        #geometries += all_geometries[72:73] # GB
        #for geometry in geometries:
        for i, geometry in enumerate(geometries[93:]):
            for line_string in geometry:
                multi_line_string = projection.project_geometry(line_string)


if __name__ == '__main__':
    unittest.main()
