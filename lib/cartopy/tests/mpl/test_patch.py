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

import six
import unittest

from matplotlib.path import Path

import cartopy.mpl.patch as cpatch


class Test_path_to_geos(unittest.TestCase):
    def test_empty_polyong(self):
        p = Path([[0, 0], [0, 0], [0, 0], [0, 0],
                  [0, 0], [0, 0], [0, 0], [0, 0]],
                 codes=[1, 2, 2, 79,
                        1, 2, 2, 79])
        geom = cpatch.path_to_geos(p)
        self.assertEqual(len(geom), 0)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
