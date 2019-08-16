# (C) British Crown Copyright 2018, Met Office
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
Helpers for Cartopy CRS subclass tests.

"""

from __future__ import (absolute_import, division, print_function)


def check_proj_params(name, crs, other_args):
    expected = other_args | {'proj=' + name, 'no_defs'}
    proj_params = set(crs.proj4_init.lstrip('+').split(' +'))
    assert expected == proj_params
