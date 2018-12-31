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
Pytest configuration for Cartopy tests.

"""

from __future__ import (absolute_import, division, print_function)

import os

import pytest


@pytest.fixture(scope='module')
def vcr_cassette_dir(request):
    module_path = request.module.__name__[len('cartopy.tests.'):].split('.')
    return os.path.join(os.path.dirname(__file__),
                        'cassettes',
                        *module_path)
