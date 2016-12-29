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

import contextlib
import functools
import re
import tempfile
import shutil
import types


@contextlib.contextmanager
def temp_dir(suffix=None):
    if suffix is None:
        suffix = ''
    dname = tempfile.mkdtemp(suffix=suffix)
    try:
        yield dname
    finally:
        shutil.rmtree(dname)


def not_a_nose_fixture(function):
    """
    Provides a decorator to mark a function as not a nose fixture.

    """
    @functools.wraps(function)
    def setup(app):
        if isinstance(app, types.ModuleType):
            return
        return function(app)
    return setup
