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

import contextlib
import functools
import tempfile
import shutil
import types

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


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


def show(projection, geometry):
    if geometry.type == 'MultiPolygon' and 1:
        multi_polygon = geometry
        for polygon in multi_polygon:
            import cartopy.mpl.patch as patch
            paths = patch.geos_to_path(polygon)
            for pth in paths:
                patch = mpatches.PathPatch(pth, edgecolor='none',
                                           lw=0, alpha=0.2)
                plt.gca().add_patch(patch)
            line_string = polygon.exterior
            plt.plot(*list(zip(*line_string.coords)),
                     marker='+', linestyle='-')
    elif geometry.type == 'MultiPolygon':
        multi_polygon = geometry
        for polygon in multi_polygon:
            line_string = polygon.exterior
            plt.plot(*list(zip(*line_string.coords)),
                     marker='+', linestyle='-')

    elif geometry.type == 'MultiLineString':
        multi_line_string = geometry
        for line_string in multi_line_string:
            plt.plot(*list(zip(*line_string.coords)),
                     marker='+', linestyle='-')

    elif geometry.type == 'LinearRing':
        plt.plot(*list(zip(*geometry.coords)), marker='+', linestyle='-')

    if 1:
        # Whole map domain
        plt.autoscale()
    elif 0:
        # The left-hand triangle
        plt.xlim(-1.65e7, -1.2e7)
        plt.ylim(0.3e7, 0.65e7)
    elif 0:
        # The tip of the left-hand triangle
        plt.xlim(-1.65e7, -1.55e7)
        plt.ylim(0.3e7, 0.4e7)
    elif 1:
        # The very tip of the left-hand triangle
        plt.xlim(-1.632e7, -1.622e7)
        plt.ylim(0.327e7, 0.337e7)
    elif 1:
        # The tip of the right-hand triangle
        plt.xlim(1.55e7, 1.65e7)
        plt.ylim(0.3e7, 0.4e7)

    plt.plot(*list(zip(*projection.boundary.coords)), marker='o',
             scalex=False, scaley=False, zorder=-1)

    plt.show()
