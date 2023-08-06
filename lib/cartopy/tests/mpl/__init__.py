# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import packaging.version


MPL_VERSION = packaging.version.parse(mpl.__version__)


def show(projection, geometry):
    orig_backend = mpl.get_backend()
    plt.switch_backend('tkagg')

    if geometry.geom_type == 'MultiPolygon' and 1:
        multi_polygon = geometry
        for polygon in multi_polygon:
            import cartopy.mpl.patch as patch
            paths = patch.geos_to_path(polygon)
            for pth in paths:
                patch = mpatches.PathPatch(pth, edgecolor='none',
                                           lw=0, alpha=0.2)
                plt.gca().add_patch(patch)
            line_string = polygon.exterior
            plt.plot(*zip(*line_string.coords),
                     marker='+', linestyle='-')
    elif geometry.geom_type == 'MultiPolygon':
        multi_polygon = geometry
        for polygon in multi_polygon:
            line_string = polygon.exterior
            plt.plot(*zip(*line_string.coords),
                     marker='+', linestyle='-')

    elif geometry.geom_type == 'MultiLineString':
        multi_line_string = geometry
        for line_string in multi_line_string:
            plt.plot(*zip(*line_string.coords),
                     marker='+', linestyle='-')

    elif geometry.geom_type == 'LinearRing':
        plt.plot(*zip(*geometry.coords), marker='+', linestyle='-')

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

    plt.plot(*zip(*projection.boundary.coords), marker='o',
             scalex=False, scaley=False, zorder=-1)

    plt.show()
    plt.switch_backend(orig_backend)
