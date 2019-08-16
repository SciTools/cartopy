# (C) British Crown Copyright 2011 - 2019, Met Office
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

from matplotlib.contour import QuadContourSet
import matplotlib.path as mpath
import numpy as np


class GeoContourSet(QuadContourSet):
    """
    A contourset designed to handle things like contour labels.

    """
    # nb. No __init__ method here - most of the time a GeoContourSet will
    # come from GeoAxes.contour[f]. These methods morph a ContourSet by
    # fiddling with instance.__class__.

    def clabel(self, *args, **kwargs):
        # nb: contour labelling does not work very well for filled
        # contours - it is recommended to only label line contours.
        # This is especially true when inline=True.

        # This wrapper exist because mpl does not properly transform
        # paths. Instead it simply assumes one path represents one polygon
        # (not necessarily the case), and it assumes that
        # transform(path.verts) is equivalent to transform_path(path).
        # Unfortunately there is no way to easily correct this error,
        # so we are forced to pre-transform the ContourSet's paths from
        # the source coordinate system to the axes' projection.
        # The existing mpl code then has a much simpler job of handling
        # pre-projected paths (which can now effectively be transformed
        # naively).

        for col in self.collections:
            # Snaffle the collection's path list. We will change the
            # list in-place (as the contour label code does in mpl).
            paths = col.get_paths()

            data_t = self.ax.transData
            # Define the transform that will take us from collection
            # coordinates through to axes projection coordinates.
            col_to_data = col.get_transform() - data_t

            # Now that we have the transform, project all of this
            # collection's paths.
            new_paths = [col_to_data.transform_path(path) for path in paths]
            new_paths = [path for path in new_paths if path.vertices.size >= 1]

            # The collection will now be referenced in axes projection
            # coordinates.
            col.set_transform(self.ax.transData)

            # Clear the now incorrectly referenced paths.
            del paths[:]

            for path in new_paths:
                if path.vertices.size == 0:
                    # Don't persist empty paths. Let's get rid of them.
                    continue

                # Split the path if it has multiple MOVETO statements.
                codes = np.array(
                    path.codes if path.codes is not None else [0])
                moveto = codes == mpath.Path.MOVETO
                if moveto.sum() <= 1:
                    # This is only one path, so add it to the collection.
                    paths.append(path)
                else:
                    # The first MOVETO doesn't need cutting-out.
                    moveto[0] = False
                    split_locs = np.flatnonzero(moveto)

                    split_verts = np.split(path.vertices, split_locs)
                    split_codes = np.split(path.codes, split_locs)

                    for verts, codes in zip(split_verts, split_codes):
                        # Add this path to the collection's list of paths.
                        paths.append(mpath.Path(verts, codes))

        # Now that we have prepared the collection paths, call on
        # through to the underlying implementation.
        super(GeoContourSet, self).clabel(*args, **kwargs)
