# (C) British Crown Copyright 2013 - 2018, Met Office
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

import warnings

import matplotlib.path as mpath
import numpy as np


def intersection_point(p0, p1, p2, p3):
    """
    Returns
    -------
    x, y
        The intersection point of the two infinite lines that pass through
        point p0->p1 and p2->p3 respectively.

    """
    x_1, y_1 = p0
    x_2, y_2 = p1
    x_3, y_3 = p2
    x_4, y_4 = p3

    div = (x_1 - x_2) * (y_3 - y_4) - (y_1 - y_2) * (x_3 - x_4)

    if div == 0:
        raise ValueError('Lines are parallel and cannot '
                         'intersect at any one point.')

    x = ((x_1 * y_2 - y_1 * x_2) * (x_3 - x_4) - (x_1 - x_2) * (x_3 *
         y_4 - y_3 * x_4)) / div
    y = ((x_1 * y_2 - y_1 * x_2) * (y_3 - y_4) - (y_1 - y_2) * (x_3 *
         y_4 - y_3 * x_4)) / div

    return x, y


# Provide a clip_path function which clips the given path to the given Bbox.
# There is inbuilt mpl functionality with v1.2.1 and beyond, but we provide
# a shim here for older mpl versions.
def clip_path(subject, clip_bbox):
    """
    Clip the given path to the given bounding box.

    """
    warnings.warn("This method has been deprecated. "
                  "You can replace ``clip_path(subject, clip_bbox)`` by "
                  "``subject.clip_to_bbox(clip_bbox)``. "
                  "See the \"What's new\" section for v0.17.")
    return subject.clip_to_bbox(clip_bbox)


def lines_intersect(p0, p1, p2, p3):
    """
    Returns
    -------
    bool
        Boolean indicating whether the two lines defined by p0->p1 and p2->p3
        intersect.
    """
    x_1, y_1 = p0
    x_2, y_2 = p1
    x_3, y_3 = p2
    x_4, y_4 = p3

    return (x_1 - x_2) * (y_3 - y_4) - (y_1 - y_2) * (x_3 - x_4) != 0
    cp1 = np.cross(p1 - p0, p2 - p0)
    cp2 = np.cross(p1 - p0, p3 - p0)
    return np.sign(cp1) == np.sign(cp2) and cp1 != 0


def bbox_to_path(bbox):
    """
    Turn the given :class:`matplotlib.transforms.Bbox` instance into
    a :class:`matplotlib.path.Path` instance.

    """
    verts = np.array([[bbox.x0, bbox.y0], [bbox.x1, bbox.y0],
                      [bbox.x1, bbox.y1], [bbox.x0, bbox.y1],
                      [bbox.x0, bbox.y0]])
    return mpath.Path(verts)
