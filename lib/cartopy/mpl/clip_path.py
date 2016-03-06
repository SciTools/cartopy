# (C) British Crown Copyright 2013 - 2016, Met Office
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

import matplotlib.path as mpath
import numpy as np


def clip_path_python(subject, clip, point_inside_clip_path):
    """
    Clip the subject path with the given clip path using the
    Sutherland-Hodgman polygon clipping algorithm.

    Args:

    * subject - The subject path to be clipped. Must be a simple, single
                polygon path with straight line segments only.
    * clip - The clip path to use. Must be a simple, single
             polygon path with straight line segments only.
    * point_inside_clip_path - a point which can be found inside the clip path
                               polygon.

    """
    inside_pt = point_inside_clip_path

    output_verts = subject.vertices

    for i in xrange(clip.vertices.shape[0] - 1):
        clip_edge = clip.vertices[i:i + 2, :]
        input_verts = output_verts
        output_verts = []
        inside = np.cross(clip_edge[1, :] - clip_edge[0, :],
                          inside_pt - clip_edge[0, :])

        try:
            s = input_verts[-1]
        except IndexError:
            break

        for e in input_verts:
            e_clip_cross = np.cross(clip_edge[1, :] - clip_edge[0, :],
                                    e - clip_edge[0, :])
            s_clip_cross = np.cross(clip_edge[1, :] - clip_edge[0, :],
                                    s - clip_edge[0, :])

            if np.sign(e_clip_cross) == np.sign(inside):
                if np.sign(s_clip_cross) != np.sign(inside):
                    p = intersection_point(clip_edge[0, :], clip_edge[1, :],
                                           e, s)
                    output_verts.append(p)
                output_verts.append(e)
            elif np.sign(s_clip_cross) == np.sign(inside):
                p = intersection_point(clip_edge[0, :], clip_edge[1, :],
                                       e, s)
                output_verts.append(p)
            s = e

    if output_verts == []:
        path = mpath.Path([[0, 0]], codes=[mpath.Path.MOVETO])
    else:
        # If the subject polygon was closed, then the return should be too.
        if np.all(subject.vertices[0, :] == subject.vertices[-1, :]):
            output_verts.append(output_verts[0])
        path = mpath.Path(np.array(output_verts))
    return path


def intersection_point(p0, p1, p2, p3):
    """
    Return the intersection point of the two infinite lines that pass through
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
if hasattr(mpath.Path, 'clip_to_bbox'):
    def clip_path(subject, clip_bbox):
        """
        Clip the given path to the given bounding box.

        """
        return subject.clip_to_bbox(clip_bbox)
else:
    def clip_path(subject, clip_bbox):
        """
        Clip the given path to the given bounding box.

        """
        # A shim on clip_path_python to support Bbox path clipping.

        bbox_patch = bbox_to_path(clip_bbox)
        bbox_center = ((clip_bbox.x0 + clip_bbox.x1) / 2,
                       (clip_bbox.y0 + clip_bbox.y1) / 2)
        return clip_path_python(subject, bbox_patch, bbox_center)


def lines_intersect(p0, p1, p2, p3):
    """
    Return whether the two lines defined by p0->p1 and p2->p3 intersect.
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
