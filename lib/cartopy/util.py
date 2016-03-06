# (C) British Crown Copyright 2014 - 2016, Met Office
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
This module contains utilities that are useful in conjunction with
cartopy.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
import numpy.ma as ma


def add_cyclic_point(data, coord=None, axis=-1):
    """
    Add a cyclic point to an array and optionally a corresponding
    coordinate.

    Args:

    * data:
        An n-dimensional array of data to add a cyclic point to.

    Kwargs:

    * coord:
        A 1-dimensional array which specifies the coordinate values for
        the dimension the cyclic point is to be added to. The coordinate
        values must be regularly spaced.

    * axis:
        Specifies the axis of the data array to add the cyclic point to.
        Defaults to the right-most axis.

    Returns:

    * cyclic_data:
        The data array with a cyclic point added.

    * cyclic_coord:
        The coordinate with a cyclic point, only returned if the coord
        keyword was supplied.

    Examples:

    Adding a cyclic point to a data array, where the cyclic dimension is
    the right-most dimension

    >>> import numpy as np
    >>> data = np.ones([5, 6]) * np.arange(6)
    >>> cyclic_data = add_cyclic_point(data)
    >>> print(cyclic_data)
    [[ 0.  1.  2.  3.  4.  5.  0.]
     [ 0.  1.  2.  3.  4.  5.  0.]
     [ 0.  1.  2.  3.  4.  5.  0.]
     [ 0.  1.  2.  3.  4.  5.  0.]
     [ 0.  1.  2.  3.  4.  5.  0.]]

    Adding a cyclic point to a data array and an associated coordinate

    >>> lons = np.arange(0, 360, 60)
    >>> cyclic_data, cyclic_lons = add_cyclic_point(data, coord=lons)
    >>> print(cyclic_data)
    [[ 0.  1.  2.  3.  4.  5.  0.]
     [ 0.  1.  2.  3.  4.  5.  0.]
     [ 0.  1.  2.  3.  4.  5.  0.]
     [ 0.  1.  2.  3.  4.  5.  0.]
     [ 0.  1.  2.  3.  4.  5.  0.]]
    >>> print(cyclic_lons)
    [  0  60 120 180 240 300 360]

    """
    if coord is not None:
        if coord.ndim != 1:
            raise ValueError('The coordinate must be 1-dimensional.')
        if len(coord) != data.shape[axis]:
            raise ValueError('The length of the coordinate does not match '
                             'the size of the corresponding dimension of '
                             'the data array: len(coord) = {}, '
                             'data.shape[{}] = {}.'.format(
                                 len(coord), axis, data.shape[axis]))
        delta_coord = np.diff(coord)
        if not np.allclose(delta_coord, delta_coord[0]):
            raise ValueError('The coordinate must be equally spaced.')
        new_coord = ma.concatenate((coord, coord[-1:] + delta_coord[0]))
    slicer = [slice(None)] * data.ndim
    try:
        slicer[axis] = slice(0, 1)
    except IndexError:
        raise ValueError('The specified axis does not correspond to an '
                         'array dimension.')
    new_data = ma.concatenate((data, data[slicer]), axis=axis)
    if coord is None:
        return_value = new_data
    else:
        return_value = new_data, new_coord
    return return_value
