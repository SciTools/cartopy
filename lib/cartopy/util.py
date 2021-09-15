# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
"""
This module contains utilities that are useful in conjunction with
cartopy.

"""

import numpy as np
import numpy.ma as ma


def add_cyclic_point(data, coord=None, axis=-1):
    """
    Add a cyclic point to an array and optionally a corresponding
    coordinate.

    Parameters
    ----------
    data
        An n-dimensional array of data to add a cyclic point to.
    coord: optional
        A 1-dimensional array which specifies the coordinate values for
        the dimension the cyclic point is to be added to. The coordinate
        values must be regularly spaced. Defaults to None.
    axis: optional
        Specifies the axis of the data array to add the cyclic point to.
        Defaults to the right-most axis.

    Returns
    -------
    cyclic_data
        The data array with a cyclic point added.
    cyclic_coord
        The coordinate with a cyclic point, only returned if the coord
        keyword was supplied.

    Examples
    --------
    Adding a cyclic point to a data array, where the cyclic dimension is
    the right-most dimension

    >>> import numpy as np
    >>> data = np.ones([5, 6]) * np.arange(6)
    >>> cyclic_data = add_cyclic_point(data)
    >>> print(cyclic_data)  # doctest: +NORMALIZE_WHITESPACE
    [[0. 1. 2. 3. 4. 5. 0.]
     [0. 1. 2. 3. 4. 5. 0.]
     [0. 1. 2. 3. 4. 5. 0.]
     [0. 1. 2. 3. 4. 5. 0.]
     [0. 1. 2. 3. 4. 5. 0.]]

    Adding a cyclic point to a data array and an associated coordinate

    >>> lons = np.arange(0, 360, 60)
    >>> cyclic_data, cyclic_lons = add_cyclic_point(data, coord=lons)
    >>> print(cyclic_data)  # doctest: +NORMALIZE_WHITESPACE
    [[0. 1. 2. 3. 4. 5. 0.]
     [0. 1. 2. 3. 4. 5. 0.]
     [0. 1. 2. 3. 4. 5. 0.]
     [0. 1. 2. 3. 4. 5. 0.]
     [0. 1. 2. 3. 4. 5. 0.]]
    >>> print(cyclic_lons)
    [  0  60 120 180 240 300 360]

    """
    if coord is not None:
        if coord.ndim != 1:
            raise ValueError('The coordinate must be 1-dimensional.')
        if len(coord) != data.shape[axis]:
            raise ValueError(f'The length of the coordinate does not match '
                             f'the size of the corresponding dimension of '
                             f'the data array: len(coord) = {len(coord)}, '
                             f'data.shape[{axis}] = {data.shape[axis]}.')
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
    new_data = ma.concatenate((data, data[tuple(slicer)]), axis=axis)
    if coord is None:
        return_value = new_data
    else:
        return_value = new_data, new_coord
    return return_value
