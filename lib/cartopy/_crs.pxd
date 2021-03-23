# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.


cdef class CRS:
    """
    Defines a Coordinate Reference System using proj.

    """

    cdef readonly proj4_init
    cdef proj4_params

    cpdef is_geodetic(self)
