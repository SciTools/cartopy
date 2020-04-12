# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from ._proj4 cimport projPJ


cdef class CRS:
    """
    Defines a Coordinate Reference System using proj.

    """

    cdef projPJ proj4
    cdef readonly proj4_init
    cdef proj4_params

    cpdef is_geodetic(self)
