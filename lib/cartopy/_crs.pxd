# (C) British Crown Copyright 2010 - 2012, Met Office
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


cdef extern from "proj_api.h":
    ctypedef void *projPJ


cdef class CRS:
    """
    Defines a Coordinate Reference System using proj.4.

    """

    cdef projPJ proj4
    cdef readonly proj4_init
    cdef proj4_params

    cpdef is_geodetic(self)
