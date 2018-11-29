# (C) British Crown Copyright 2018, Met Office
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
This file declares the Proj API, version 4.

"""


cdef extern from "proj_api.h":
    ctypedef void *projPJ
    projPJ pj_init_plus(char *) nogil
    void pj_free(projPJ) nogil
    int pj_transform(projPJ, projPJ, long, int, double *, double *, double *) nogil
    int pj_is_latlong(projPJ) nogil
    char *pj_strerrno(int) nogil
    int *pj_get_errno_ref() nogil
    char *pj_get_release() nogil
    cdef double DEG_TO_RAD
    cdef double RAD_TO_DEG
