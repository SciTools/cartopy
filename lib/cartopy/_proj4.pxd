# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

"""
This file declares the Proj API, version 4.

"""


cdef extern from "proj_api.h":
    ctypedef void *projPJ
    ctypedef struct projLP:
        double u
        double v

    projPJ pj_init_plus(char *) nogil
    void pj_free(projPJ) nogil
    void pj_get_spheroid_defn(projPJ, double *, double *) nogil
    int pj_transform(projPJ, projPJ, long, int, double *, double *, double *) nogil
    int pj_is_latlong(projPJ) nogil
    char *pj_strerrno(int) nogil
    int *pj_get_errno_ref() nogil
    char *pj_get_release() nogil
    cdef double DEG_TO_RAD
    cdef double RAD_TO_DEG
