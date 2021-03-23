# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

"""
This file declares the Proj API, version 4.

"""


cdef extern from "proj.h":
    ctypedef void *PJ_CONTEXT;
    ctypedef void *PJ;
    ctypedef void *PJ_AREA;
    ctypedef struct PJ_INFO:
        int     major
        int     minor
        int     patch
        char   *release
        char   *version
        char   *searchpath
    ctypedef enum PJ_DIRECTION:
        PJ_FWD   =  1
        PJ_IDENT =  0
        PJ_INV   = -1

    PJ *proj_create(PJ_CONTEXT *, const char *) nogil;
    PJ *proj_create_crs_to_crs(PJ_CONTEXT *, const char *, const char *, PJ_AREA *) nogil;
    PJ *proj_destroy(PJ *) nogil;
    size_t proj_trans_generic(PJ *, PJ_DIRECTION, double *, size_t, size_t, double *, size_t, size_t, double *, size_t, size_t, double *, size_t, size_t) nogil;
    int proj_angular_output(PJ *, PJ_DIRECTION) nogil;
    int proj_context_errno(PJ_CONTEXT *) nogil;
    int proj_errno(PJ *) nogil;
    char *proj_errno_string(int) nogil;
    PJ_INFO proj_info() nogil;
