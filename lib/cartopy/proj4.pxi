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

cdef extern from "proj.h":

    cdef enum:
        PROJ_VERSION_MAJOR
        PROJ_VERSION_MINOR
        PROJ_VERSION_PATCH

    # projPJ  has been replaced by PJ *
    ctypedef struct PJ
    ctypedef struct PJ_CONTEXT
    int proj_errno(const PJ *P)
    int proj_context_errno (PJ_CONTEXT *ctx)
    const char * proj_errno_string (int err)
    PJ *proj_create (PJ_CONTEXT *ctx, const char *definition)
    void proj_destroy(PJ *obj)

    cdef struct PJ_PROJ_INFO:
        const char  *id
        const char  *description
        const char  *definition
        int         has_inverse #1 if an inverse mapping exists, 0 otherwise              */
        double      accuracy

    PJ_PROJ_INFO proj_pj_info(PJ *P)

    PJ *proj_crs_get_geodetic_crs(PJ_CONTEXT *ctx, const PJ *crs)

    ctypedef struct PJ_XYZT:
        double   x,   y,  z, t
    ctypedef struct PJ_UVWT:
        double   u,   v,  w, t
    ctypedef struct PJ_LPZT:
        double lam, phi,  z, t
    ctypedef struct PJ_OPK:
        double o, p, k
    ctypedef struct PJ_ENU:
        double e, n, u
    ctypedef struct PJ_GEOD:
        double s, a1, a2

    ctypedef struct PJ_UV:
        double   u,   v
    ctypedef struct PJ_XY:
        double   x,   y
    ctypedef struct PJ_LP:
        double lam, phi

    ctypedef struct PJ_XYZ:
        double   x,   y,  z
    ctypedef struct PJ_UVW:
        double   u,   v,  w
    ctypedef struct PJ_LPZ:
        double lam, phi,  z


    cdef union PJ_COORD:
        double v[4];
        PJ_XYZT xyzt;
        PJ_UVWT uvwt;
        PJ_LPZT lpzt;
        PJ_GEOD geod;
        PJ_OPK opk;
        PJ_ENU enu;
        PJ_XYZ xyz;
        PJ_UVW uvw;
        PJ_LPZ lpz;
        PJ_XY xy;
        PJ_UV uv;
        PJ_LP lp;

    cdef enum PJ_DIRECTION:
        PJ_FWD   =  1 # Forward
        PJ_IDENT =  0 # Do nothing
        PJ_INV   = -1 # Inverse

    PJ_COORD proj_trans (PJ *P, PJ_DIRECTION direction, PJ_COORD coord)
    size_t proj_trans_generic (
        PJ *P,
        PJ_DIRECTION direction,
        double *x, size_t sx, size_t nx,
        double *y, size_t sy, size_t ny,
        double *z, size_t sz, size_t nz,
        double *t, size_t st, size_t nt
    );
    ctypedef struct PJ_AREA
    PJ *proj_create_crs_to_crs(PJ_CONTEXT *ctx, const char *source_crs, const char *target_crs, PJ_AREA *area);

    ctypedef enum PJ_TYPE:
        PJ_TYPE_UNKNOWN
        PJ_TYPE_ELLIPSOID
        PJ_TYPE_GEODETIC_REFERENCE_FRAME
        PJ_TYPE_DYNAMIC_GEODETIC_REFERENCE_FRAME
        PJ_TYPE_VERTICAL_REFERENCE_FRAME
        PJ_TYPE_DYNAMIC_VERTICAL_REFERENCE_FRAME
        PJ_TYPE_DATUM_ENSEMBLE
        PJ_TYPE_CRS
        PJ_TYPE_GEODETIC_CRS
        PJ_TYPE_GEOCENTRIC_CRS
        PJ_TYPE_GEOGRAPHIC_CRS
        PJ_TYPE_GEOGRAPHIC_2D_CRS
        PJ_TYPE_GEOGRAPHIC_3D_CRS
        PJ_TYPE_VERTICAL_CRS
        PJ_TYPE_PROJECTED_CRS
        PJ_TYPE_COMPOUND_CRS
        PJ_TYPE_TEMPORAL_CRS
        PJ_TYPE_ENGINEERING_CRS
        PJ_TYPE_BOUND_CRS
        PJ_TYPE_OTHER_CRS
        PJ_TYPE_CONVERSION
        PJ_TYPE_TRANSFORMATION
        PJ_TYPE_CONCATENATED_OPERATION
        PJ_TYPE_OTHER_COORDINATE_OPERATION

    PJ_TYPE proj_get_type(const PJ *obj);

    PJ *proj_get_ellipsoid(PJ_CONTEXT *ctx, const PJ *obj)
    int proj_ellipsoid_get_parameters(PJ_CONTEXT *ctx,
                                      const PJ *ellipsoid,
                                      double *out_semi_major_metre,
                                      double *out_semi_minor_metre,
                                      int    *out_is_semi_minor_computed,
                                      double *out_inv_flattening);

    double proj_torad (double angle_in_degrees)
    double proj_todeg (double angle_in_radians)

    ctypedef enum PJ_LOG_LEVEL:
        PJ_LOG_NONE  = 0
        PJ_LOG_ERROR = 1
        PJ_LOG_DEBUG = 2
        PJ_LOG_TRACE = 3
        PJ_LOG_TELL  = 4
    ctypedef void (*PJ_LOG_FUNCTION)(void *, int, const char *)
    void proj_log_func (PJ_CONTEXT *ctx, void *app_data, PJ_LOG_FUNCTION logf)
