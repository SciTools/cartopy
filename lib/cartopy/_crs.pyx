# (C) British Crown Copyright 2011 - 2012, Met Office
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


import numpy as np

cimport numpy as np


from cython.operator cimport dereference as deref


cdef extern from "proj_api.h":
    ctypedef void *projPJ
    projPJ pj_init_plus(char *)
    void pj_free(projPJ)
    int pj_transform(projPJ, projPJ, long, int, double *, double *, double *)
    int pj_is_latlong(projPJ)
    char *pj_strerrno(int)
    int *pj_get_errno_ref()
    double DEG_TO_RAD
    double RAD_TO_DEG


cdef double NAN = float('nan')


class Proj4Error(Exception):
    def __init__(self):
        cdef int status
        status = deref(pj_get_errno_ref())
        msg = 'Error from proj.4: {}'.format(pj_strerrno(status))
        self.status = status
        Exception.__init__(self, msg)


cdef class CRS:
    """
    Defines a Coordinate Reference System using proj.4.

    """
    def __init__(self, proj4_params):
        self.proj4_params = dict(proj4_params)
        # Use WGS84 ellipse if one is not specified in proj4_params
        if 'ellps' not in self.proj4_params:
            self.proj4_params['ellps'] = 'WGS84'
        init_items = ['+{}={}'.format(k, v) for
                      k, v in self.proj4_params.iteritems()]
        self.proj4_init = ' '.join(init_items)
        self.proj4 = pj_init_plus(self.proj4_init)
        if not self.proj4:
            raise Proj4Error()

    def __reduce__(self):
        return self.__class__, tuple()
    
    def __getstate__(self):
        return {'proj4_params': self.proj4_params}
    
    def __setstate__(self, state):
        self.__init__(self, **state)
    
    # TODO
    #def __str__
    #def _geod(self): # to return the pyproj.Geod

    def _as_mpl_transform(self, axes=None):
        # XXX This has been replicated in the crs.py Projection class, needs to be consolidated? 
        import cartopy.mpl_integration.geoaxes as geoaxes
        if not isinstance(axes, geoaxes.GenericProjectionAxes):
            raise ValueError('Axes should be an instance of GenericProjectionAxes, got %s' % type(axes))
        return geoaxes.InterProjectionTransform(self, axes.projection) + axes.transData

    property proj4_params:
        def __get__(self):
            return dict(self.proj4_params)

    def as_geocentric(self):
        """
        Returns a new Geocentric CRS with the same ellipse/datum as this
        CRS.

        """
        params = self.proj4_params
        keywords = {}
        for param_name, keyword in {'ellps': 'ellipse', 'datum': 'datum'}.iteritems():
            if param_name in params:
                keywords[keyword] = params[param_name]
        return Geocentric(**keywords)

    def as_geodetic(self):
        """
        Returns a new Geodetic CRS with the same ellipse/datum as this
        CRS.

        """
        params = self.proj4_params
        keywords = {}
        for param_name, keyword in {'ellps': 'ellipse', 'datum': 'datum'}.iteritems():
            if param_name in params:
                keywords[keyword] = params[param_name]
        return Geodetic(**keywords)

    cpdef is_geodetic(self):
        return bool(pj_is_latlong(self.proj4))

    def transform_point(self, double x, double y, CRS src_crs not None):
        cdef:
            double cx, cy
            int status
        cx = x
        cy = y
        if src_crs.is_geodetic():
            cx *= DEG_TO_RAD
            cy *= DEG_TO_RAD
        status = pj_transform(src_crs.proj4, self.proj4, 1, 1, &cx, &cy, NULL);
        #if status == -14 or status == -20:
            # -14 => "latitude or longitude exceeded limits"
            # -20 => "tolerance condition error"
        #    cx = cy = NAN
        #elif status != 0:
        #    raise Proj4Error()
        if self.is_geodetic():
            cx *= RAD_TO_DEG
            cy *= RAD_TO_DEG
        return (cx, cy)

    def transform_points(self, CRS src_crs not None, 
                                np.ndarray[np.double_t, ndim=1] x not None, 
                                np.ndarray[np.double_t, ndim=1] y not None, 
                                np.ndarray[np.double_t, ndim=1] z=None):

        cdef np.ndarray[np.double_t, ndim=2] result
        
        if z is None:
            if x.shape[0] != y.shape[0]:
                raise ValueError('x and y arrays must have the same length')
        elif not x.shape[0] == y.shape[0] == z.shape[0]:
            raise ValueError('x, y, and z arrays must have the same length')

        npts = x.shape[0]

        result = np.empty([npts, 3], dtype=np.double)
        result[:, 0] = x
        result[:, 1] = y
        if z is None:
            result[:, 2] = 0
        else:
            result[:, 2] = z

        if src_crs.is_geodetic():
            result = np.deg2rad(result)

        status = pj_transform(src_crs.proj4, self.proj4, npts, 3,
                              &result[0, 0], &result[0, 1], &result[0, 2]);
                              
        if self.is_geodetic():
            result = np.rad2deg(result)
        #if status:
        #    raise Proj4Error()

        return result


class Geodetic(CRS):
    """
    Defines a latitude/longitude coordinate system with spherical topology
    and geographical distance.

    NB. Coordinates are measured in degrees.

    """
    # XXX Providing a default datum is bad. Providing the ellipse on its own is sufficient to define the ellipse, 
    # and in some cases, can overwrite the desired, well defined ellipse.
    def __init__(self, ellipse='WGS84', datum='WGS84'):
        """
        Create a Geodetic CRS.
        
        Kwargs:
        
            * ellipse      - Ellipsoid definiton.
            * datum        - Datum definiton.
        
        """
        proj4_params = {'proj': 'lonlat', 'ellps': ellipse, 'datum': datum}
        super(Geodetic, self).__init__(proj4_params)        
        
    # XXX Implement fwd such as Basemap's Geod. Would be used in the tissot example.
    # Could come from http://geographiclib.sourceforge.net
      

class Geocentric(CRS):
    """
    Defines a Geocentric coordinate system, where x, y, z are Cartesian coordinates from the center of the Earth.
    
    """
    def __init__(self, ellipse='WGS84', datum='WGS84'):
        proj4_params = {'proj': 'geocent', 'ellps': ellipse, 'datum': datum}
        super(Geocentric, self).__init__(proj4_params)
