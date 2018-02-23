# (C) British Crown Copyright 2014 - 2018, Met Office
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
Provide support for converting PROJ.4 strings to Projection instances.

"""

from __future__ import (absolute_import, division, print_function)

import cartopy.crs as ccrs
import shapely.geometry as sgeom

_GLOBE_PARAMS = {'datum': 'datum',
                 'ellps': 'ellipse',
                 'a': 'semimajor_axis',
                 'b': 'semiminor_axis',
                 'f': 'flattening',
                 'rf': 'inverse_flattening',
                 'towgs84': 'towgs84',
                 'nadgrids': 'nadgrids'}
# Map PROJ.4 'proj' parameter to CRS class
PROJ_TO_CRS = {}


def get_proj4_dict(proj4_terms):
    """Convert a PROJ.4 string to a dictionary.

    Parameters
    ----------
    proj4_terms: (str, dict, or iterable of key-value pairs)

    Returns
    -------
    get_proj4_dict
        All PROJ.4 parameters in a dictionary. Any keys with no value
        are set to `None`.

    """
    if isinstance(proj4_terms, dict):
        return proj4_terms
    elif isinstance(proj4_terms, str):
        terms = []
        for term in proj4_terms.split(' '):
            parts = term.strip('+').split('=')
            if len(parts) == 1:
                terms.append((parts[0], None))
            else:
                terms.append(tuple(parts[:2]))
    else:
        # assume list of key value pairs
        terms = proj4_terms

    return dict(terms)


def _split_globe_parameters(proj4_dict):
    projection_terms = {}
    globe_terms = {}
    for name, value in proj4_dict.items():
        if name in _GLOBE_PARAMS:
            globe_terms[name] = value
        else:
            projection_terms[name] = value
    return projection_terms, globe_terms


def _globe_from_proj4(globe_terms):
    """Create a `Globe` object from PROJ.4 parameters."""
    globe = ccrs.Globe(**{_GLOBE_PARAMS[name]: value for name, value in
                          globe_terms.items()})
    return globe


class _PROJ4Projection(ccrs.Projection):
    def __init__(self, proj4_terms, globe=None, bounds=None):
        terms = get_proj4_dict(proj4_terms)
        projection_terms, globe_terms = _split_globe_parameters(terms)
        if globe is None:
            globe = _globe_from_proj4(globe_terms)

        super(_PROJ4Projection, self).__init__(projection_terms, globe)

        # FIXME: Can we guess at the bounds if not provided? Maybe transform
        #        an array of points and take the min/max of the result?
        #        It looks like that's what LambertConformal does.
        self.bounds = bounds

    def __repr__(self):
        return '_PROJ4Projection({})'.format(self.proj4_init)

    @property
    def boundary(self):
        x0, x1, y0, y1 = self.bounds
        return sgeom.LineString([(x0, y0), (x0, y1), (x1, y1), (x1, y0),
                                 (x0, y0)])

    @property
    def x_limits(self):
        x0, x1, y0, y1 = self.bounds
        return (x0, x1)

    @property
    def y_limits(self):
        x0, x1, y0, y1 = self.bounds
        return (y0, y1)

    @property
    def threshold(self):
        x0, x1, y0, y1 = self.bounds
        return min(x1 - x0, y1 - y0) / 100.


def _all_subclasses(cls):
    return cls.__subclasses__() + [g for s in cls.__subclasses__()
                                   for g in _all_subclasses(s)]


def from_proj4(proj4_terms):
    proj4_dict = get_proj4_dict(proj4_terms)

    if not PROJ_TO_CRS:
        # initialize this here instead of at import
        for crs_class in _all_subclasses(ccrs.CRS):
            cls_proj = getattr(crs_class, '_proj4_proj', None)
            if cls_proj is not None and cls_proj not in PROJ_TO_CRS:
                PROJ_TO_CRS[cls_proj] = crs_class

    if 'proj' not in proj4_dict:
        raise ValueError("Missing PROJ.4 parameter: proj")

    proj = proj4_dict['proj']
    crs_class = PROJ_TO_CRS.get(proj)

    # couldn't find a known CRS class
    if crs_class is None:
        # we don't want to allow non-CRS/generic Projection classes
        raise ValueError("Projection '{}' is not implemented yet.".format(
            proj))

    projection_dict, globe_dict = _split_globe_parameters(proj4_dict)
    globe = _globe_from_proj4(globe_dict)
    return crs_class.from_proj4(projection_dict, globe=globe)
