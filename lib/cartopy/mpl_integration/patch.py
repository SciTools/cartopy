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


"""
Provides geos <-> mpl path support.

Please note: This code needs a significant spruce and review.
  
"""
import shapely.geometry as shpgeom

import numpy as np
from matplotlib.path import Path


def geos_to_path(shape):
    """Return a list of paths created from this shape."""
    if isinstance(shape, (list, tuple)):
        r = []
        for shp in shape:
            r.extend(geos_to_path(shp))
        return r 
    if isinstance(shape, (shpgeom.linestring.LineString, shpgeom.point.Point)):
        return [Path(np.vstack(shape.xy).T)]
    elif isinstance(shape, (shpgeom.multipolygon.MultiPolygon)):
        r = []
        for shp in shape:
            r.extend(geos_to_path(shp))
        return r
        
    elif isinstance(shape, (shpgeom.polygon.Polygon)):
        def poly_codes(poly):
            r = np.ones(len(poly.xy[0])) * Path.LINETO
            r[0] = Path.MOVETO
            return r
        
        vertices = np.concatenate([np.array(shape.exterior.xy)] +
                                  [np.array(ring.xy) for ring in shape.interiors], 1).T
        codes = np.concatenate(
                [poly_codes(shape.exterior)]
                + [poly_codes(ring) for ring in shape.interiors])
        return [Path(vertices, codes)]
    elif isinstance(shape, (MultiPolygon, GeometryCollection, MultiLineString, MultiPoint)):
        r = []
        for geom in shape.geoms:
            r.extend(geos_to_path(geom))
        return r
    elif isinstance(shape, (GeometryCollection, MultiLineString)):
        print type(shape)
        return [Path(np.vstack(line.xy).T) for line in shape]
    elif hasattr(shape, '_as_mpl_path'):
        verts, codes = shape._as_mpl_path()
        return [Path(verts, codes)]
    else:
        raise ValueError('Unexpected GEOS type. Got %s' % type(shape))




from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.cm


from matplotlib.path import Path
from shapely.geometry.collection import GeometryCollection
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.multilinestring import MultiLineString
from shapely.geometry import polygon, linestring, point
from shapely.geometry.multipoint import MultiPoint




def path_segments(path, transform=None, remove_nans=False, clip=None,
                      quantize=False, simplify=False, curves=False, 
                      stroke_width=1.0, snap=False):
    """See path.iter_segments. Vectorised version of that."""
    # XXX assigned to avoid a ValueError inside the mpl C code...
    a = transform, remove_nans, clip, quantize, simplify, curves
    
    vertices, codes = mpath.cleanup_path(path, transform, remove_nans, clip,
                                       snap, stroke_width, simplify, curves)
    
    # cut the final vertices which end with 0
    return vertices[:-1, :], codes[:-1]
    
    
#    print vertices, codes
#    for v, c in path.iter_segments(transform, simplify=simplify, curves=curves):
#        print 'iter: ', v, c
#    vertices = vertices[np.where(codes != 0)[0], :]
#    codes = codes[np.where(codes != 0)[0]]
    return vertices, codes
#    return path.get_verts_codes(transform, remove_nans, clip,
#                                         quantize, simplify, curves)





def path_interpolation(path, n_steps):
#            path_verts, path_codes =  zip(*list(path.iter_segments(curves=False)))
    path_verts, path_codes = path_segments(path, curves=False)
    path_verts = np.array(path_verts)
    path_codes = np.array(path_codes)
    verts_split_inds = np.where(path_codes == Path.MOVETO)[0]
    verts_split = np.split(path_verts, verts_split_inds, 0)
    codes_split = np.split(path_codes, verts_split_inds, 0)
    
    v_collection = []
    c_collection = []
    for path_verts, path_codes in zip(verts_split, codes_split):
        if len(path_verts) == 0:
            continue
        import matplotlib.cbook
#                print path_verts.shape
        verts = matplotlib.cbook.simple_linear_interpolation(path_verts, n_steps)
        v_collection.append(verts)
#                print verts.shape
        codes = np.ones(verts.shape[0]) * Path.LINETO
        codes[0] = Path.MOVETO
        c_collection.append(codes)
        
    return Path(np.concatenate(v_collection), np.concatenate(c_collection))


def path_to_geos(path):
    """
    """
    import time
    start_time = time.time()
    log = []
    
    DEBUG = False
#    path_verts, path_codes =  zip(*list(path.iter_segments(curves=False)))
    path_verts, path_codes = path_segments(path, curves=False)
    path_verts = np.array(path_verts)
    path_codes = np.array(path_codes)
    
    if DEBUG: print 'codes:', path_codes
    verts_split_inds = np.where(path_codes == Path.MOVETO)[0]
    verts_split = np.split(path_verts, verts_split_inds, 0)
    codes_split = np.split(path_codes, verts_split_inds, 0)
    
    if DEBUG: print 'vs: ', `verts_split`
    if DEBUG: print 'cs: ', `codes_split`
    
    log.append('split done %s' % (time.time() - start_time))
    
    collection = []
    for path_verts, path_codes in zip(verts_split, codes_split):
        if len(path_verts) == 0:
            continue
        # XXX A path can be given which does not end with close poly, in that situation, we have to guess?
        if DEBUG: print 'pv: ', path_verts
        # XXX Implement a point
                
        if path_verts.shape[0] > 2 and (path_codes[-1] == Path.CLOSEPOLY or all(path_verts[0, :] == path_verts[-1, :])):
            if path_codes[-1] == Path.CLOSEPOLY:
                ipath2 = polygon.Polygon(path_verts[:-1, :])
            else:
                ipath2 = polygon.Polygon(path_verts)
        else:
            ipath2 = linestring.LineString(path_verts)
            
        if (len(collection) > 0 and 
                 isinstance(collection[-1][0], polygon.Polygon) and
                 isinstance(ipath2, polygon.Polygon) and
                 collection[-1][0].contains(ipath2.exterior)):
            collection[-1][1].append(ipath2.exterior)
        else:
            # collection is a list of [exernal_poly, list_of_internal_polys]
            collection.append([ipath2, []])
    
    log.append('collection done before while %s.' % (time.time() - start_time))
    
    log.append('Len of collection %s.' % (len(collection)))

    res = []
    
    for external_poly, internal_polys in collection:
#        print external_poly
        if len(internal_polys) > 0:
#            print internal_polys
            # XXX worry about islands within lakes
            poly = polygon.Polygon(external_poly.exterior, internal_polys)
        else:
            poly = external_poly
        res.append(poly)
    collection = res
#    if len(collection) > 1:
#        i = 0
#        while i<len(collection)-1:
#            poly = collection[i]
#            poly2 = collection[i+1]
#                        
#            # TODO Worry about islands within lakes
#            if isinstance(poly, polygon.Polygon) and isinstance(poly2, polygon.Polygon):  
#                if poly.contains(poly2):
#                    # XXX This is the slow bit!
##                    collection[i] = polygon.Polygon(poly.exterior, list(poly.interiors) + [poly2.exterior])                    
#                    collection.pop(i+1)
#                    continue
#            else:
#                res.append([poly])
#            i+=1
#    
#        log.append('Post len of collection %s.' % (len(collection)))
#        log.append('collection done after while %s' % (time.time() - start_time))
        
    if len(collection) == 1:
        result = collection
    else:
        if all([isinstance(geom, linestring.LineString) for geom in collection]):
            result = [MultiLineString(collection)]
        else:
            result = collection
#            if DEBUG: print 'geom: ', collection, type(collection)
#            raise NotImplementedError('The path given was not a collection of line strings, ' 
#                                      'nor a single polygon with interiors.')
    if (time.time() - start_time) > 1:
        print 'geos time %s' % (time.time() - start_time)
        print '\n'.join(log)
    
    return result
