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


def _xml_children(element, name=None):
    for child in element.childNodes:
        if name is None or (hasattr(child, "tagName") and name in child.tagName):
            yield child


def _xml_grand_children(feature, name=None):
    for feature_child in _xml_children(feature):
        for child in _xml_children(feature_child, name):
            yield child


def _split_coord_list(coord_list):
    # Convert gml coord list into sequence of (x,y) tuples.
    coords = coord_list.split()
    coords = [(i.split(",")) for i in coords]
    x = [float(i[0]) for i in coords]
    y = [float(i[1]) for i in coords]
    return zip(x, y)


def _gml_to_shapely(gml):

    import xml.dom.minidom
    import shapely.geometry

    doc = xml.dom.minidom.parseString(gml)
    geoms = []

    feature_collection = _xml_children(doc, "wfs:FeatureCollection").next()
    if feature_collection:
        
        # All feature members.
        for feature_member in _xml_children(feature_collection, "gml:featureMember"):
            # Probably only one child, but loop just in case.
            for feature in _xml_children(feature_member):
            
                try:
                    name = _xml_children(feature, "name").next()
                    name = name.firstChild.data
                except:
                    name = "(no name)"
                
                # get anything with a geometry child
                for linestring in _xml_grand_children(feature, "gml:LineString"):
                    # Probably only one 'coordinates' child, but loop just in case.
                    for coords in _xml_children(linestring, "gml:coordinates"):
                        coord_list = coords.firstChild.data
                        xy = _split_coord_list(coord_list)
                        geoms.append(shapely.geometry.LineString(xy))
    
                for point in _xml_grand_children(feature, "gml:Point"):
                    # Probably only one 'coordinates' child, but loop just in case.
                    for coords in _xml_children(point, "gml:coordinates"):
                        coords = coords.firstChild.data
                        coords = coords.split(",")
                        x, y = float(coords[0]), float(coords[1])
                        geoms.append(shapely.geometry.Point(x, y))
    
                for multilinestring in _xml_grand_children(feature, "gml:MultiLineString"):
                    for linestring in _xml_grand_children(multilinestring, "gml:LineString"):
                        for coords in _xml_children(linestring, "gml:coordinates"):
                            coord_list = coords.firstChild.data
                            xy = _split_coord_list(coord_list)
                            geoms.append(shapely.geometry.LineString(xy))
    
                for multilinestring in _xml_grand_children(feature, "gml:MultiPolygon"):
                    for polygon in _xml_grand_children(multilinestring, "gml:Polygon"):
                        for linearring in _xml_grand_children(polygon, "gml:LinearRing"):
                            for coords in _xml_children(linearring, "gml:coordinates"):
                                coord_list = coords.firstChild.data
                                xy = _split_coord_list(coord_list)
                                geoms.append(shapely.geometry.Polygon(xy))

    else:
        raise ValueError("No FeatureCollection")

    return geoms


def wfs(url, version=None, features=None, max_features=None, bbox=None, srs=None):
    """
    Retrieve WFS features.
    
    Returns shapely geometry.
    
    """
    import urllib

    # Check params and construct request string.    
    request = url
    if not request.endswith("?"):
        request += "?"
        
    if version is None:
        version = "1.0.0"
        
    if features is None:
        raise ValueError("No features were requested")
        
    request += "service=WFS&version={}&request=getFeature&typename={}".format(version, features)

    if max_features is not None:
        request += "&MAXFEATURES={}".format(max_features)
    
    if bbox is not None:
        request += "&bbox={}".format(bbox)

    if srs is not None:
        request += "&srs={}".format(srs)
    
    # Get the gml and turn it into shapely geometry.
    gml = urllib.urlopen(request).read()
    return _gml_to_shapely(gml)


if __name__ == "__main__":
    
    url = "http://exxvmgpmdev0.meto.gov.uk:8399/arcgis/services/JM_test/DMMS_test_global/MapServer/WFSServer?"
    version = "1.0.0"
    features = "DMMS_test_global:Global_land_polygons"

    geoms = wfs(url, version, features, srs="EPSG:4326")

    import matplotlib.pyplot as plt    
    import cartopy.crs as ccrs
    
    plt.axes(projection=ccrs.PlateCarree())
    plt.gca().add_geometries(geoms, ccrs.PlateCarree())
    plt.show()
    
    