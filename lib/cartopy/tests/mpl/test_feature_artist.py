# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

from unittest import mock

import matplotlib.colors as mcolors
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import pyproj
import pytest
import shapely

import cartopy.crs as ccrs
from cartopy.feature import Feature, ShapelyFeature
from cartopy.mpl.feature_artist import (
    _PARTS_KEY,
    FeatureArtist,
    _freeze,
    _GeomKey,
)
from cartopy.mpl.path import shapely_to_path


@pytest.mark.parametrize("source, expected", [
    [{1: 0}, frozenset({(1, 0)})],
    [[1, 2], (1, 2)],
    [[1, {}], (1, frozenset())],
    [[1, {'a': [1, 2, 3]}], (1, frozenset([('a', (1, 2, 3))]))],
    [{'edgecolor': 'face', 'zorder': -1,
      'facecolor': np.array([0.9375, 0.9375, 0.859375])},
     frozenset([('edgecolor', 'face'), ('zorder', -1),
                ('facecolor', (0.9375, 0.9375, 0.859375))])],
])
def test_freeze(source, expected):
    assert _freeze(source) == expected


@pytest.fixture
def feature():
    circle1 = shapely.Point(0, 0).buffer(1)
    circle2 = shapely.Point(0, 0).buffer(10)
    square = shapely.Polygon([(30, 0), (50, 0), (50, 20), (30, 20), (30, 0)])
    geoms = [circle1, circle2, square]
    feature = ShapelyFeature(geoms, ccrs.PlateCarree())
    return feature


def cached_paths(geom, target_projection):
    # Use the cache in FeatureArtist to get back the projected path
    # for the given geometry.
    geom_cache = FeatureArtist._geom_key_to_path_cache.get(_GeomKey(geom), {})
    return geom_cache.get(target_projection, None)


@pytest.mark.natural_earth
@pytest.mark.mpl_image_compare(filename='feature_artist.png')
@pytest.mark.parametrize(
    'method',
    ['default', 'facecolor_list', 'cmap', 'styled_feature', 'styler'])
def test_feature_artist(feature, method):
    # Set up a common map for the image tests.
    # The extent is chosen to include only the square geometry from `feature`.
    # This means that we can check that `array` or a list of facecolors remains 1-to-1
    # with the list of geometries.
    prj_crs = ccrs.Robinson()
    fig, ax = plt.subplots(subplot_kw={'projection': prj_crs})
    ax.set_extent([20, 180, -90, 90])
    ax.coastlines()

    match method:
        case 'default':
            ax.add_feature(feature, facecolor='blue')

        case 'facecolor_list':
            ax.add_feature(feature, facecolor=['red', 'green', 'blue'])

        case 'cmap':
            cmap = mcolors.ListedColormap(['red', 'gray', 'blue'])
            ax.add_feature(feature, cmap=cmap, array=[0, 0, 1])

        case 'styled_feature':
            geoms = list(feature.geometries())
            styled_feature = ShapelyFeature(geoms, crs=ccrs.PlateCarree(),
                                            facecolor='blue')
            ax.add_feature(styled_feature)

        case 'styler':
            geoms = list(feature.geometries())

            def styler(geom):
                if geom == geoms[1]:
                    return {'facecolor': 'red'}
                else:
                    return {'facecolor': 'blue'}

            ax.add_feature(feature, facecolor='grey', styler=styler)

        case _:
            raise ValueError(f'Unknown feature artist draw method {method!r}')

    return fig


def test_feature_artist_geom_single_path(feature):
    plot_crs = ccrs.PlateCarree(central_longitude=180)
    fig, ax = plt.subplots(subplot_kw={'projection': plot_crs})
    ax.add_feature(feature)

    fig.draw_without_rendering()

    # Circles get split into two geometries across the dateline, but should still be
    # plotted as one compound path to ensure style consistency.
    for geom in feature.geometries():
        assert isinstance(cached_paths(geom, plot_crs), mpath.Path)


@pytest.fixture
def multipart_feature():
    # One MultiPolygon with parts spread around the globe, like the records
    # of the Natural Earth land dataset.
    squares = [shapely.box(x, 0, x + 10, 10) for x in (-170, 0, 30, 160)]
    return ShapelyFeature([shapely.MultiPolygon(squares)], ccrs.PlateCarree())


def test_feature_artist_multipart_projects_visible_parts(multipart_feature):
    plot_crs = ccrs.LambertConformal(central_longitude=15)
    fig, ax = plt.subplots(subplot_kw={'projection': plot_crs})
    ax.set_extent([-10, 50, -5, 20], crs=ccrs.PlateCarree())
    artist = ax.add_feature(multipart_feature)

    fig.draw_without_rendering()

    [geom] = multipart_feature.geometries()
    [(yielded_geom, path)] = list(artist._get_geoms_paths())
    # The styler and array handling still see the original geometry.
    assert yielded_geom is geom
    # Only the parts in view are projected, the whole geometry is not.
    assert cached_paths(geom, plot_crs) is None
    parts = FeatureArtist._geom_key_to_path_cache[_GeomKey(geom)][_PARTS_KEY]
    projected = [cached_paths(part, plot_crs) is not None for part in parts]
    assert projected == [False, True, True, False]

    # The combined path is the same as projecting the visible parts together.
    expected = shapely_to_path(plot_crs.project_geometry(
        shapely.MultiPolygon(list(parts[1:3])), ccrs.PlateCarree()))
    np.testing.assert_array_equal(path.vertices, expected.vertices)
    np.testing.assert_array_equal(path.codes, expected.codes)


def test_feature_artist_multipart_in_view_not_split(multipart_feature):
    plot_crs = ccrs.Robinson()
    fig, ax = plt.subplots(subplot_kw={'projection': plot_crs})
    ax.set_global()
    ax.add_feature(multipart_feature)

    fig.draw_without_rendering()

    [geom] = multipart_feature.geometries()
    assert isinstance(cached_paths(geom, plot_crs), mpath.Path)
    assert _PARTS_KEY not in FeatureArtist._geom_key_to_path_cache[_GeomKey(geom)]


def test_feature_artist_multipart_styler(multipart_feature):
    seen = []

    def styler(geom):
        seen.append(geom)
        return {'facecolor': 'red'}

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-10, 50, -5, 20])
    ax.add_feature(multipart_feature, styler=styler)

    fig.draw_without_rendering()

    assert seen == list(multipart_feature.geometries())


# The contents of the .prj files of the Natural Earth shapefiles.
NATURAL_EARTH_WKT = (
    'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'
    'SPHEROID["WGS_1984",6378137.0,298.257223563]],'
    'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')


class CRSOnReadFeature(Feature):
    """
    A feature that replaces its CRS with an equivalent one when reading its
    geometries, like NaturalEarthFeature does with the CRS of the shapefile.

    """

    def __init__(self, geoms):
        super().__init__(ccrs.PlateCarree())
        self._geoms = geoms

    def geometries(self):
        self._crs = ccrs.Projection(pyproj.CRS.from_wkt(NATURAL_EARTH_WKT))
        return iter(self._geoms)


def test_feature_artist_crs_changed_by_reading_geometries(multipart_feature):
    # Geometries are read lazily while drawing. The CRS of the feature from
    # before reading them has to be used, so that no projection is needed
    # if it matches the projection of the axes.
    # The CRS after reading is equivalent, but does not compare equal.
    assert ccrs.PlateCarree() != ccrs.Projection(
        pyproj.CRS.from_wkt(NATURAL_EARTH_WKT))
    feature = CRSOnReadFeature(list(multipart_feature.geometries()))

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-10, 50, -5, 20])
    ax.add_feature(feature)

    with mock.patch.object(ccrs.PlateCarree, 'project_geometry') as project:
        fig.draw_without_rendering()
    project.assert_not_called()


@pytest.mark.parametrize('autolim', [False, True])
def test_feature_artist_autolim(autolim):
    plot_crs = ccrs.PlateCarree(central_longitude=180)
    fig, ax = plt.subplots(subplot_kw={'projection': plot_crs})

    square = shapely.Polygon([(30, 0), (50, 0), (50, 20), (30, 20), (30, 0)])
    ax.add_geometries([square], crs=ccrs.PlateCarree(), autolim=autolim)

    if autolim:
        expected = [-150, 0, 20, 20]  # Fit to square projected 180 degrees
    else:
        expected = [-180, -90, 360, 180]  # Leave at default of whole globe

    np.testing.assert_allclose(ax.dataLim.bounds, expected)
