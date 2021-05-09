# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from unittest import mock

import numpy as np
import pytest
import shapely.geometry as sgeom
from matplotlib.transforms import IdentityTransform

import cartopy.crs as ccrs
import cartopy.mpl.geoaxes as geoaxes
from cartopy.feature import ShapelyFeature
from cartopy.mpl.feature_artist import FeatureArtist, _freeze, _GeomKey
from cartopy.mpl import style


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
    unit_circle = sgeom.Point(0, 0).buffer(0.5)
    unit_square = unit_circle.envelope
    geoms = [unit_circle, unit_square]
    feature = ShapelyFeature(geoms, ccrs.PlateCarree())
    return feature


def mocked_axes(extent, projection=ccrs.PlateCarree()):
    return mock.MagicMock(
        get_extent=mock.Mock(return_value=extent),
        projection=projection,
        spec=geoaxes.GeoAxes,
        transData=IdentityTransform(),
        patch=mock.sentinel.patch,
        figure=mock.sentinel.figure)


# Need to initialize private renderer properties on the sentinal
mock.sentinel.renderer._raster_depth = 0
mock.sentinel.renderer._rasterizing = False


def style_from_call(call):
    args, kwargs = call
    # Drop the transform keyword.
    kwargs.pop('transform')
    return kwargs


def cached_paths(geom, target_projection):
    # Use the cache in FeatureArtist to get back the projected path
    # for the given geometry.
    geom_cache = FeatureArtist._geom_key_to_path_cache.get(_GeomKey(geom), {})
    return geom_cache.get(target_projection, None)


@mock.patch('matplotlib.collections.PathCollection')
def test_feature_artist_draw(path_collection_cls, feature):
    geoms = list(feature.geometries())

    fa = FeatureArtist(feature, facecolor='red')
    prj_crs = ccrs.Robinson()
    fa.axes = mocked_axes(extent=[-10, 10, -10, 10], projection=prj_crs)
    fa.draw(mock.sentinel.renderer)

    transform = prj_crs._as_mpl_transform(fa.axes)
    expected_paths = (cached_paths(geoms[0], prj_crs) +
                      cached_paths(geoms[1], prj_crs), )
    expected_style = {'facecolor': 'red'}

    args, kwargs = path_collection_cls.call_args_list[0]
    assert transform == kwargs.pop('transform', None)
    assert kwargs == expected_style
    assert args == expected_paths

    path_collection_cls().set_clip_path.assert_called_once_with(fa.axes.patch)
    path_collection_cls().set_figure.assert_called_once_with(fa.axes.figure)
    path_collection_cls().draw(mock.sentinel.renderer)


@mock.patch('matplotlib.collections.PathCollection')
def test_feature_artist_draw_styler(path_collection_cls, feature):
    geoms = list(feature.geometries())
    style1 = {'facecolor': 'blue', 'edgecolor': 'white'}
    style2 = {'color': 'black', 'linewidth': 1}
    style2_finalized = style.finalize(style.merge(style2))

    def styler(geom):
        if geom == geoms[0]:
            return style1
        else:
            return style2

    fa = FeatureArtist(feature, styler=styler, linewidth=2)
    extent = [-10, 10, -10, 10]
    prj_crs = ccrs.Robinson()
    fa.axes = mocked_axes(extent, projection=prj_crs)
    fa.draw(mock.sentinel.renderer)

    transform = prj_crs._as_mpl_transform(fa.axes)

    calls = [{'paths': (cached_paths(geoms[0], prj_crs), ),
              'style': dict(linewidth=2, **style1)},
             {'paths': (cached_paths(geoms[1], prj_crs), ),
              'style': style2_finalized}]

    assert path_collection_cls.call_count == 2
    for expected_call, (actual_args, actual_kwargs) in \
            zip(calls, path_collection_cls.call_args_list):
        assert expected_call['paths'] == actual_args
        assert transform == actual_kwargs.pop('transform')
        assert expected_call['style'] == actual_kwargs
