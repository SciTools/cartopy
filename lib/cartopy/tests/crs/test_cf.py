# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
"""
Tests for creating CRS from CF attributes
"""
import numpy as np
import pytest

import cartopy.crs as ccrs


# All taken from examples of the CF-Convenvtions, section 5.6
CRS = [
    (
        {'grid_mapping_name': 'latitude_longitude'},
        ccrs.PlateCarree()
    ),
        (
        {
            'grid_mapping_name': 'latitude_longitude',
            'semi_major_axis': 6371000.0,
            'inverse_flattening': 0,
        },
        ccrs.PlateCarree(globe=ccrs.Globe(
            semimajor_axis=6371000.0,
            ellipse='sphere'
        ))
    ),
    (
        {
            'grid_mapping_name': 'rotated_latitude_longitude',
            'grid_north_pole_latitude': 36.0,
            'grid_north_pole_longitude': 74.0
        },
        ccrs.RotatedPole(pole_latitude=36, pole_longitude=74)
    ),
    (
        {
            'grid_mapping_name': 'lambert_conformal_conic',
            'standard_parallel': 25.0,
            'longitude_of_central_meridian': 265.0,
            'latitude_of_projection_origin': 25.0
        },
        ccrs.LambertConformal(
            central_longitude=265.0,
            central_latitude=25.0,
            standard_parallels=(25,)
        )
    ),
]


@pytest.mark.parametrize("cfattrs,exp", CRS)
def test_from_cf(cfattrs, exp):
    crs = ccrs.from_cf(**cfattrs)
    assert crs == exp


def test_from_cf_OSGB():
    # These attrs are from the cf-conventions examples
    # We can't compare the TransverseMercator created with those directly with OSGB
    # So we compare values from a transform
    crs = ccrs.from_cf(**
        {
            "grid_mapping_name": "transverse_mercator",
            "semi_major_axis": 6377563.396,
            "inverse_flattening": 299.3249646,
            "longitude_of_prime_meridian": 0.0,
            "latitude_of_projection_origin": 49.0,
            "longitude_of_central_meridian": -2.0,
            "scale_factor_at_central_meridian": 0.9996012717,
            "false_easting": 400000.0,
            "false_northing": -100000.0,
            "unit": "metre",
        }
    )
    exp = ccrs.OSGB()

    # Arbitrary points, all within the british isles
    ptsLon = np.array([-4.256495470173718, -0.23201298137905724, -4.5156834308392035])
    ptsLat = np.array([58.3199544202207960, 52.59337104346588, 50.3827602907713])

    cfpts = crs.transform_points(ccrs.PlateCarree(), ptsLon, ptsLat)
    exppts = exp.transform_points(ccrs.PlateCarree(), ptsLon, ptsLat)

    np.testing.assert_array_equal(cfpts, exppts)
