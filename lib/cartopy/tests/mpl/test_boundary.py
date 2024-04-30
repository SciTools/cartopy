# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
from itertools import chain

import matplotlib.pyplot as plt
from matplotlib.path import Path
import numpy as np
import pytest

import cartopy.crs as ccrs

circle_verts = np.array([
    (21.33625905034713, 41.90051020408163),
    (20.260167503134653, 36.31721708143521),
    (17.186058185078757, 31.533809612693403),
    (12.554340705353228, 28.235578526280886),
    (7.02857405747471, 26.895042042418392),
    (1.4004024147599825, 27.704250959518802),
    (-3.5238590865749853, 30.547274662874656),
    (-7.038740387110195, 35.0168098237746),
    (-8.701144846177232, 41.42387784534415),
    (-7.802741418346137, 47.03850217758886),
    (-4.224119444187849, 52.606946662776494),
    (0.5099122186645388, 55.75656240544797),
    (6.075429329647158, 56.92110282927732),
    (11.675092899771812, 55.933731058974054),
    (16.50667197715324, 52.93590205611499),
    (19.87797456957319, 48.357097196578444),
    (21.33625905034713, 41.90051020408163)
    ])

rectangle_verts = np.array([
    (55.676020408163225, 36.16071428571428),
    (130.29336734693877, 36.16071428571428),
    (130.29336734693877, -4.017857142857167),
    (55.676020408163225, -4.017857142857167),
    (55.676020408163225, 36.16071428571428)
    ])

@pytest.mark.mpl_image_compare(filename='single_path_boundary_circle.png')
def test_single_path_boundary_circle():
    path = circle_verts
    bnds = [*map(min,zip(*path)), *map(max,zip(*path))]

    f, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
    ax.set_extent((bnds[0], bnds[2], bnds[1], bnds[3]))
    ax.coastlines()
    ax.set_boundary(Path(path))
    return f

@pytest.mark.mpl_image_compare(filename='single_path_boundary_rectangle.png')
def test_single_path_boundary_rectangle():
    path = rectangle_verts
    bnds = [*map(min,zip(*path)), *map(max,zip(*path))]

    f, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
    ax.set_extent((bnds[0], bnds[2], bnds[1], bnds[3]))
    ax.coastlines()
    ax.set_boundary(Path(path))
    return f

@pytest.mark.mpl_image_compare(filename='multi_path_boundary.png')
def test_multi_path_boundary():
    offsets = np.array([[30, 30], [70, 30], [110, 20]])
    circles = [circle_verts + o for o in offsets]

    circle_codes = [
        Path.MOVETO,
        *[Path.LINETO]*(len(circle_verts)),
        Path.CLOSEPOLY
        ]

    rectangle_codes = [
        Path.MOVETO,
        *[Path.LINETO]*(len(rectangle_verts)),
        Path.CLOSEPOLY
        ]

    path = [*chain(*[[c[0], *c, c[-1]] for c in circles]),
            rectangle_verts[0], *rectangle_verts, rectangle_verts[-1]]

    codes = [*chain(*[circle_codes] * len(circles)), *rectangle_codes]

    bnds = [*map(min,zip(*path)), *map(max,zip(*path))]

    f, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
    ax.set_extent((bnds[0], bnds[2], bnds[1], bnds[3]))
    ax.coastlines()
    ax.set_boundary(Path(path, codes))

    return f
