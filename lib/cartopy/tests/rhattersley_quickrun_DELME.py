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


import itertools
import math
import os.path
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import shapely.geometry as sgeom

import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader

COASTLINE_PATH = shapereader.natural_earth(name='coastline')
LAND_PATH = shapereader.natural_earth(name='land')


def _arrows(projection, geometry):
    coords = list(geometry.coords)
    for i, xyxy in enumerate(list(zip(coords[:-1], coords[1:]))):
        if i % 11 == 0:
            start, end = xyxy
            dx = end[0] - start[0]
            dy = end[1] - start[1]
            s = 4
            #width = math.sqrt(dx*dx+dy*dy)*0.5*s
            #arrow = mpatches.Arrow(start[0], start[1], dx*s,
            #                       dy*s, width=width, alpha=0.4)
            mag = math.sqrt(dx * dx + dy * dy)
            length = projection.threshold * s
            try:
                dx = length * dx / mag
                dy = length * dy / mag
                width = projection.threshold * s * 0.5
                arrow = mpatches.Arrow(start[0], start[1], dx, dy,
                                       width=width, alpha=0.4)
                plt.gca().add_patch(arrow)
            except ZeroDivisionError:
                pass


def draw_line_string(projection, line_string, color='black', linestyle='-'):
    multi_line_string = projection.project_geometry(line_string)
    for line_string in multi_line_string:
        plt.plot(*list(zip(*line_string.coords)),
                 marker='', color=color, linestyle=linestyle)
        #_arrows(projection, line_string)


def draw_polygon(projection, polygon, color=None):
    multi_polygon = projection.project_geometry(polygon)
    for polygon in multi_polygon:
        #plt.plot(*zip(*polygon.exterior.coords), marker='+', color=color)
        #_arrows(projection, polygon.exterior)
        #continue
        import cartopy.mpl.patch as patch
        paths = patch.geos_to_path(polygon)
        for pth in paths:
            patch = mpatches.PathPatch(pth, edgecolor='none',
                                       alpha=0.5, facecolor=color, lw=0)
            plt.gca().add_patch(patch)
        #plt.fill(*zip(*polygon.exterior.coords), edgecolor='none',
        #         alpha=0.5, facecolor=color)


def wave_data():
    import numpy as np
    # make up some data on a regular lat/lon grid.
    nlats = 73
    nlons = 145
    delta = 2. * np.pi / (nlons - 1)
    lats = (0.5 * np.pi - delta * np.indices((nlats, nlons))[0, :, :])
    lons = (delta * np.indices((nlats, nlons))[1, :, :])
    wave = 0.75 * (np.sin(2. * lats) ** 8 * np.cos(4. * lons))
    mean = 0.5 * np.cos(2. * lats) * ((np.sin(2. * lats)) ** 2 + 2.)
    lats = np.rad2deg(lats)
    lons = np.rad2deg(lons)
    data = wave + mean
    return lons, lats, data


def test(projections):
    coords = [(-0.08, 51.53), (132.00, 43.17)]  # London to Vladivostock
    orig_line_string = sgeom.LineString(coords)

    n_rows = math.ceil(math.sqrt(len(projections)))
    n_cols = math.ceil(len(projections) / n_rows)

    figure, axes_grid = plt.subplots(int(n_rows), int(n_cols))
    if n_rows == 1 and n_cols == 1:
        axes_list = [axes_grid]
    else:
        axes_list = axes_grid.flat

    for projection, axes in zip(projections, axes_list):
        plt.sca(axes)

        colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

        bits = (
            #'contour',
            #'contourf',
            'boundary',
            #'line',
            #'grid',
            'coastline',
            #'polygons',
            'continents',
        )

        if 'contour' in bits:
            # Contours - placeholder for MPL integration
            cs = plt.contour(*wave_data())
            plt.cla()
            for c in cs.collections:
            #for c in cs.collections[2:3]:
                for p in c.get_paths():
                #for p in c.get_paths()[1:]:
                    xy = [segment[0] for segment in p.iter_segments()]
                    line_string = sgeom.LineString(xy)
                    #line_string = sgeom.LineString(xy[:3])
                    draw_line_string(projection, line_string,
                                     color=c.get_color()[0])

        if 'contourf' in bits:
            # Filled contours - placeholder for MPL integration
            cs = plt.contourf(*wave_data())
            plt.cla()
            for c in cs.collections:
            #for i, c in enumerate(cs.collections[2:3]):
                for p in c.get_paths():
                #for j, p in enumerate(c.get_paths()[1:2]):
                    xy = [segment[0] for segment in p.iter_segments()]
                    xy = [xy for xy in xy if xy[1] > -90]
                    polygon = sgeom.Polygon(xy)
                    #polygon = sgeom.Polygon(xy[53:56])
                    draw_polygon(projection, polygon,
                                 color=c.get_facecolor()[0])

        if 'boundary' in bits:
            plt.plot(*list(zip(*projection.boundary.coords)), marker='')

        if 'line' in bits:
            draw_line_string(projection, orig_line_string, color='red')
            polygon = sgeom.LineString([(-50, -80), (90, -80)])
            draw_line_string(projection, polygon)

        if 'grid' in bits:
            # Grid lines
            step = 15
            lons = list(range(0, 360, step))
            for lon in lons:
                line_string = sgeom.LineString(
                    [(lon, -75), (lon, 0), (lon, 75)])
                draw_line_string(projection, line_string, linestyle=':')
            lons = lons + [lons[0]]
            lats = list(range(-90 + step, 90, step))
            for lat in lats:
                line_string = sgeom.LineString([(lon, lat) for lon in lons])
                draw_line_string(projection, line_string, linestyle=':')

        if 'coastline' in bits:
            reader = shapereader.Reader(COASTLINE_PATH)
            print('Reading coastline ...')
            all_geometries = list(reader.geometries())
            print('   ... done.')
            geometries = []
            geometries += all_geometries
            #geometries += all_geometries[48:52] # Aus & Taz
            #geometries += all_geometries[72:73] # GB
            #for geometry in geometries:
            for i, geometry in enumerate(geometries):
                for line_string in geometry:
                    try:
                        draw_line_string(projection, line_string)
                    except ValueError:
                        print(i)
                        print(geometry)
                        raise
                import sys
                sys.stdout.write('.')
                sys.stdout.flush()

        if 'polygons' in bits:
            # Square over pole (CW)
            polygon = sgeom.Polygon(
                [(0, 75), (-90, 75), (-180, 75), (-270, 75)])
            draw_polygon(projection, polygon)

            # Square (CW)
            polygon = sgeom.Polygon(
                [(150, 75), (-150, 75), (-150, 55), (150, 55)])
            draw_polygon(projection, polygon)

            # Wedge - demonstrates removal of interior when split (CW)
            polygon = sgeom.Polygon([(-5, 10), (20, 0), (-5, -10), (10, 0)])
            draw_polygon(projection, polygon)

            # "Antarctica" (incl. non-physical boundary segments) (CW)
            polygon = sgeom.Polygon([(-50, -80), (90, -80), (160, -70),
                                     (160, -90), (-160, -90), (-160, -70)])
            draw_polygon(projection, polygon)

            # Wedge
            polygon = sgeom.Polygon([(-10, 30), (10, 60), (10, 50)])
            draw_polygon(projection, polygon)

        if 'continents' in bits:
            reader = shapereader.Reader(LAND_PATH)
            print('Reading continents ...')
            all_geometries = list(reader.geometries())
            print('   ... done.')
            geometries = []
            geometries += all_geometries

            #geometries += all_geometries[7:8] # Antarctica
            #geometries += all_geometries[16:17] # Some E-equatorial island
            #geometries += all_geometries[93:94] # Some NE island
            #geometries += all_geometries[112:113] # Africa & Asia

            #geometries += all_geometries[95:96] # North and South America
            #geometries += all_geometries[126:] # Greenland

            #geometries += all_geometries[0:7]
            #geometries += all_geometries[8:]
            #geometries += all_geometries[8:16]
            #geometries += all_geometries[17:93]
            #geometries += all_geometries[94:112]
            #geometries += all_geometries[113:]

            for i, multi_polygon in enumerate(geometries):
                for polygon in multi_polygon:
                    polygon = sgeom.Polygon([xy for xy in polygon.exterior.coords if xy[1] > -90])
                    draw_polygon(projection, polygon, color=next(colors))
                    #draw_line_string(projection, polygon)
                import sys
                sys.stdout.write('.')
                sys.stdout.flush()

        plt.title(type(projection).__name__)
        plt.xlim(projection.x_limits)
        plt.ylim(projection.y_limits)


if __name__ == '__main__':
    if 'anim' not in sys.argv:
        projections = [
            ccrs.PlateCarree(-105),
            #ccrs.PlateCarree(20),
            #ccrs.TransverseMercator(central_longitude=-90),
            #ccrs.NorthPolarStereo(),
            #ccrs.Mercator(),
            #ccrs.LambertCylindrical(),
            #ccrs.Robinson(),
            #ccrs.Robinson(170.5),
            #ccrs.Miller(),
            #ccrs.Mollweide(),
            #ccrs.Stereographic(),
            #ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5),
            #ccrs.OSGB(),

            # Incorrect lat/lon grid lines
            #ccrs.Orthographic(central_longitude=-90, central_latitude=45),

            # Not fully implemented (e.g. missing boundary definition)
            #ccrs.InterruptedGoodeHomolosine(),
        ]
        test(projections)
        plt.show()
    else:
        plt.figure(figsize=(15, 6))
        plt.ion()
        while True:
            for lon in range(0, 360, 10):
                print(lon)
                projections = [
                    ccrs.PlateCarree(lon),
                    ccrs.NorthPolarStereo(),
                    ccrs.Robinson(lon),
                ]
                plt.clf()
                plt.suptitle(lon)
                test(projections)
                plt.draw()
                import time
                #time.sleep(1)
