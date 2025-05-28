"""
Tube Stations
-------------

Produces a map showing London Underground station locations with high
resolution background imagery provided by OpenStreetMap.

"""

from matplotlib.path import Path
import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
from cartopy.io.img_tiles import OSM


def tube_locations():
    """
    Return an (n, 2) array of selected London Tube locations in Ordnance
    Survey GB coordinates.

    Source: https://www.doogal.co.uk/london_stations.php

    """
    return np.array(
        [
            [531738.0, 180890.0],
            [532379.0, 179734.0],
            [531096.0, 181642.0],
            [530234.0, 180492.0],
            [531688.0, 181150.0],
            [530242.0, 180982.0],
            [531940.0, 179144.0],
            [530406.0, 180380.0],
            [529012.0, 180283.0],
            [530553.0, 181488.0],
            [531165.0, 179489.0],
            [529987.0, 180812.0],
            [532347.0, 180962.0],
            [529102.0, 181227.0],
            [529612.0, 180625.0],
            [531566.0, 180025.0],
            [529629.0, 179503.0],
            [532105.0, 181261.0],
            [530995.0, 180810.0],
            [529774.0, 181354.0],
            [528941.0, 179131.0],
            [531050.0, 179933.0],
            [530240.0, 179718.0],
        ]
    )


def main():
    imagery = OSM()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=imagery.crs)
    ax.set_extent([-0.14, -0.1, 51.495, 51.515], ccrs.PlateCarree())

    # Construct concentric circles and a rectangle,
    # suitable for a London Underground logo.
    theta = np.linspace(0, 2 * np.pi, 100)
    circle_verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    concentric_circle = Path.make_compound_path(
        Path(circle_verts[::-1]), Path(circle_verts * 0.6)
    )

    rectangle = Path([[-1.1, -0.2], [1, -0.2], [1, 0.3], [-1.1, 0.3]])

    # Add the imagery to the map.
    ax.add_image(imagery, 14)

    # Plot the locations twice, first with the red concentric circles,
    # then with the blue rectangle.
    xs, ys = tube_locations().T
    ax.plot(
        xs,
        ys,
        transform=ccrs.OSGB(approx=False),
        marker=concentric_circle,
        color='red',
        markersize=9,
        linestyle='',
    )
    ax.plot(
        xs,
        ys,
        transform=ccrs.OSGB(approx=False),
        marker=rectangle,
        color='blue',
        markersize=11,
        linestyle='',
    )

    ax.set_title('London underground locations')
    plt.show()


if __name__ == '__main__':
    main()
