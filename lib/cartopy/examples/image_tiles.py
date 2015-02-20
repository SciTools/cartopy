__tags__ = ['Web services']
"""
Web tile imagery
----------------

This example demonstrates how imagery from a tile
providing web service can be accessed.

"""
import matplotlib.pyplot as plt

from cartopy.io.img_tiles import StamenTerrain


def main():
    tiler = StamenTerrain()
    mercator = tiler.crs
    ax = plt.axes(projection=mercator)
    ax.set_extent([-90, -73, 22, 34])

    ax.add_image(tiler, 6)

    ax.coastlines('10m')
    plt.show()


if __name__ == '__main__':
    main()
