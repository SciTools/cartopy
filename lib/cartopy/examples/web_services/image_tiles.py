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

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=mercator)
    ax.set_extent([-90, -73, 22, 34])

    ax.add_image(tiler, 6)

    ax.coastlines('10m')
    plt.show()


if __name__ == '__main__':
    main()
