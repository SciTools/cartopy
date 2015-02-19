__tags__ = ['Scalar data']
"""
This example illustrates the automatic download of
STRM data, gap filling (using gdal) and adding shading
to create a so-called "Shaded Relief SRTM".

Originally contributed by Thomas Lecocq (http://geophysique.be).

"""
import cartopy.crs as ccrs
from cartopy.io import srtm
import matplotlib.pyplot as plt

from cartopy.io import PostprocessedRasterSource, LocatedImage
from cartopy.io.srtm import SRTM3Source


def fill_and_shade(located_elevations):
    """
    Given an array of elevations in a LocatedImage, fill any holes in
    the data and add a relief (shadows) to give a realistic 3d appearance.

    """
    new_elevations = srtm.fill_gaps(located_elevations.image, max_distance=15)
    new_img = srtm.add_shading(new_elevations, azimuth=135, altitude=15)
    return LocatedImage(new_img, located_elevations.extent)


def main():
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Define a raster source which uses the SRTM3 data and applies the
    # fill_and_shade function when the data is retrieved.
    shaded_srtm = PostprocessedRasterSource(SRTM3Source(), fill_and_shade)

    # Add the shaded SRTM source to our map with a grayscale colormap.
    ax.add_raster(shaded_srtm, cmap='Greys')

    # This data is high resolution, so pick a small area which has some
    # interesting orography.
    ax.set_extent([12, 13, 47, 48])

    plt.title("SRTM Shaded Relief Map")

    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_left = False

    plt.show()


if __name__ == '__main__':
    main()
