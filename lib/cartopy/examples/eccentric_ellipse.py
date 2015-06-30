__tags__ = ['Miscellanea']
"""
Eccentric Ellipsoids
--------------------

This example demonstrates how :class:`cartopy.crs.Globe` can be used
to define a highly eccentric elliptical model of a geoid. This globe
definition is used in defining a Geostationary projection. The projection
is then visualised with a Natural Earth image and coastlines, which have both
been reprojected on the fly.

"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def main():
    # We define the semimajor and semiminor axes, but must also tell the
    # globe not to use the WGS84 ellipse, which is its default behaviour.
    eccentric_globe = ccrs.Globe(semimajor_axis=1000, semiminor_axis=500,
                                 ellipse=None)
    geostationary = ccrs.Geostationary(globe=eccentric_globe)

    ax = plt.axes(projection=geostationary)
    ax.stock_img()
    ax.coastlines()
    plt.show()


if __name__ == '__main__':
    main()
