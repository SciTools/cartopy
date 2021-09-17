"""
Contour transform options
=========================

This example demonstrates the difference between transforming
the points before/after generating the contours. It uses the
**transform_first** keyword argument to indicate that Cartopy should
transform the points before calling the contouring algorithm,
which can have a significant impact on speed (it is much faster
to transform points than it is to transform patches). This does
have a negative impact on the wrapped coordinates as one can see in the
second axes that the data does not extend to the full global extent.
"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from waves import sample_data


def main():

    # Use the waves example to provide some sample data, but make it
    # more dependent on y for more interesting contours.
    x, y, z = sample_data((20, 40))
    z = z * -1.5 * y

    # Setup a global EckertIII map with faint coastlines.
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1, projection=ccrs.EckertIII())
    ax1.set_title("transform_first=False")
    ax2 = fig.add_subplot(2, 1, 2, projection=ccrs.EckertIII())
    ax2.set_title("transform_first=True")

    for ax, transform_first in zip([ax1, ax2], [False, True]):
        ax.set_global()
        ax.coastlines('110m', alpha=0.1)

        # Add colourful filled contours.
        filled_c = ax.contourf(x, y, z, transform=ccrs.PlateCarree(),
                               transform_first=transform_first)

        # And black line contours.
        ax.contour(x, y, z, levels=filled_c.levels,
                   colors=['black'],
                   transform=ccrs.PlateCarree(),
                   transform_first=transform_first)

    plt.show()


if __name__ == '__main__':
    main()
