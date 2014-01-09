import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs


def gen_data():
    # Generate some data
    x = np.arange(-60, 45, 5)
    y = np.arange(30, 80, 5)

    x2d, y2d = np.meshgrid(x, y)
    u = 40 * np.cos(np.deg2rad(y2d))
    v = 40 * np.cos(2. * np.deg2rad(x2d))
    mag = (u**2 + v**2)**.5
    return dict(x=x, y=y, u=u, v=v, magnitude=mag)


def main():
    # Get some data with x, y, u, v, components
    wind = gen_data()

    # Setup our figure
    plot_extent = [-60, 40, 30, 70]

    # Stream plots
    ax = plt.axes(projection=ccrs.RotatedPole(
        pole_latitude=45, pole_longitude=180))
    ax.set_extent(plot_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.stock_img()
    plt.title('Stream plot, North Polar Stereographic projection')
    ax.streamplot(wind['x'], wind['y'], wind['u'], wind['v'], density=1,
                  arrowsize=1, transform=ccrs.PlateCarree(),
                  linewidth=3 * wind['magnitude']/wind['magnitude'].max(),
                  color=wind['magnitude'])

    plt.show()


if __name__ == '__main__':
    main()
