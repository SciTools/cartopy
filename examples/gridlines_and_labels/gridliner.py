"""
Gridlines and tick labels
-------------------------

These examples demonstrate how to quickly add longitude
and latitude gridlines and tick labels on a non-rectangular projection.

As you can see on the first example,
longitude labels may be drawn on left and right sides,
and latitude labels may be drawn on bottom and top sides.
Thanks to the ``dms`` keyword, minutes are used when appropriate
to display fractions of degree.

In the second example, labels are still drawn at the map edges
despite its complexity, and some others are also drawn within the map
boundary.

In the third example, labels are drawn only on the left and bottom sides.
"""
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import matplotlib.pyplot as plt


def main():

    rotated_crs = ccrs.RotatedPole(pole_longitude=120.0, pole_latitude=70.0)
    ax0 = plt.axes(projection=rotated_crs)
    ax0.set_extent([-6, 1, 47.5, 51.5], crs=ccrs.PlateCarree())
    ax0.add_feature(cfeature.LAND.with_scale('110m'))
    ax0.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

    plt.figure(figsize=(6.9228, 3))
    ax1 = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
    ax1.coastlines(resolution='110m')
    ax1.gridlines(draw_labels=True)

    plt.figure(figsize=(7, 3))
    ax2 = plt.axes(projection=ccrs.PlateCarree())
    ax2.coastlines(resolution='110m')
    gl = ax2.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    plt.show()

    plt.figure(figsize=(7, 3))
    ax3 = plt.axes(projection=ccrs.PlateCarree())
    ax3.set_extent([-65, -40, -15, 10])

    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    ax3.add_feature(states_provinces, edgecolor='gray')

    ax3.coastlines(resolution='110m')
    ax3.coastlines(resolution='110m')
    gl = ax3.gridlines(draw_labels=True)

    gl.change_gridline_tick_decimal_separator('{0:.3f}',
                                              axis='both')

    gl.set_latitude_hemisphere_str('Norte', 'Sul')

    gl.set_longitude_hemisphere_str('O', 'L')

    gl.top_labels = False
    gl.right_labels = False
    plt.show()

    return plt.gcf().get_axes(), gl


def gridliner_with_custom_changes_in_its_ticklabels():
    import cartopy.feature as cfeature
    plt.figure(figsize=(7, 3))
    ax3 = plt.axes(projection=ccrs.PlateCarree())
    ax3.set_extent([-65, 40, -15, 10])

    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    ax3.add_feature(states_provinces, edgecolor='gray')

    ax3.coastlines(resolution='110m')
    ax3.coastlines(resolution='110m')
    gl = ax3.gridlines(draw_labels=True)

    gl.change_gridline_tick_decimal_separator('{0:.2f}',
                                              axis='both',
                                              decimal_separator=',')

    gl.set_latitude_hemisphere_str(' - Norte', ' - Sul')

    gl.set_longitude_hemisphere_str('Oeste', 'Leste')

    gl.top_labels = False
    gl.right_labels = False
    plt.show()
    plt.close('all')


if __name__ == '__main__':

    gridliner_with_custom_changes_in_its_ticklabels()

    axes, gl = main()

    print('N° of axes ', len(axes))

    ax = axes[0]
