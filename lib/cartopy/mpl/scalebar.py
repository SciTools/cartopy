import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import font_manager as mfonts
import matplotlib.patches as patches


def _axes_to_lonlat(ax, coords):
    """description:
        Transform the axes coordinates into (lon, lat)

       Returns tuple(n,2): in the (lon,lat) format"""
    display = ax.transAxes.transform(coords)
    data = ax.transData.inverted().transform(display)
    lonlat = ccrs.PlateCarree().transform_point(*data, ax.projection)

    return lonlat


def geodesy_distance_between_points(p_start, p_end):
    geodesic = cgeo.Geodesic()

    distances, start_azimuth, end_azimuth = geodesic.inverse(
        p_start, p_end).base.T

    return distances


def _point_along_line(ax, start, distance, angle=-90):
    """Point at a given distance from start at a given angle.

    Args:
        ax:       CartoPy axes.
        start:    Starting point for the line in axes coordinates.
        distance: Positive physical distance to travel in meters.
        angle:    Anti-clockwise angle for the bar, in degrees. Default: 0

    Returns:
        (lon,lat) coords of a point (a (2, 1)-shaped NumPy array)
    """
    # Direction vector of the line in axes coordinates.
    geodesic = cgeo.Geodesic()

    start_coords = _axes_to_lonlat(ax, start)

    Direct_R = geodesic.direct(start_coords, angle, distance)

    longitudes, latitudes, forw_azi = Direct_R.base.T

    # print('Distance', distance)
    # print('start_coords: ', start_coords)
    # print('longitudes: ', longitudes)

    target_point = (longitudes, latitudes)

    # actual_dist = geodesic.inverse(start_coords, 
    #                                target_point).base.ravel()[0]
    # print('Starting point', start_coords)

    # print('Ending point', target_point)

    # print('Expected distance between points: ', distance)

    # print('Actual distance between points: ', actual_dist)

    return start_coords, target_point


def _add_bbox(ax, list_of_patches, paddings={}, bbox_kwargs={}):
    '''
    Description:
        This helper function adds a box behind the scalebar:
            Code inspired by:
    https://stackoverflow.com/questions/17086847/box-around-text-in-matplotlib


    '''

    zorder = list_of_patches[0].get_zorder() - 1

    xmin = min([t.get_window_extent().xmin for t in list_of_patches])
    xmax = max([t.get_window_extent().xmax for t in list_of_patches])
    ymin = min([t.get_window_extent().ymin for t in list_of_patches])
    ymax = max([t.get_window_extent().ymax for t in list_of_patches])

    xmin, ymin = ax.transAxes.inverted().transform((xmin, ymin))
    xmax, ymax = ax.transAxes.inverted().transform((xmax, ymax))

    xmin = xmin - ((xmax - xmin) * paddings['xmin'])
    ymin = ymin - ((ymax - ymin) * paddings['ymin'])

    xmax = xmax + ((xmax - xmin) * paddings['xmax'])
    ymax = ymax + ((ymax - ymin) * paddings['ymax'])

    width = (xmax - xmin)
    height = (ymax - ymin)

    # Setting xmin according to height

    rect = patches.Rectangle((xmin, ymin),
                             width,
                             height,
                             facecolor=bbox_kwargs['facecolor'],
                             edgecolor=bbox_kwargs['edgecolor'],
                             alpha=bbox_kwargs['alpha'],
                             transform=ax.transAxes,
                             fill=True,
                             clip_on=False,
                             zorder=zorder)

    ax.add_patch(rect)
    return ax, rect


def fancy_scalebar(ax,
                   location,
                   length,
                   unit_name='km',
                   tol=0.01,
                   angle=90,
                   dy=5,

                   max_stripes=5,
                   ytick_label_margins=0.25,
                   fontsize=8,
                   font_weight='bold',
                   rotation=45,
                   zorder=999,
                   paddings={'xmin': 0.3,
                               'xmax': 0.3,
                               'ymin': 0.3,
                               'ymax': 0.3},

                   bbox_kwargs={'facecolor': 'w',
                                'edgecolor': 'k',
                                'alpha': 0.7},
                   numeric_scale_bar=True,
                   numeric_scale_bar_kwgs={'x_text_offset': 0,
                                           'y_text_offset': -40,
                                           'box_x_coord': 0.5,
                                           'box_y_coord': 0.01}
                   ):

    # Convert all units and types.
    location = np.asarray(location)  # For vector addition.

    # End-point of bar in lon/lat coords.
    start, end = _point_along_line(ax, location, length, angle=angle)

    # choose exact X points as sensible grid ticks with Axis 'ticker' helper
    xcoords = np.empty(max_stripes + 1)
    xlabels = [0]

    xcoords[0] = start[0]

    for i in range(0, max_stripes):

        startp, endp = _point_along_line(
            ax, location, length * (i + 1), angle=angle)

        xcoords[i + 1] = endp[0]

        label = length * (i + 1)

        xlabels.append(label)

    print('xcoords: ', xcoords)
    print('xlabels', xlabels)

    # grab min+max for limits
    xl0, xl1 = xcoords.min(), xcoords.max()

    print('Min - Max coords: ', xl0, xl1)

    # calculate Axes Y coordinates of box top+bottom

    yl0 = start[1]

    yl1 = yl0 + dy

    y_margin = (yl1 - yl0) * ytick_label_margins

    # Setting offset transformer
    transform = ax.transData

    # fill black/white 'stripes' and draw their boundaries

    fill_colors = ['black', 'white']
    i_color = 0

    filled_boxs = []
    for xi0, xi1 in zip(xcoords[:-1], xcoords[1:]):

        # fill region
        filled_box = plt.fill(
            (xi0, xi1, xi1, xi0, xi0),
            (yl0, yl0, yl1, yl1, yl0),

            fill_colors[i_color],
            transform=transform,
            clip_on=False,
            zorder=zorder
        )

        filled_boxs.append(filled_box[0])

        # draw boundary
        plt.plot((xi0, xi1, xi1, xi0, xi0),
                 (yl0, yl0, yl1, yl1, yl0),
                 'black',
                 clip_on=False,
                 transform=transform,
                 zorder=zorder)

        i_color = 1 - i_color

    # adding boxes

    ax, rect = _add_bbox(ax,
                         filled_boxs,
                         bbox_kwargs=bbox_kwargs,
                         paddings=paddings)

    # add short tick lines
    for x in xcoords:
        plt.plot((x, x),
                 (yl0, yl0 - y_margin),
                 'black',
                 transform=transform,
                 zorder=zorder,
                 clip_on=False)

    # add a scale legend 'Km'
    font_props = mfonts.FontProperties(size=fontsize,
                                       weight=font_weight)

    plt.text(0.5 * (xl0 + xl1),
             yl1 + y_margin,
             unit_name,
             color='k',
             verticalalignment='bottom',
             horizontalalignment='center',
             fontproperties=font_props,
             transform=transform,
             clip_on=False,
             zorder=zorder)

    # add numeric labels

    if unit_name == 'km':
        divider = 1e-3

    for x, xlabel in zip(xcoords, xlabels):
        plt.text(x,
                 yl0 - 2 * y_margin,
                 '{:g}'.format((xlabel * divider)),
                 verticalalignment='top',
                 horizontalalignment='center',
                 fontproperties=font_props,
                 transform=transform,
                 rotation=rotation,
                 clip_on=False,
                 zorder=zorder + 1,
                 # bbox=dict(facecolor='red', alpha=0.5) # this would add a box
                 # only around the xticks
                 )

    # Adjusting figure borders to ensure that the scalebar is within its limits

    # ax.get_figure().canvas.draw()
    # ax.get_figure().tight_layout()

    # get rectangle background bbox

    if numeric_scale_bar:

        add_numeric_scale_bar(ax,
                              rect,
                              numeric_scale_bar_kwgs,
                              fontprops=font_props)


def _add_numeric_scale_bar(ax, inches_to_cm=1 / 2.54):
    '''
    Description:
        This function adds a text object, which contains
        the respective numeric scale of the map.


    Parameters:
        ax (cartopy geoaxes)

        inches_to_cm (float): the standard ration of inches to cm
                              Standard value: inches_to_cm=1/2.54
    Returns
        dx_fig(float): the figure relative

        dx_mapa/10 (float): geographical distance of 1cm
                            in the map in respect to ground
    '''

    fig = ax.get_figure()

    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

    bbox_in_data_coords = (ax.get_window_extent()
                           .transformed(ax.transData.inverted())
                           )

    dx_fig = bbox.width * inches_to_cm  # width in cms

    # Getting distance:
    x0, x1, y0, y1 = ax.get_extent()

    proj4_params = ax.projection.proj4_params

    units = proj4_params.get('units', None)
    # if ax projection is a projected crs:
    if units is not None:
        dx_mapa = x1 - x0

    # in case it is not a projected crs (i.e.: PlateCarree):
    else:

        lon_min = x0
        lat_mean = np.mean([y0, y1])

        # Define starting point.
        start = (lon_min, lat_mean)

        delta_x = bbox_in_data_coords.width  # in degrees

        end = (lon_min + delta_x, lat_mean)

        dx_mapa = geodesy_distance_between_points(start, end)

        # meters to cm
        dx_mapa = dx_mapa * 1e2

    # updating dx_mapa, so that it will always be [1 in fig cm: dx_mapa cm]
    dx_mapa = dx_mapa / dx_fig

    # dividing by 10... It fix the error found by comparing with Qgis (why?)

    return dx_fig, dx_mapa / 10


def add_numeric_scale_bar(ax, patch, numeric_scale_bar_kwgs, fontprops=None):
    '''
    Description:
        This function adds a text object surrounded by a patch.

        The Text objection contains the respective numeric scale of the map.


    Parameters:
        ax (cartopy geoaxes)

        patch (matplotlib.patch.box): the patch that will be used as background
        of the numeric scalebar

    Returns
        None
    '''

    if fontprops is None:
        fontprops = mfonts.FontProperties(size=8,
                                          weight='bold')

    dx_fig, dx_mapa = _add_numeric_scale_bar(ax)

    rx, ry = patch.get_xy()

    # cy = ry + patch.get_height() / 2.0

    xytext = (numeric_scale_bar_kwgs['x_text_offset'],
              numeric_scale_bar_kwgs['y_text_offset']
              )

    xy = (numeric_scale_bar_kwgs['box_x_coord'],
          numeric_scale_bar_kwgs['box_y_coord']
          )

    ax.annotate('1:{0:.0f}'.format(dx_mapa[0]),
                xy=xy,

                xytext=xytext,
                color='black',
                weight='bold',
                zorder=patch.zorder + 1,
                xycoords=patch,
                textcoords='offset points',
                font_properties=fontprops,
                ha='center',
                va='center')
