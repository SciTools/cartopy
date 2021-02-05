# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import cartopy.crs as ccrs
import cartopy.geodesic as cgeo
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import (AuxTransformBox, VPacker, HPacker,
                                  TextArea)
from matplotlib.offsetbox import AnchoredOffsetbox
import matplotlib.transforms as transforms


def sbs_to_patch(sbs, transform, unit, padding=2,
                 bbox_transform='axes fraction',
                 bbox_to_anchor=(0.2, 0.3)):

    # First create a single Patch for all Sbs
    of1 = HPacker(width=2,
                  height=1,
                  pad=1,
                  sep=0,
                  align="center",
                  mode="expand", children=sbs)

    t = AnchoredOffsetbox("upper left",
                          pad=0.4, frameon=False,
                          bbox_transform=bbox_transform,
                          bbox_to_anchor=bbox_to_anchor,
                          child=of1)

    return t


def get_unit_converter(unit):

    lookuptable = {'km': 1000,
                   'mi': 1.60934 * 1000,
                   'dm': 1e-1,
                   'cm': 1e-2,
                   'mm': 1e-3,
                   'um': 1e-6,
                   'nm': 1e-9}  # Miles to Km

    return lookuptable.get(unit, 'km')


def _point_along_line(ax, start, distance, projected=False, verbose=False):
    """Point at a given distance from start at a given angle.

    Args:
        ax:       CartoPy axes.
        start:    Starting point for the line in data coordinates.
        distance: Positive physical distance to travel in meters.
        angle:    Anti-clockwise angle for the bar, in degrees. Default: 0

    Returns:
        (lon,lat) coords of a point (a (2, 1)-shaped NumPy array)
    """

    # Direction vector of the line in axes coordinates.

    if not projected:

        geodesic = cgeo.Geodesic()

        Direct_R = geodesic.direct(start, 90, distance)

        target_longitude, target_latitude, forw_azi = Direct_R.base.T

        target_point = ([target_longitude[0], target_latitude[0]])

        actual_dist = geodesic.inverse(start,
                                       target_point).base.ravel()[0]
        if verbose:

            print('Starting point', start)

            print('target point', target_point)
            print('Expected distance between points: ', distance)

            print('Actual distance between points: ', actual_dist)

    if projected:

        longitude, latitude = start

        target_longitude = longitude + distance

        target_point = (target_longitude, latitude)

        if verbose:
            print('Axes is projected? ', projected)
            print('Expected distance between points: ', distance)

            print('Actual distance between points: ',
                  target_longitude - longitude)

    return start, target_point


class AnchoredScaleBar(AnchoredOffsetbox):
    def __init__(self, ax,
                 transform,
                 width,
                 height,
                 zorder,
                 xlabel,
                 fc,
                 ylabels=None,
                 loc=4,
                 fontsize=5,
                 pad=0.1,
                 borderpad=0.1,
                 sep=2,
                 prop=None,
                 add_ruler=False,
                 ruler_unit='Km',
                 ruler_unit_fontsize=7,
                 ruler_fontweight='bold',
                 tick_fontweight='light',
                 **kwargs):
        """
        Draw a horizontal and/or vertical  bar with the size in
        data coordinate of the give axes. A label will be drawn
        underneath (center-aligned).

        - transform : the coordinate frame (typically axes.transData)

        - sizex,sizey : width of x,y bar, in data units. 0 to omit

        - labelx,labely : labels for x,y bars; None to omit

        - loc : position in containing axes

        - pad, borderpad : padding, in fraction of the legend
        font size (or prop)

        - sep : separation between labels and bars in points.

        - **kwargs : additional arguments passed to base class

        constructor
        """

        if ruler_unit_fontsize is None:
            ruler_unit_fontsize = fontsize * 1.5

        Rect = Rectangle((0, 0),
                         width, height, fc=fc,
                         edgecolor='k',
                         zorder=zorder,)

        ATB = AuxTransformBox(transform)

        ATB.add_artist(Rect)

        Txt_xlabel = TextArea(xlabel,
                              textprops=dict(fontsize=fontsize,
                                             fontweight=tick_fontweight),
                              minimumdescent=True)

        # vertically packing a single stripe with respective label

        child = VPacker(children=[Txt_xlabel,
                                  ATB],
                        align="right", pad=5, sep=0)

        if add_ruler:

            Text = TextArea(ruler_unit,
                            textprops=dict(fontsize=ruler_unit_fontsize,
                                           fontweight=ruler_fontweight))

            child = VPacker(children=[child, Text],
                            align="center", pad=5, sep=0)

        else:

            Text = TextArea('',
                            textprops=dict(fontsize=ruler_unit_fontsize))

            child = VPacker(children=[child, Text],
                            align="right", pad=5, sep=0)

        # horizontally packing all child packs in a single offsetBox

        AnchoredOffsetbox.__init__(self,
                                   loc='center left',
                                   borderpad=borderpad,
                                   child=child,
                                   prop=prop,
                                   frameon=False,
                                   **kwargs)


def _add_scalebar(ax,
                  projected,
                  xcoords,
                  height,
                  xlabels=None,
                  ylabels=None,
                  loc=4,
                  bbox_to_anchor=(0.2, 0.5),
                  bbox_transform='axes fraction',
                  frameon=False,
                  fontsize=5,
                  tick_fontweight='light',
                  ruler_unit_fontsize=None,
                  ruler_unit='Km',
                  ruler_fontweight='bold',
                  **kwargs):
    """ Add scalebars to axes
    Adds a set of scale bars to *ax*, matching the size
    to the ticks of the plot
    and optionally hiding the x and y axes
    - ax : the axis to attach ticks to
    - matchx,matchy : if True, set size of scale bars to spacing
    between ticks
                    if False, size should be set using sizex and
                    sizey params

    - hidex,hidey : if True, hide x-axis and y-axis of parent

    - **kwargs : additional arguments passed to AnchoredScaleBars

    Returns
        created scalebar object
    """

    if not projected:
        proj = ax.projection._as_mpl_transform(ax)
    else:
        proj = ax.transData
    blended_transform = transforms.blended_transform_factory(
        proj, ax.get_figure().dpi_scale_trans)

    SBs = []

    average_nticks = min(len(xlabels) - 1, len(xcoords) - 1)

    for enum, (xcor, xlabel) in enumerate(zip(xcoords[1:],
                                              xlabels)):
        width = xcor - xcoords[0]
        if enum % 2 == 0:
            fc = 'white'
        else:
            fc = 'black'

        xlabel = int(xlabel)

        zorder = 999 - enum

        if enum == average_nticks:
            add_ruler = True

        else:  # odd number
            add_ruler = False

        sb = AnchoredScaleBar(ax,
                              blended_transform,
                              width,
                              height,
                              zorder,
                              xlabel,
                              fc=fc,
                              ylabels=ylabels,
                              loc=loc,
                              fontsize=fontsize,
                              bbox_transform=ax.transAxes,
                              bbox_to_anchor=bbox_to_anchor,
                              add_ruler=add_ruler,
                              ruler_unit=ruler_unit,
                              ruler_unit_fontsize=ruler_unit_fontsize,
                              ruler_fontweight=ruler_fontweight,
                              tick_fontweight=tick_fontweight,
                              **kwargs)

        sb.set_clip_on(False)
        sb.set_zorder(zorder)
        SBs.append(sb)

    # SB = sbs_to_patch(SBs,
    #                   blended_transform,
    #                   ruler_unit, padding=2,
    #                   bbox_transform=ax.transAxes,
    #                   bbox_to_anchor=bbox_to_anchor,)
    # ax.add_artist(SB)
    for sb in SBs:
        ax.add_artist(sb)

    return sb


def add_scalebar(ax,
                 bbox_to_anchor,
                 length,
                 ruler_unit='km',
                 ruler_fontweight='bold',
                 tick_fontweight='light',
                 ruler_unit_fontsize=10,
                 dy=5,
                 max_stripes=5,
                 ytick_label_margins=0.25,
                 fontsize=8,

                 frameon=False,
                 font_weight='bold',
                 rotation=45,
                 zorder=999,
                 paddings={'xmin': 0.1,
                             'xmax': 0.1,
                             'ymin': 0.3,
                             'ymax': 0.8},

                 bbox_kwargs={'facecolor': 'w',
                              'edgecolor': 'k',
                              'alpha': 0.7},
                 numeric_scale_bar=True,
                 numeric_scale_bar_kwgs={'x_text_offset': 0,
                                         'y_text_offset': -40,
                                         'box_x_coord': 0.5,
                                         'box_y_coord': 0.01},
                 verbose=False):
    '''
    Description

    ----------
        This function draws a scalebar in the given geoaxes.

    Parameters
    ----------
        ax (geoaxes):

        location (length 2 tuple):
            It sets where the scalebar will be drawn
            in axes fraction units.

        length (float):
            The distance in geodesic meters that will be used
            for generating the scalebar.

        unit_name (str):
            Standard (km).


        angle (int or float): in azimuth degrees.
            The angle that will be used for evaluating the scalebar.

            If 90 (degrees), the distance between each tick in the
            scalebar will be evaluated in respect to the longitude
            of the map.

            If 0 (degrees), the ticks will be evaluated in accordance
            to variation in the latitude of the map.

        dy (int or float):
            The hight of the scalebar in axes fraction.

        max_stripes (int):
            The number of stripes present in the scalebar.

        ytick_label_margins (int or float):
            The size of the margins for drawing the scalebar ticklabels.

        fontsize (int or float):
            The fontsize used for drawing the scalebar ticklabels.

        font_weight (str):
            the fontweight used for drawing the scalebar ticklabels.

        rotation (int or float):
            the rotation used for drawing the scalebar ticklabels.

        zorder(int):
            The zorder used for drawing the scalebar.

        paddings (dict):
            A dictionary defining the padding to draw a background box
            around the scalebar.

            Example of allowed arguments for padding:
                {'xmin': 0.3,
                 'xmax': 0.3,
                 'ymin': 0.3,
                 'ymax': 0.3}

        bbox_kwargs (dict):
            A dictionary defining the background box
            around the scalebar.

            Example of allowed arguments for padding:
                {'facecolor': 'w',
                 'edgecolor': 'k',
                 'alpha': 0.7}

        numeric_scale_bar(bool):
            whether or not to draw a number scalebar along side the
            graphic scalebar. Notice that this option can drastically
            vary in value, depending on the geoaxes projection used.

        numeric_scale_bar_kwgs (dict):
            A dictionary defining the numeric scale bar.

            Example of allowed arguments:
                {'x_text_offset': 0,
                 'y_text_offset': -40,
                 'box_x_coord': 0.5,
                 'box_y_coord': 0.01}

    Returns
    ----------
    None
    '''

    proj_units = ax.projection.proj4_params.get('units', 'degrees')
    if proj_units.startswith('deg'):
        projected = False

    elif proj_units.startswith('m'):
        projected = True

    # getting the basic unit converter for labeling the xticks
    unit_converter = get_unit_converter(ruler_unit)

    if verbose:
        print('Axes is projected? ', projected)

    # Map central XY data coordinates
    x0, x1, y0, y1 = ax.get_extent()

    central_coord_map = np.mean([[x0, x1], [y0, y1]], axis=1).tolist()

    # End-point of bar in lon/lat coords.
    start, end = _point_along_line(ax,
                                   central_coord_map,
                                   length,
                                   projected=projected,
                                   verbose=verbose)

    # choose exact X points as sensible grid ticks with Axis 'ticker' helper
    xcoords = np.empty(max_stripes + 1)
    xlabels = []

    xcoords[0] = start[0]

    ycoords = np.empty_like(xcoords)

    for i in range(0, max_stripes):

        startp, endp = _point_along_line(ax, central_coord_map,
                                         length * (i + 1),
                                         projected=projected)
        # to ensure that the scalebar will not be so long as to
        # cause errors in the plot.
        if i > 0:
            if endp[0] < xcoords[i]:
                break

        xcoords[i + 1] = endp[0]

        ycoords[i + 1] = end[1]

        label = round(length * (i + 1) / unit_converter)

        xlabels.append(label)

    # Stacking data coordinates (the target ticks of the scalebar) in a list

    _add_scalebar(ax,
                  projected,
                  xcoords,
                  dy,
                  xlabels=xlabels,
                  ylabels=None,
                  loc=4,

                  frameon=frameon,
                  bbox_to_anchor=bbox_to_anchor,
                  fontsize=fontsize,
                  tick_fontweight=tick_fontweight,
                  ruler_unit_fontsize=ruler_unit_fontsize,
                  ruler_unit=ruler_unit,
                  ruler_fontweight=ruler_fontweight,
                  )


if '__main__' == __name__:

    def test_scalebar():
        """Test"""

        fig, axes = plt.subplots(1, 2,
                                 subplot_kw={'projection':
                                             ccrs.Mercator()})

        projections = [ccrs.Mercator(), ccrs.PlateCarree()]

        axes = axes.ravel()

        for proj, ax in zip(projections, axes):

            ax.projection = proj

            ax.set_title(ax.projection.__class__.__name__)

            add_scalebar(ax,
                         bbox_to_anchor=(0.1, 0.2),
                         length=10_000_000,
                         ruler_unit='km',
                         max_stripes=3,
                         fontsize=8,
                         frameon=True,
                         ruler_unit_fontsize=15,
                         ruler_fontweight='bold',
                         tick_fontweight='bold',
                         dy=0.085)

            ax.gridlines(draw_labels=True)
            ax.stock_img()
            ax.coastlines()

        return axes

    axes = test_scalebar()
