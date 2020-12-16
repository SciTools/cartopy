import cartopy.crs as ccrs
from cartopy import feature as cfeature
import cartopy.geodesic as cgeo
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import font_manager as mfonts
import matplotlib.ticker as mticker
import matplotlib.patches as patches

from geopy.distance import distance, lonlat
import geopy

def _axes_to_lonlat(ax, coords):
    """(lon, lat) from axes coordinates."""
    display = ax.transAxes.transform(coords)
    data = ax.transData.inverted().transform(display)
    lonlat = ccrs.PlateCarree().transform_point(*data, ax.projection)

    return lonlat


def _upper_bound(start, direction, distance, dist_func):
    """A point farther than distance from start, in the given direction.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        direction  Nonzero (2, 1)-shaped array, a direction vector.
        distance:  Positive distance to go past.
        dist_func: A two-argument function which returns distance.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    if distance <= 0:
        raise ValueError(f"Minimum distance is not positive: {distance}")

    if np.linalg.norm(direction) == 0:
        raise ValueError("Direction vector must not be zero.")

    # Exponential search until the distance between start and end is
    # greater than the given limit.
    length = 0.1
    end = start + length * direction

    while dist_func(start, end) < distance:
        length *= 2
        end = start + length * direction

    return end


def _distance_along_line(start, end, distance, dist_func, tol):
    """Point at a distance from start on the segment  from start to end.

    It doesn't matter which coordinate system start is given in, as long
    as dist_func takes points in that coordinate system.

    Args:
        start:     Starting point for the line.
        end:       Outer bound on point's location.
        distance:  Positive distance to travel.
        dist_func: Two-argument function which returns distance.
        tol:       Relative error in distance to allow.

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    initial_distance = dist_func(start, end)
    if initial_distance < distance:
        raise ValueError(f"End is closer to start ({initial_distance}) than "
                         f"given distance ({distance}).")

    if tol <= 0:
        raise ValueError(f"Tolerance is not positive: {tol}")

    # Binary search for a point at the given distance.
    left = start
    right = end

    while not np.isclose(dist_func(start, right), distance, rtol=tol):
        midpoint = (left + right) / 2

        # If midpoint is too close, search in second half.
        if dist_func(start, midpoint) < distance:
            left = midpoint
        # Otherwise the midpoint is too far, so search in first half.
        else:
            right = midpoint

    return right


def _point_along_line(ax, start, distance, angle=0, tol=0.01):
    """Point at a given distance from start at a given angle.

    Args:
        ax:       CartoPy axes.
        start:    Starting point for the line in axes coordinates.
        distance: Positive physical distance to travel.
        angle:    Anti-clockwise angle for the bar, in radians. Default: 0
        tol:      Relative error in distance to allow. Default: 0.01

    Returns:
        Coordinates of a point (a (2, 1)-shaped NumPy array).
    """
    # Direction vector of the line in axes coordinates.
    direction = np.array([np.cos(angle), np.sin(angle)])

    geodesic = cgeo.Geodesic()

    # Physical distance between points.
    def dist_func(a_axes, b_axes):
        a_phys = _axes_to_lonlat(ax, a_axes)
        b_phys = _axes_to_lonlat(ax, b_axes)

        # Geodesic().inverse returns a NumPy MemoryView like [[distance,
        # start azimuth, end azimuth]].
        return geodesic.inverse(a_phys, b_phys).base[0, 0]

    end = _upper_bound(start, direction, distance, dist_func)

    return _distance_along_line(start, end, distance, dist_func, tol)


def scale_bar(ax, location, length, metres_per_unit=1000, unit_name='km',
              tol=0.01, angle=0, color='black', linewidth=3, text_offset=0.005,
              ha='center', va='bottom', plot_kwargs=None, text_kwargs=None,
              **kwargs):
    """Add a scale bar to CartoPy axes.

    For angles between 0 and 90 the text and line may be plotted at
    slightly different angles for unknown reasons. To work around this,
    override the 'rotation' keyword argument with text_kwargs.

    Args:
        ax:              CartoPy axes.
        location:        Position of left-side of bar in axes coordinates.
        length:          Geodesic length of the scale bar.
        metres_per_unit: Number of metres in the given unit. Default: 1000
        unit_name:       Name of the given unit. Default: 'km'
        tol:             Allowed relative error in length of bar. Default: 0.01
        angle:           Anti-clockwise rotation of the bar.
        color:           Color of the bar and text. Default: 'black'
        linewidth:       Same argument as for plot.
        text_offset:     Perpendicular offset for text in axes coordinates.
                         Default: 0.005
        ha:              Horizontal alignment. Default: 'center'
        va:              Vertical alignment. Default: 'bottom'
        **plot_kwargs:   Keyword arguments for plot, overridden by **kwargs.
        **text_kwargs:   Keyword arguments for text, overridden by **kwargs.
        **kwargs:        Keyword arguments for both plot and text.
    """
    # Setup kwargs, update plot_kwargs and text_kwargs.
    if plot_kwargs is None:
        plot_kwargs = {}
    if text_kwargs is None:
        text_kwargs = {}

    plot_kwargs = {'linewidth': linewidth, 'color': color, **plot_kwargs,
                   **kwargs}
    text_kwargs = {'ha': ha, 'va': va, 'rotation': angle, 'color': color,
                   **text_kwargs, **kwargs}

    # Convert all units and types.
    location = np.asarray(location)  # For vector addition.
    length_metres = length * metres_per_unit
    angle_rad = angle * np.pi / 180

    # End-point of bar.
    end = _point_along_line(ax, location, length_metres, angle=angle_rad,
                            tol=tol)

    # Coordinates are currently in axes coordinates, so use transAxes to
    # put into data coordinates. *zip(a, b) produces a list of x-coords,
    # then a list of y-coords.
    ax.plot(*zip(location, end), transform=ax.transAxes, **plot_kwargs)
    
    # Push text away from bar in the perpendicular direction.
    midpoint = (location + end) / 2
    offset = text_offset * np.array([-np.sin(angle_rad), np.cos(angle_rad)])
    text_location = midpoint + offset

    # 'rotation' keyword argument is in text_kwargs.
    ax.text(*text_location, f"{length} {unit_name}", rotation_mode='anchor',
            transform=ax.transAxes, **text_kwargs)




def _add_bbox(ax, list_of_patches, paddings={}, bbox_kwargs={}):
    
    '''
    Description:
        This helper function adds a box behind the scalebar:
            Code inspired by: https://stackoverflow.com/questions/17086847/box-around-text-in-matplotlib
    
    
    '''
    
    zorder = list_of_patches[0].get_zorder() - 1
    
    xmin = min([t.get_window_extent().xmin for t in list_of_patches])
    xmax = max([t.get_window_extent().xmax for t in list_of_patches])
    ymin = min([t.get_window_extent().ymin for t in list_of_patches])
    ymax = max([t.get_window_extent().ymax for t in list_of_patches])
    

    xmin, ymin = ax.transAxes.inverted().transform((xmin, ymin))
    xmax, ymax = ax.transAxes.inverted().transform((xmax, ymax))

    
    xmin = xmin - ( (xmax-xmin) * paddings['xmin'])
    ymin = ymin - ( (ymax-ymin) * paddings['ymin'])
    
    xmax = xmax + ( (xmax-xmin) * paddings['xmax'])
    ymax = ymax + ( (ymax-ymin) * paddings['ymax'])
    
    width = (xmax-xmin)
    height = (ymax-ymin)
    
    # Setting xmin according to height
    
    
    rect = patches.Rectangle((xmin,ymin),
                              width,
                              height, 
                              facecolor=bbox_kwargs['facecolor'], 
                              edgecolor =bbox_kwargs['edgecolor'],
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
                   
                   metres_per_unit=1000, 
                   unit_name='km',
                   tol=0.01, 
                   angle=0,
                   dy = 0.05,
                   
                   max_stripes=5,
                   ytick_label_margins = 0.25,
                   fontsize= 8,
                   font_weight='bold',
                   rotation = 45,
                   zorder=999,
                   paddings = {'xmin':0.3,
                             'xmax':0.3,
                             'ymin':0.3,
                             'ymax':0.3},
    
                 bbox_kwargs = {'facecolor':'w',
                                'edgecolor':'k',
                                'alpha':0.7},
                 add_numeric_scale_bar=True,
                 numeric_scale_bar_kwgs={'x_text_offset':0,
                                         'y_text_offset':-20,
                                         'box_x_coord':0.5,
                                         'box_y_coord':0.01}
                 ):
    
    
   
    
    
    
    # Convert all units and types.
    location = np.asarray(location)  # For vector addition.
    length_metres = length * metres_per_unit
    angle_rad = angle * np.pi / 180

    # End-point of bar.
    end = _point_along_line(ax, location, length_metres, angle=angle_rad,
                            tol=tol)
    
    
    x0 = location[0]
    x1 = end[0]
    ycoord = location[1]
    
    dx = x1 - x0
    
    
    # choose exact X points as sensible grid ticks with Axis 'ticker' helper
    xcoords = []
    ycoords = []
    xlabels = []
    
    for i in range(0 , 1+ max_stripes):
        dlength = (dx * i) + x0
        xlabels.append(  (length_metres * i) )
        
        xcoords.append(dlength)
        ycoords.append(ycoord)
        
    # Convertin x_vals to axes fraction data:
    xcoords = np.asanyarray(xcoords)
    ycoords = np.asanyarray(ycoords)
    
    
    # grab min+max for limits
    xl0, xl1 = xcoords[0], xcoords[-1]
    
    
    # calculate Axes Y coordinates of box top+bottom
    
    yl0, yl1 = ycoord, ycoord + dy 

    
    # calculate Axes Y distance of ticks + label margins
    y_margin = (yl1-yl0)*ytick_label_margins
    
    
    transform=ax.transAxes
    
    
    # fill black/white 'stripes' and draw their boundaries
    fill_colors = ['black', 'white']
    i_color = 0
    
    filled_boxs = []
    for xi0, xi1 in zip(xcoords[:-1],xcoords[1:]):
        
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
             bbox_kwargs = bbox_kwargs ,
             paddings =paddings)
    
    
    
    # add short tick lines
    for x in xcoords:
        plt.plot((x, x), (yl0, yl0-y_margin), 'black', 
                 transform=transform,
                 zorder=zorder,
                 clip_on=False)
    
    
    
    # add a scale legend 'Km'
    font_props = mfonts.FontProperties(size=fontsize, 
                                       weight=font_weight)
    
    plt.text(
        0.5 * (xl0 + xl1),
        yl1 + y_margin,
        'Km',
        color='k',
        verticalalignment='bottom',
        horizontalalignment='center',
        fontproperties=font_props,
        transform=transform,
        clip_on=False,
        zorder=zorder)

    # add numeric labels
    for x, xlabel in zip(xcoords, xlabels):
        plt.text(x,
                 yl0 - 2 * y_margin,
                 '{:g}'.format((xlabel) * 0.001),
                 verticalalignment='top',
                 horizontalalignment='center',
                 fontproperties=font_props,
                 transform=transform,
                 rotation=rotation,
                 clip_on=False,
                 zorder=zorder+1,
                #bbox=dict(facecolor='red', alpha=0.5) # this would add a box only around the xticks
                )
    
    
    # Adjusting figure borders to ensure that the scalebar is within its limits

    #ax.get_figure().canvas.draw()
    #ax.get_figure().tight_layout() 

    # get rectangle background bbox
    
    
    
    
    if add_numeric_scale_bar:
    
        _add_numeric_scale_bar(ax, rect, numeric_scale_bar_kwgs, fontprops=font_props)