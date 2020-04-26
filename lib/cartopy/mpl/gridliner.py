# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from __future__ import (absolute_import, division, print_function)

import operator
import warnings

import matplotlib
import matplotlib.collections as mcollections
import matplotlib.ticker as mticker
import matplotlib.transforms as mtrans
import matplotlib.path as mpath
import numpy as np
import shapely.geometry as sgeom

import cartopy
from cartopy.crs import Projection, _RectangularProjection
from cartopy.mpl.ticker import (
    LongitudeLocator, LatitudeLocator,
    LongitudeFormatter, LatitudeFormatter)

degree_locator = mticker.MaxNLocator(nbins=9, steps=[1, 1.5, 1.8, 2, 3, 6, 10])
classic_locator = mticker.MaxNLocator(nbins=9)
classic_formatter = mticker.ScalarFormatter

_DEGREE_SYMBOL = u'\u00B0'
_X_INLINE_PROJS = (
    cartopy.crs.InterruptedGoodeHomolosine,
    cartopy.crs.LambertConformal,
    cartopy.crs.Mollweide,
    cartopy.crs.Sinusoidal,
    cartopy.crs.RotatedPole,
)
_POLAR_PROJS = (
    cartopy.crs.NorthPolarStereo,
    cartopy.crs.SouthPolarStereo,
    cartopy.crs.Stereographic
)


def _fix_lons(lons):
    """
    Fix the given longitudes into the range ``[-180, 180]``.

    """
    lons = np.array(lons, copy=False, ndmin=1)
    fixed_lons = ((lons + 180) % 360) - 180
    # Make the positive 180s positive again.
    fixed_lons[(fixed_lons == -180) & (lons > 0)] *= -1
    return fixed_lons


def _lon_hemisphere(longitude):
    """Return the hemisphere (E, W or '' for 0) for the given longitude."""
    longitude = _fix_lons(longitude)
    if longitude > 0:
        hemisphere = 'E'
    elif longitude < 0:
        hemisphere = 'W'
    else:
        hemisphere = ''
    return hemisphere


def _lat_hemisphere(latitude):
    """Return the hemisphere (N, S or '' for 0) for the given latitude."""
    if latitude > 0:
        hemisphere = 'N'
    elif latitude < 0:
        hemisphere = 'S'
    else:
        hemisphere = ''
    return hemisphere


def _east_west_formatted(longitude, num_format='g'):
    fmt_string = u'{longitude:{num_format}}{degree}{hemisphere}'
    return fmt_string.format(longitude=abs(longitude), num_format=num_format,
                             hemisphere=_lon_hemisphere(longitude),
                             degree=_DEGREE_SYMBOL)


def _north_south_formatted(latitude, num_format='g'):
    fmt_string = u'{latitude:{num_format}}{degree}{hemisphere}'
    return fmt_string.format(latitude=abs(latitude), num_format=num_format,
                             hemisphere=_lat_hemisphere(latitude),
                             degree=_DEGREE_SYMBOL)


#: A formatter which turns longitude values into nice longitudes such as 110W
LONGITUDE_FORMATTER = mticker.FuncFormatter(lambda v, pos:
                                            _east_west_formatted(v))
#: A formatter which turns longitude values into nice longitudes such as 45S
LATITUDE_FORMATTER = mticker.FuncFormatter(lambda v, pos:
                                           _north_south_formatted(v))


class Gridliner(object):
    # NOTE: In future, one of these objects will be add-able to a GeoAxes (and
    # maybe even a plain old mpl axes) and it will call the "_draw_gridliner"
    # method on draw. This will enable automatic gridline resolution
    # determination on zoom/pan.
    def __init__(self, axes, crs, draw_labels=False, xlocator=None,
                 ylocator=None, collection_kwargs=None,
                 xformatter=None, yformatter=None, dms=False,
                 x_inline=None, y_inline=None, auto_inline=True):
        """
        Object used by :meth:`cartopy.mpl.geoaxes.GeoAxes.gridlines`
        to add gridlines and tick labels to a map.

        Parameters
        ----------
        axes
            The :class:`cartopy.mpl.geoaxes.GeoAxes` object to be drawn on.
        crs
            The :class:`cartopy.crs.CRS` defining the coordinate system that
            the gridlines are drawn in.
        draw_labels: optional
            Toggle whether to draw labels. For finer control, attributes of
            :class:`Gridliner` may be modified individually. Defaults to False.
        xlocator: optional
            A :class:`matplotlib.ticker.Locator` instance which will be used
            to determine the locations of the gridlines in the x-coordinate of
            the given CRS. Defaults to None, which implies automatic locating
            of the gridlines.
        ylocator: optional
            A :class:`matplotlib.ticker.Locator` instance which will be used
            to determine the locations of the gridlines in the y-coordinate of
            the given CRS. Defaults to None, which implies automatic locating
            of the gridlines.
        xformatter: optional
            A :class:`matplotlib.ticker.Formatter` instance to format labels
            for x-coordinate gridlines. It defaults to None, which implies the
            use of a :class:`cartopy.mpl.ticker.LongitudeFormatter` initiated
            with the ``dms`` argument, if the crs is of
            :class:`~cartopy.crs.PlateCarree` type.
        yformatter: optional
            A :class:`matplotlib.ticker.Formatter` instance to format labels
            for y-coordinate gridlines. It defaults to None, which implies the
            use of a :class:`cartopy.mpl.ticker.LatitudeFormatter` initiated
            with the ``dms`` argument, if the crs is of
            :class:`~cartopy.crs.PlateCarree` type.
        collection_kwargs: optional
            Dictionary controlling line properties, passed to
            :class:`matplotlib.collections.Collection`. Defaults to None.
        dms: bool
            When default locators and formatters are used,
            ticks are able to stop on minutes and seconds if minutes
            is set to True, and not fraction of degrees.
        x_inline: optional
            Toggle whether the x labels drawn should be inline.
        y_inline: optional
            Toggle whether the y labels drawn should be inline.
        auto_inline: optional
            Set x_inline and y_inline automatically based on projection.

        Notes
        -----
        The "x" and "y" labels for locators and formatters do not necessarily
        correspond to X and Y, but to the first and second coordinates of the
        specified CRS. For the common case of PlateCarree gridlines, these
        correspond to longitudes and latitudes. Depending on the projection
        used for the map, meridians and parallels can cross both the X axis and
        the Y axis.
        """
        self.axes = axes

        #: The :class:`~matplotlib.ticker.Locator` to use for the x
        #: gridlines and labels.
        if xlocator is not None:
            if not isinstance(xlocator, mticker.Locator):
                xlocator = mticker.FixedLocator(xlocator)
            self.xlocator = xlocator
        elif isinstance(crs, cartopy.crs.PlateCarree):
            self.xlocator = LongitudeLocator(dms=dms)
        else:
            self.xlocator = classic_locator

        #: The :class:`~matplotlib.ticker.Locator` to use for the y
        #: gridlines and labels.
        if ylocator is not None:
            if not isinstance(ylocator, mticker.Locator):
                ylocator = mticker.FixedLocator(ylocator)
            self.ylocator = ylocator
        elif isinstance(crs, cartopy.crs.PlateCarree):
            self.ylocator = LatitudeLocator(dms=dms)
        else:
            self.ylocator = classic_locator

        if xformatter is None:
            if isinstance(crs, cartopy.crs.PlateCarree):
                xformatter = LongitudeFormatter(dms=dms)
            else:
                xformatter = classic_formatter()
        #: The :class:`~matplotlib.ticker.Formatter` to use for the lon labels.
        self.xformatter = xformatter

        if yformatter is None:
            if isinstance(crs, cartopy.crs.PlateCarree):
                yformatter = LatitudeFormatter(dms=dms)
            else:
                yformatter = classic_formatter()
        #: The :class:`~matplotlib.ticker.Formatter` to use for the lat labels.
        self.yformatter = yformatter

        #: Whether to draw labels on the top of the map.
        self.top_labels = draw_labels

        #: Whether to draw labels on the bottom of the map.
        self.bottom_labels = draw_labels

        #: Whether to draw labels on the left hand side of the map.
        self.left_labels = draw_labels

        #: Whether to draw labels on the right hand side of the map.
        self.right_labels = draw_labels

        if auto_inline:
            if isinstance(self.axes.projection, _X_INLINE_PROJS):
                self.x_inline = True
                self.y_inline = False
            elif isinstance(self.axes.projection, _POLAR_PROJS):
                self.x_inline = False
                self.y_inline = True
            else:
                self.x_inline = False
                self.y_inline = False

        # overwrite auto_inline if necessary
        if x_inline is not None:
            #: Whether to draw x labels inline
            self.x_inline = x_inline
        elif not auto_inline:
            self.x_inline = False

        if y_inline is not None:
            #: Whether to draw y labels inline
            self.y_inline = y_inline
        elif not auto_inline:
            self.y_inline = False

        #: Whether to draw the x gridlines.
        self.xlines = True

        #: Whether to draw the y gridlines.
        self.ylines = True

        #: A dictionary passed through to ``ax.text`` on x label creation
        #: for styling of the text labels.
        self.xlabel_style = {}

        #: A dictionary passed through to ``ax.text`` on y label creation
        #: for styling of the text labels.
        self.ylabel_style = {}

        #: The padding from the map edge to the x labels in points.
        self.xpadding = 5

        #: The padding from the map edge to the y labels in points.
        self.ypadding = 5

        #: Allow the rotation of labels.
        self.rotate_labels = True

        # Current transform
        self.crs = crs

        # if the user specifies tick labels at this point, check if they can
        # be drawn. The same check will take place at draw time in case
        # public attributes are changed after instantiation.
        if draw_labels and not (x_inline or y_inline or auto_inline):
            self._assert_can_draw_ticks()

        #: The number of interpolation points which are used to draw the
        #: gridlines.
        self.n_steps = 100

        #: A dictionary passed through to
        #: ``matplotlib.collections.LineCollection`` on grid line creation.
        self.collection_kwargs = collection_kwargs

        #: The x gridlines which were created at draw time.
        self.xline_artists = []

        #: The y gridlines which were created at draw time.
        self.yline_artists = []

        # Plotted status
        self._plotted = False

        # Check visibility of labels at each draw event
        # (or once drawn, only at resize event ?)
        self.axes.figure.canvas.mpl_connect('draw_event', self._draw_event)

    @property
    def xlabels_top(self):
        warnings.warn('The .xlabels_top attribute is deprecated. Please '
                      'use .top_labels to toggle visibility instead.')
        return self.top_labels

    @xlabels_top.setter
    def xlabels_top(self, value):
        warnings.warn('The .xlabels_top attribute is deprecated. Please '
                      'use .top_labels to toggle visibility instead.')
        self.top_labels = value

    @property
    def xlabels_bottom(self):
        warnings.warn('The .xlabels_bottom attribute is deprecated. Please '
                      'use .bottom_labels to toggle visibility instead.')
        return self.bottom_labels

    @xlabels_bottom.setter
    def xlabels_bottom(self, value):
        warnings.warn('The .xlabels_bottom attribute is deprecated. Please '
                      'use .bottom_labels to toggle visibility instead.')
        self.bottom_labels = value

    @property
    def ylabels_left(self):
        warnings.warn('The .ylabels_left attribute is deprecated. Please '
                      'use .left_labels to toggle visibility instead.')
        return self.left_labels

    @ylabels_left.setter
    def ylabels_left(self, value):
        warnings.warn('The .ylabels_left attribute is deprecated. Please '
                      'use .left_labels to toggle visibility instead.')
        self.left_labels = value

    @property
    def ylabels_right(self):
        warnings.warn('The .ylabels_right attribute is deprecated. Please '
                      'use .right_labels to toggle visibility instead.')
        return self.right_labels

    @ylabels_right.setter
    def ylabels_right(self, value):
        warnings.warn('The .ylabels_right attribute is deprecated. Please '
                      'use .right_labels to toggle visibility instead.')
        self.right_labels = value

    def _draw_event(self, event):
        if self.has_labels():
            self._update_labels_visibility(event.renderer)

    def has_labels(self):
        return hasattr(self, '_labels') and self._labels

    @property
    def label_artists(self):
        if self.has_labels():
            return self._labels
        return []

    def _crs_transform(self):
        """
        Get the drawing transform for our gridlines.

        Note
        ----
            The drawing transform depends on the transform of our 'axes', so
            it may change dynamically.

        """
        transform = self.crs
        if not isinstance(transform, mtrans.Transform):
            transform = transform._as_mpl_transform(self.axes)
        return transform

    @staticmethod
    def _round(x, base=5):
        if np.isnan(base):
            base = 5
        return int(base * round(float(x) / base))

    def _find_midpoints(self, lim, ticks):
        # Find the center point between each lat gridline.
        if len(ticks) > 1:
            cent = np.diff(ticks).mean() / 2
        else:
            cent = np.nan
        if isinstance(self.axes.projection, _POLAR_PROJS):
            lq = 90
            uq = 90
        else:
            lq = 25
            uq = 75
        midpoints = (self._round(np.percentile(lim, lq), cent),
                     self._round(np.percentile(lim, uq), cent))
        return midpoints

    def _draw_gridliner(self, nx=None, ny=None, renderer=None):
        """Create Artists for all visible elements and add to our Axes."""
        # Check status
        if self._plotted:
            return
        self._plotted = True

        # Inits
        lon_lim, lat_lim = self._axes_domain(nx=nx, ny=ny)

        transform = self._crs_transform()
        rc_params = matplotlib.rcParams
        n_steps = self.n_steps
        crs = self.crs

        # Get nice ticks within crs domain
        lon_ticks = self.xlocator.tick_values(lon_lim[0], lon_lim[1])
        lat_ticks = self.ylocator.tick_values(lat_lim[0], lat_lim[1])
        lon_ticks = [value for value in lon_ticks
                     if value >= max(lon_lim[0], crs.x_limits[0]) and
                     value <= min(lon_lim[1], crs.x_limits[1])]
        lat_ticks = [value for value in lat_ticks
                     if value >= max(lat_lim[0], crs.y_limits[0]) and
                     value <= min(lat_lim[1], crs.y_limits[1])]

        #####################
        # Gridlines drawing #
        #####################

        collection_kwargs = self.collection_kwargs
        if collection_kwargs is None:
            collection_kwargs = {}
        collection_kwargs = collection_kwargs.copy()
        collection_kwargs['transform'] = transform
        # XXX doesn't gracefully handle lw vs linewidth aliases...
        collection_kwargs.setdefault('color', rc_params['grid.color'])
        collection_kwargs.setdefault('linestyle', rc_params['grid.linestyle'])
        collection_kwargs.setdefault('linewidth', rc_params['grid.linewidth'])

        # Meridians
        lat_min, lat_max = lat_lim
        if lat_ticks:
            lat_min = min(lat_min, min(lat_ticks))
            lat_max = max(lat_max, max(lat_ticks))
        lon_lines = np.empty((len(lon_ticks), n_steps, 2))
        lon_lines[:, :, 0] = np.array(lon_ticks)[:, np.newaxis]
        lon_lines[:, :, 1] = np.linspace(lat_min, lat_max,
                                         n_steps)[np.newaxis, :]

        if self.xlines:
            nx = len(lon_lines) + 1
            # XXX this bit is cartopy specific. (for circular longitudes)
            # Purpose: omit plotting the last x line,
            # as it may overlap the first.
            if (isinstance(crs, Projection) and
                    isinstance(crs, _RectangularProjection) and
                    abs(np.diff(lon_lim)) == abs(np.diff(crs.x_limits))):
                nx -= 1
            lon_lc = mcollections.LineCollection(lon_lines,
                                                 **collection_kwargs)
            self.xline_artists.append(lon_lc)
            self.axes.add_collection(lon_lc, autolim=False)

        # Parallels
        lon_min, lon_max = lon_lim
        if lon_ticks:
            lon_min = min(lon_min, min(lon_ticks))
            lon_max = max(lon_max, max(lon_ticks))
        lat_lines = np.empty((len(lat_ticks), n_steps, 2))
        lat_lines[:, :, 0] = np.linspace(lon_min, lon_max,
                                         n_steps)[np.newaxis, :]
        lat_lines[:, :, 1] = np.array(lat_ticks)[:, np.newaxis]
        if self.ylines:
            lat_lc = mcollections.LineCollection(lat_lines,
                                                 **collection_kwargs)
            self.yline_artists.append(lat_lc)
            self.axes.add_collection(lat_lc, autolim=False)

        #################
        # Label drawing #
        #################

        self.bottom_label_artists = []
        self.top_label_artists = []
        self.left_label_artists = []
        self.right_label_artists = []
        if not (self.left_labels or self.right_labels or
                self.bottom_labels or self.top_labels):
            return
        self._assert_can_draw_ticks()

        # Get the real map boundaries
        map_boundary_vertices = self.axes.patch.get_path().vertices
        map_boundary = sgeom.Polygon(map_boundary_vertices)

        self._labels = []

        if self.x_inline:
            y_midpoints = self._find_midpoints(lat_lim, lat_ticks)
        if self.y_inline:
            x_midpoints = self._find_midpoints(lon_lim, lon_ticks)

        for lonlat, lines, line_ticks, formatter, label_style in (
                ('lon', lon_lines, lon_ticks,
                 self.xformatter, self.xlabel_style),
                ('lat', lat_lines, lat_ticks,
                 self.yformatter, self.ylabel_style)):

            formatter.set_locs(line_ticks)

            for line, tick_value in zip(lines, line_ticks):
                # Intersection of line with map boundary
                line = self.axes.projection.transform_points(
                    crs, line[:, 0], line[:, 1])[:, :2]
                infs = np.isinf(line).any(axis=1)
                line = line.compress(~infs, axis=0)
                if line.size == 0:
                    continue
                line = sgeom.LineString(line)
                if line.intersects(map_boundary):
                    intersection = line.intersection(map_boundary)
                    del line
                    if intersection.is_empty:
                        continue
                    if isinstance(intersection, sgeom.MultiPoint):
                        if len(intersection) < 2:
                            continue
                        tails = [[(pt.x, pt.y) for pt in intersection[:2]]]
                        heads = [[(pt.x, pt.y)
                                  for pt in intersection[-1:-3:-1]]]
                    elif isinstance(intersection, (sgeom.LineString,
                                                   sgeom.MultiLineString)):
                        if isinstance(intersection, sgeom.LineString):
                            intersection = [intersection]
                        elif len(intersection) > 4:
                            # Gridline and map boundary are parallel
                            # and they intersect themselves too much
                            # it results in a multiline string
                            # that must be converted to a single linestring.
                            # This is an empirical workaround for a problem
                            # that can probably be solved in a cleaner way.
                            xy = np.append(intersection[0], intersection[-1],
                                           axis=0)
                            intersection = [sgeom.LineString(xy)]
                        tails = []
                        heads = []
                        for inter in intersection:
                            if len(inter.coords) < 2:
                                continue
                            tails.append(inter.coords[:2])
                            heads.append(inter.coords[-1:-3:-1])
                        if not tails:
                            continue
                    elif isinstance(intersection,
                                    sgeom.collection.GeometryCollection):
                        # This is a collection of Point and LineString that
                        # represent the same gridline.
                        # We only consider the first geometries, merge their
                        # coordinates and keep first two points to get only one
                        # tail ...
                        xy = []
                        for geom in intersection.geoms:
                            for coord in geom.coords:
                                xy.append(coord)
                                if len(xy) == 2:
                                    break
                            if len(xy) == 2:
                                break
                        tails = [xy]
                        # ... and the last geometries, merge their coordinates
                        # and keep last two points to get only one head.
                        xy = []
                        for geom in reversed(intersection.geoms):
                            for coord in reversed(geom.coords):
                                xy.append(coord)
                                if len(xy) == 2:
                                    break
                            if len(xy) == 2:
                                break
                        heads = [xy]
                    else:
                        warnings.warn(
                            'Unsupported intersection geometry for gridline '
                            'labels: ' + intersection.__class__.__name__)
                        continue
                    del intersection

                    # Loop on head and tail and plot label by extrapolation
                    for tail, head in zip(tails, heads):
                        for i, (pt0, pt1) in enumerate([tail, head]):
                            kw, angle, loc = self._segment_to_text_specs(
                                pt0, pt1, lonlat)
                            if not getattr(self, loc+'_labels'):
                                continue
                            kw.update(label_style,
                                      bbox={'pad': 0, 'visible': False})
                            text = formatter(tick_value)

                            if self.y_inline and lonlat == 'lat':
                                # 180 degrees isn't formatted with a
                                # suffix and adds confusion if it's inline
                                if abs(tick_value) == 180:
                                    continue
                                x = x_midpoints[i]
                                y = tick_value
                                kw.update(clip_on=True)
                                y_set = True
                            else:
                                x = pt0[0]
                                y_set = False

                            if self.x_inline and lonlat == 'lon':
                                if abs(tick_value) == 180:
                                    continue
                                x = tick_value
                                y = y_midpoints[i]
                                kw.update(clip_on=True)
                            elif not y_set:
                                y = pt0[1]

                            tt = self.axes.text(x, y, text, **kw)
                            tt._angle = angle
                            priority = (((lonlat == 'lon') and
                                         loc in ('bottom', 'top')) or
                                        ((lonlat == 'lat') and
                                         loc in ('left', 'right')))
                            self._labels.append((lonlat, priority, tt))
                            getattr(self, loc + '_label_artists').append(tt)

        # Sort labels
        if self._labels:
            self._labels.sort(key=operator.itemgetter(0), reverse=True)
            self._update_labels_visibility(renderer)

    def _segment_to_text_specs(self, pt0, pt1, lonlat):
        """Get appropriate kwargs for a label from lon or lat line segment"""
        x0, y0 = pt0
        x1, y1 = pt1
        angle = np.arctan2(y0-y1, x0-x1) * 180 / np.pi
        kw, loc = self._segment_angle_to_text_specs(angle, lonlat)
        return kw, angle, loc

    def _text_angle_to_specs_(self, angle, lonlat):
        """Get specs for a rotated label from its angle in degrees"""

        angle %= 360
        if angle > 180:
            angle -= 360

        if ((self.x_inline and lonlat == 'lon') or
                (self.y_inline and lonlat == 'lat')):
            kw = {'rotation': 0, 'rotation_mode': 'anchor',
                  'ha': 'center', 'va': 'center'}
            loc = 'bottom'
            return kw, loc

        # Default options
        kw = {'rotation': angle, 'rotation_mode': 'anchor'}

        # Options that depend in which quarter the angle falls
        if abs(angle) <= 45:
            loc = 'right'
            kw.update(ha='left', va='center')

        elif abs(angle) >= 135:
            loc = 'left'
            kw.update(ha='right', va='center')
            kw['rotation'] -= np.sign(angle) * 180

        elif angle > 45:
            loc = 'top'
            kw.update(ha='center', va='bottom', rotation=angle-90)

        else:
            loc = 'bottom'
            kw.update(ha='center', va='top', rotation=angle+90)

        return kw, loc

    def _segment_angle_to_text_specs(self, angle, lonlat):
        """Get appropriate kwargs for a given text angle"""
        kw, loc = self._text_angle_to_specs_(angle, lonlat)
        if not self.rotate_labels:
            angle = {'top': 90., 'right': 0.,
                     'bottom': -90., 'left': 180.}[loc]
            del kw['rotation']

        if ((self.x_inline and lonlat == 'lon') or
                (self.y_inline and lonlat == 'lat')):
            kw.update(transform=cartopy.crs.PlateCarree())
        else:
            xpadding = (self.xpadding if self.xpadding is not None
                        else matplotlib.rc_params['xtick.major.pad'])
            ypadding = (self.ypadding if self.ypadding is not None
                        else matplotlib.rc_params['ytick.major.pad'])
            dx = ypadding * np.cos(angle * np.pi / 180)
            dy = xpadding * np.sin(angle * np.pi / 180)
            transform = mtrans.offset_copy(
                self.axes.transData, self.axes.figure,
                x=dx, y=dy, units='dots')
            kw.update(transform=transform)

        return kw, loc

    def _update_labels_visibility(self, renderer):
        """Update the visibility of each plotted label

        The following rules apply:

        - Labels are plotted and checked by order of priority,
          with a high priority for longitude labels at the bottom and
          top of the map, and the reverse for latitude labels.
        - A label must not overlap another label marked as visible.
        - A label must not overlap the map boundary.
        - When a label is about to be hidden, other angles are tried in the
          absolute given limit of max_delta_angle by increments of delta_angle
          of difference from the original angle.
        """
        if renderer is None or not self._labels:
            return
        paths = []
        outline_path = None
        delta_angle = 22.5
        max_delta_angle = 45
        axes_children = self.axes.get_children()

        def remove_path_dupes(path):
            """
            Remove duplicate points in a path (zero-length segments).

            This is necessary only for Matplotlib 3.1.0 -- 3.1.2, because
            Path.intersects_path incorrectly returns True for any paths with
            such segments.
            """
            segment_length = np.diff(path.vertices, axis=0)
            mask = np.logical_or.reduce(segment_length != 0, axis=1)
            mask = np.append(mask, True)
            path = mpath.Path(np.compress(mask, path.vertices, axis=0),
                              np.compress(mask, path.codes, axis=0))
            return path

        for lonlat, priority, artist in self._labels:

            if artist not in axes_children:
                warnings.warn('The labels of this gridliner do not belong to '
                              'the gridliner axes')

            orig_specs = {'rotation': artist.get_rotation(),
                          'ha': artist.get_ha(),
                          'va': artist.get_va()}
            # Compute angles to try
            angles = [None]
            for abs_delta_angle in np.arange(delta_angle, max_delta_angle+1,
                                             delta_angle):
                angles.append(artist._angle + abs_delta_angle)
                angles.append(artist._angle - abs_delta_angle)

            # Loop on angles until it works
            for angle in angles:
                if ((self.x_inline and lonlat == 'lon') or
                        (self.y_inline and lonlat == 'lat')):
                    angle = 0

                if angle is not None:
                    specs, _ = self._segment_angle_to_text_specs(angle, lonlat)
                    artist.update(specs)

                artist.update_bbox_position_size(renderer)
                this_patch = artist.get_bbox_patch()
                this_path = this_patch.get_path().transformed(
                    this_patch.get_transform())
                if '3.1.0' <= matplotlib.__version__ <= '3.1.2':
                    this_path = remove_path_dupes(this_path)
                center = artist.get_transform().transform_point(
                    artist.get_position())
                visible = False

                for path in paths:

                    # Check it does not overlap another label
                    if this_path.intersects_path(path):
                        break

                else:

                    # Finally check that it does not overlap the map
                    if outline_path is None:
                        outline_path = (self.axes.patch.get_path()
                                        .transformed(self.axes.transData))
                        if '3.1.0' <= matplotlib.__version__ <= '3.1.2':
                            outline_path = remove_path_dupes(outline_path)
                    # Inline must be within the map.
                    if ((lonlat == 'lon' and self.x_inline) or
                            (lonlat == 'lat' and self.y_inline)):
                        # TODO: When Matplotlib clip path works on text, this
                        # clipping can be left to it.
                        if outline_path.contains_point(center):
                            visible = True
                    # Non-inline must not run through the outline.
                    elif not outline_path.intersects_path(this_path):
                        visible = True

                    # Good
                    if visible:
                        break

                if ((self.x_inline and lonlat == 'lon') or
                        (self.y_inline and lonlat == 'lat')):
                    break

            # Action
            artist.set_visible(visible)
            if not visible:
                artist.update(orig_specs)
            else:
                paths.append(this_path)

    def _assert_can_draw_ticks(self):
        """
        Check to see if ticks can be drawn. Either returns True or raises
        an exception.

        """
        # Check labelling is supported, currently a limited set of options.
        if not isinstance(self.crs, cartopy.crs.PlateCarree):
            raise TypeError('Cannot label {crs.__class__.__name__} gridlines.'
                            ' Only PlateCarree gridlines are currently '
                            'supported.'.format(crs=self.crs))
        return True

    def _axes_domain(self, nx=None, ny=None):
        """Return lon_range, lat_range"""
        DEBUG = False

        transform = self._crs_transform()

        ax_transform = self.axes.transAxes
        desired_trans = ax_transform - transform

        nx = nx or 100
        ny = ny or 100
        x = np.linspace(1e-9, 1 - 1e-9, nx)
        y = np.linspace(1e-9, 1 - 1e-9, ny)
        x, y = np.meshgrid(x, y)

        coords = np.column_stack((x.ravel(), y.ravel()))

        in_data = desired_trans.transform(coords)

        ax_to_bkg_patch = self.axes.transAxes - self.axes.patch.get_transform()

        # convert the coordinates of the data to the background patches
        # coordinates
        background_coord = ax_to_bkg_patch.transform(coords)
        ok = self.axes.patch.get_path().contains_points(background_coord)

        if DEBUG:
            import matplotlib.pyplot as plt
            plt.plot(coords[ok, 0], coords[ok, 1], 'or',
                     clip_on=False, transform=ax_transform)
            plt.plot(coords[~ok, 0], coords[~ok, 1], 'ob',
                     clip_on=False, transform=ax_transform)

        inside = in_data[ok, :]

        # If there were no data points in the axes we just use the x and y
        # range of the projection.
        if inside.size == 0:
            lon_range = self.crs.x_limits
            lat_range = self.crs.y_limits
        else:
            # np.isfinite must be used to prevent np.inf values that
            # not filtered by np.nanmax for some projections
            lat_max = np.compress(np.isfinite(inside[:, 1]),
                                  inside[:, 1])
            if lat_max.size == 0:
                lon_range = self.crs.x_limits
                lat_range = self.crs.y_limits
            else:
                lat_max = lat_max.max()
                lon_range = np.nanmin(inside[:, 0]), np.nanmax(inside[:, 0])
                lat_range = np.nanmin(inside[:, 1]), lat_max

        # XXX Cartopy specific thing. Perhaps make this bit a specialisation
        # in a subclass...
        crs = self.crs
        if isinstance(crs, Projection):
            lon_range = np.clip(lon_range, *crs.x_limits)
            lat_range = np.clip(lat_range, *crs.y_limits)

            # if the limit is >90% of the full x limit, then just use the full
            # x limit (this makes circular handling better)
            prct = np.abs(np.diff(lon_range) / np.diff(crs.x_limits))
            if prct > 0.9:
                lon_range = crs.x_limits

        return lon_range, lat_range
