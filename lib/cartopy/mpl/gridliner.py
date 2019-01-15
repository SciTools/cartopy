# (C) British Crown Copyright 2011 - 2019, Met Office
#
# This file is part of cartopy.
#
# cartopy is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cartopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cartopy.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division, print_function)

import operator
from array import array

import matplotlib
import matplotlib.collections as mcollections
import matplotlib.ticker as mticker
import matplotlib.transforms as mtrans
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry as sgeom
from warnings import warn

import cartopy
from cartopy.crs import Projection, _RectangularProjection
from cartopy.mpl.ticker import (
        LongitudeLocator, LatitudeLocator,
        LongitudeFormatter, LatitudeFormatter)

degree_locator = mticker.MaxNLocator(nbins=9, steps=[1, 1.5, 1.8, 2, 3, 6, 10])
classic_locator = mticker.MaxNLocator(nbins=9)

_DEGREE_SYMBOL = u'\u00B0'


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


def _text_angle_to_specs_(angle):
    """Get appropriate kwargs for a rotated label from its angle in degrees"""

    if matplotlib.__version__ >= '3.1':
        # rotation_mode='anchor' and va_align_center='center_baseline'
        # are incompatible before mpl-3.1
        va_align_center = 'center_baseline'
    else:
        va_align_center = 'center'

    angle %= 360
    if angle > 180:
        angle -= 360

    # Default options
    kw = {'rotation': angle, 'rotation_mode': 'anchor'}

    # Options that depend in which quarter the angle falls
    if abs(angle) <= 45:

        loc = 'right'
        kw.update(ha='left', va=va_align_center)

    elif abs(angle) >= 135:

        loc = 'left'
        kw.update(ha='right', va=va_align_center)
        kw['rotation'] -= np.sign(angle) * 180

    elif angle > 45:

        loc = 'top'
        kw.update(ha='center', va='bottom', rotation=angle-90)

    else:

        loc = 'bottom'
        kw.update(ha='center', va='top', rotation=angle+90)

    return kw, loc


class Gridliner(object):
    # NOTE: In future, one of these objects will be add-able to a GeoAxes (and
    # maybe even a plain old mpl axes) and it will call the "_draw_gridliner"
    # method on draw. This will enable automatic gridline resolution
    # determination on zoom/pan.
    def __init__(self, axes, crs, draw_labels=False, xlocator=None,
                 ylocator=None, collection_kwargs=None,
                 xformatter=None, yformatter=None,
                 minutes=False):
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
            A :class:`matplotlib.ticker.Formatter` instance to format
            longitude labels. Defaults to
            :class:`cartopy.mpl.ticker.LongitudeFormatter`.
        yformatter: optional
            A :class:`matplotlib.ticker.Formatter` instance to format
            latitude labels. Defaults to
            :class:`cartopy.mpl.ticker.LatitudeFormatter`.
        collection_kwargs: optional
            Dictionary controlling line properties, passed to
            :class:`matplotlib.collections.Collection`. Defaults to None.
        minutes: bool
            When default locators and formatters are used,
            ticks are able to stop on minutes and seconds if minutes
            is set to True, and not fraction of degrees.

        .. note:: Note that the "x" and "y" labels for locators and
             formatters do not necessarily correspond to X and Y,
             but to longitudes and latitudes: indeed, according to
             geographical projections, meridians and parallels can
             cross both the X axis and the Y axis.

        """
        self.axes = axes

        #: The :class:`~matplotlib.ticker.Locator` to use for the x
        #: gridlines and labels.
        if xlocator is not None:
            if not isinstance(xlocator, mticker.Locator):
                xlocator = mticker.FixedLocator(xlocator)
            self.xlocator = xlocator
        elif isinstance(crs, cartopy.crs.PlateCarree):
            self.xlocator = LongitudeLocator(minutes=minutes)
        else:
            self.xlocator = classic_locator

        #: The :class:`~matplotlib.ticker.Locator` to use for the y
        #: gridlines and labels.
        if ylocator is not None:
            if not isinstance(ylocator, mticker.Locator):
                ylocator = mticker.FixedLocator(ylocator)
            self.ylocator = ylocator
        elif isinstance(crs, cartopy.crs.PlateCarree):
            self.ylocator = LatitudeLocator(minutes=minutes)
        else:
            self.ylocator = classic_locator

        #: The :class:`~matplotlib.ticker.Formatter` to use for the lon labels.
        self.xformatter = xformatter or LongitudeFormatter(minutes=minutes)

        #: The :class:`~matplotlib.ticker.Formatter` to use for the lat labels.
        self.yformatter = yformatter or LatitudeFormatter(minutes=minutes)

        #: Whether to draw labels on the top of the map.
        self.top_labels = draw_labels

        #: Whether to draw labels on the bottom of the map.
        self.bottom_labels = draw_labels

        #: Whether to draw labels on the left hand side of the map.
        self.left_labels = draw_labels

        #: Whether to draw labels on the right hand side of the map.
        self.right_labels = draw_labels

        #: Whether to draw the longitude gridlines (meridians).
        self.xlines = True

        #: Whether to draw the latitude gridlines (parallels).
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
        if draw_labels:
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

    def _draw_gridliner(self, nx=None, ny=None, background_patch=None,
                        renderer=None):
        """Create Artists for all visible elements and add to our Axes."""
        # Check status
        if self._plotted:
            return
        self._plotted = True

        # Inits
        lon_lim, lat_lim = self._axes_domain(
                nx=nx, ny=ny, background_patch=background_patch)

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
        lon_lines = []
        lat_min, lat_max = lat_lim
        if lat_ticks:
            lat_min = min(lat_min, min(lat_ticks))
            lat_max = max(lat_max, max(lat_ticks))
        for x in lon_ticks:
            ticks = list(zip(
                [x]*n_steps,
                np.linspace(lat_min, lat_max, n_steps)))
            lon_lines.append(ticks)

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
        lat_lines = []
        lon_min, lon_max = lon_lim
        if lon_ticks:
            lon_min = min(lon_min, min(lon_ticks))
            lon_max = max(lon_max, max(lon_ticks))
        for y in lat_ticks:
            ticks = list(zip(
                np.linspace(lon_min, lon_max, n_steps),
                [y]*n_steps))
            lat_lines.append(ticks)
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
        map_boundary_vertices = self.axes.background_patch.get_path().vertices
        map_boundary = sgeom.Polygon(map_boundary_vertices)

        self._labels = []
        for lonlat, lines, line_ticks, formatter, label_style in (
                ('lon', lon_lines, lon_ticks,
                 self.xformatter, self.xlabel_style),
                ('lat', lat_lines, lat_ticks,
                 self.yformatter, self.ylabel_style)):

            formatter.set_locs(line_ticks)

            for line, tick_value in zip(lines, line_ticks):

                # Intersection of line with map boundary
                line = np.array(line)
                line = self.axes.projection.transform_points(
                        crs, line[:, 0], line[:, 1])[:, :2]
                infs = np.isinf(line)
                if infs.any():
                    if infs.all():
                        continue
                    line = line.compress(~infs.any(axis=1), axis=0)
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
                            x, y = intersection[0].coords.xy
                            x += intersection[-1].coords.xy[0]
                            y += intersection[-1].coords.xy[1]
                            xy = np.array((x, y))
                            intersection = [sgeom.LineString(xy.T)]
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
                        # represent the same gridline. We only consider
                        # the first and last geometries, merge their
                        # coordinates and keep first or last two points
                        # to get only one tail and and one head.
                        tails = []
                        heads = []
                        for ht, slicer in [(tails, slice(0, 2)),
                                           (heads, slice(-1, -3, -1))]:
                            x = array('d', [])
                            y = array('d', [])
                            for geom in list(intersection.geoms)[slicer]:
                                x += geom.xy[0]
                                y += geom.xy[1]
                            ht.append(list(zip(x[slicer], y[slicer])))
                    else:
                        warn('Unsupported intersection geometry for gridline'
                             'labels: '+intersection.__class__)
                        continue
                    del intersection

                    # Loop on head and tail and plot label by extrapolation
                    for tail, head in zip(tails, heads):
                        for pt0, pt1 in [tail, head]:
                            kw, angle, loc = self._segment_to_text_specs(
                                    pt0, pt1)
                            if not getattr(self, loc+'_labels'):
                                continue
                            kw.update(label_style,
                                      bbox={'pad': 0, 'visible': False})
                            x = round(pt0[0], 5)
                            y = round(pt0[1], 5)
                            text = formatter(tick_value)
                            tt = self.axes.text(pt0[0], pt0[1], text, **kw)
                            tt._angle = angle
                            priority = (((lonlat == 'lon') and
                                         loc in ('bottom', 'top')) or
                                        ((lonlat == 'lat') and
                                         loc in ('left', 'right')))
                            self._labels.append((priority, tt))
                            getattr(self, loc + '_label_artists').append(tt)

        # Sort labels
        if self._labels:
            self._labels.sort(key=operator.itemgetter(0), reverse=True)
            self._update_labels_visibility(renderer)

    def _segment_to_text_specs(self, pt0, pt1):
        """Get appropriate kwargs for a label from lon or lat line segment"""
        x0, y0 = pt0
        x1, y1 = pt1
        angle = np.arctan2(y0-y1, x0-x1) * 180 / np.pi
        kw, loc = self._segment_angle_to_text_specs(angle)
        return kw, angle, loc

    def _segment_angle_to_text_specs(self, angle):
        """Get appropriate kwargs for a given text angle"""
        kw, loc = _text_angle_to_specs_(angle)
        if not self.rotate_labels:
            angle = {'top': 90., 'right': 0.,
                     'bottom': -90., 'left': 180.}[loc]
            del kw['rotation']

        xpadding = (self.xpadding if self.xpadding is not None
                    else matplotlib.rc_params['xtick.major.pad'])
        ypadding = (self.ypadding if self.ypadding is not None
                    else matplotlib.rc_params['ytick.major.pad'])
        dx = ypadding * np.cos(angle * np.pi / 180)
        dy = xpadding * np.sin(angle * np.pi / 180)
        transform = mtrans.offset_copy(self.axes.transData, self.axes.figure,
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
        for priority, artist in self._labels:

            if artist not in axes_children:
                warn('The labels of this gridliner do not belong'
                     'to the gridliner axes')

            # Compute angles to try
            orig_specs = {'rotation': artist.get_rotation(),
                          'ha': artist.get_ha(),
                          'va': artist.get_va()}
            angles = [None]
            for abs_delta_angle in np.arange(delta_angle, max_delta_angle+1,
                                             delta_angle):
                for sign_delta_angle in (1, -1):
                    angle = artist._angle + sign_delta_angle * abs_delta_angle
                    angles.append(angle)

            # Loop on angles until it works
            for angle in angles:

                if angle is not None:
                    specs, _ = self._segment_angle_to_text_specs(angle)
                    plt.setp(artist, **specs)

                artist.update_bbox_position_size(renderer)
                this_patch = artist.get_bbox_patch()
                this_path = this_patch.get_path().transformed(
                        this_patch.get_transform())
                visible = False
                for path in paths:

                    # Check it does not overlap another label
                    if this_path.intersects_path(path):
                        break

                else:

                    # Finally check that it does not overlap the map
                    if outline_path is None:
                        outline_path = self.axes.background_patch.get_path(
                            ).transformed(self.axes.transData)
                    visible = not outline_path.intersects_path(this_path)

                    # Good
                    if visible:
                        break

            # Action
            artist.set_visible(visible)
            if not visible:
                plt.setp(artist, **orig_specs)
            else:
                paths.append(this_path)
                # is this necessary?
#                artist.update_bbox_position_size(renderer)

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

    def _axes_domain(self, nx=None, ny=None, background_patch=None):
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

        ax_to_bkg_patch = self.axes.transAxes - \
            background_patch.get_transform()

        # convert the coordinates of the data to the background patches
        # coordinates
        background_coord = ax_to_bkg_patch.transform(coords)
        ok = background_patch.get_path().contains_points(background_coord)

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
                                  inside[:, 1]).max()
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
