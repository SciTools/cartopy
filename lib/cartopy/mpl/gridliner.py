# (C) British Crown Copyright 2011 - 2016, Met Office
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

import matplotlib
import matplotlib.collections as mcollections
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.transforms as mtrans
import numpy as np
import six

import cartopy
from cartopy.crs import Projection, _RectangularProjection


degree_locator = mticker.MaxNLocator(nbins=9, steps=[1, 2, 3, 6, 15, 18])

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


def _lon_heimisphere(longitude):
    """Return the hemisphere (E, W or '' for 0) for the given longitude."""
    longitude = _fix_lons(longitude)
    if longitude > 0:
        hemisphere = 'E'
    elif longitude < 0:
        hemisphere = 'W'
    else:
        hemisphere = ''
    return hemisphere


def _lat_heimisphere(latitude):
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
                             hemisphere=_lon_heimisphere(longitude),
                             degree=_DEGREE_SYMBOL)


def _north_south_formatted(latitude, num_format='g'):
    fmt_string = u'{latitude:{num_format}}{degree}{hemisphere}'
    return fmt_string.format(latitude=abs(latitude), num_format=num_format,
                             hemisphere=_lat_heimisphere(latitude),
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
                 ylocator=None, collection_kwargs=None):
        """
        Object used by :meth:`cartopy.mpl.geoaxes.GeoAxes.gridlines`
        to add gridlines and tick labels to a map.

        Args:

        * axes
            The :class:`cartopy.mpl.geoaxes.GeoAxes` object to be drawn on.

        * crs
            The :class:`cartopy.crs.CRS` defining the coordinate system that
            the gridlines are drawn in.

        * draw_labels
            Toggle whether to draw labels. For finer control, attributes of
            :class:`Gridliner` may be modified individually.

        * xlocator
            A :class:`matplotlib.ticker.Locator` instance which will be used
            to determine the locations of the gridlines in the x-coordinate of
            the given CRS. Defaults to None, which implies automatic locating
            of the gridlines.

        * ylocator
            A :class:`matplotlib.ticker.Locator` instance which will be used
            to determine the locations of the gridlines in the y-coordinate of
            the given CRS. Defaults to None, which implies automatic locating
            of the gridlines.

        * collection_kwargs
            Dictionary controlling line properties, passed to
            :class:`matplotlib.collections.Collection`.

        """
        self.axes = axes

        #: The :class:`~matplotlib.ticker.Locator` to use for the x
        #: gridlines and labels.
        self.xlocator = xlocator or degree_locator

        #: The :class:`~matplotlib.ticker.Locator` to use for the y
        #: gridlines and labels.
        self.ylocator = ylocator or degree_locator

        #: The :class:`~matplotlib.ticker.Formatter` to use for the x labels.
        self.xformatter = mticker.ScalarFormatter()
        self.xformatter.create_dummy_axis()

        #: The :class:`~matplotlib.ticker.Formatter` to use for the y labels.
        self.yformatter = mticker.ScalarFormatter()
        self.yformatter.create_dummy_axis()

        #: Whether to draw labels on the top of the map.
        self.xlabels_top = draw_labels

        #: Whether to draw labels on the bottom of the map.
        self.xlabels_bottom = draw_labels

        #: Whether to draw labels on the left hand side of the map.
        self.ylabels_left = draw_labels

        #: Whether to draw labels on the right hand side of the map.
        self.ylabels_right = draw_labels

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

        self.crs = crs

        # if the user specifies tick labels at this point, check if they can
        # be drawn. The same check will take place at draw time in case
        # public attributes are changed after instantiation.
        if draw_labels:
            self._assert_can_draw_ticks()

        #: The number of interpolation points which are used to draw the
        #: gridlines.
        self.n_steps = 30

        #: A dictionary passed through to
        #: ``matplotlib.collections.LineCollection`` on grid line creation.
        self.collection_kwargs = collection_kwargs

        #: The x gridlines which were created at draw time.
        self.xline_artists = []

        #: The y gridlines which were created at draw time.
        self.yline_artists = []

        #: The x labels which were created at draw time.
        self.xlabel_artists = []

        #: The y labels which were created at draw time.
        self.ylabel_artists = []

    def _crs_transform(self):
        """
        Get the drawing transform for our gridlines.

        .. note::
            this depends on the transform of our 'axes', so it may change
            dynamically.

        """
        transform = self.crs
        if not isinstance(transform, mtrans.Transform):
            transform = transform._as_mpl_transform(self.axes)
        return transform

    def _add_gridline_label(self, value, axis, upper_end):
        """
        Create a Text artist on our axes for a gridline label.

        Args:

        * value
            Coordinate value of this gridline.  The text contains this
            value, and is positioned centred at that point.

        * axis
            which axis the label is on: 'x' or 'y'.

        * upper_end
            If True, place at the maximum of the "other" coordinate (Axes
            coordinate == 1.0).  Else 'lower' end (Axes coord = 0.0).

        """
        transform = self._crs_transform()
        shift_dist_points = 5     # A margin from the map edge.
        if upper_end is False:
            shift_dist_points = -shift_dist_points
        if axis == 'x':
            x = value
            y = 1.0 if upper_end else 0.0
            h_align = 'center'
            v_align = 'bottom' if upper_end else 'top'
            tr_x = transform
            tr_y = self.axes.transAxes + \
                mtrans.ScaledTranslation(
                    0.0,
                    shift_dist_points * (1.0 / 72),
                    self.axes.figure.dpi_scale_trans)
            str_value = self.xformatter(value)
            user_label_style = self.xlabel_style
        elif axis == 'y':
            y = value
            x = 1.0 if upper_end else 0.0
            v_align = 'center'
            h_align = 'left' if upper_end else 'right'
            tr_y = transform
            tr_x = self.axes.transAxes + \
                mtrans.ScaledTranslation(
                    shift_dist_points * (1.0 / 72),
                    0.0,
                    self.axes.figure.dpi_scale_trans)
            str_value = self.yformatter(value)
            user_label_style = self.ylabel_style
        else:
            raise ValueError(
                "Unknown axis, {!r}, must be either 'x' or 'y'".format(axis))

        # Make a 'blended' transform for label text positioning.
        # One coord is geographic, and the other a plain Axes
        # coordinate with an appropriate offset.
        label_transform = mtrans.blended_transform_factory(
            x_transform=tr_x, y_transform=tr_y)

        label_style = {'verticalalignment': v_align,
                       'horizontalalignment': h_align,
                       }
        label_style.update(user_label_style)

        # Create and add a Text artist with these properties
        text_artist = mtext.Text(x, y, str_value,
                                 clip_on=False,
                                 transform=label_transform, **label_style)
        if axis == 'x':
            self.xlabel_artists.append(text_artist)
        elif axis == 'y':
            self.ylabel_artists.append(text_artist)
        self.axes.add_artist(text_artist)

    def _draw_gridliner(self, nx=None, ny=None, background_patch=None):
        """Create Artists for all visible elements and add to our Axes."""
        x_lim, y_lim = self._axes_domain(nx=nx, ny=ny,
                                         background_patch=background_patch)

        transform = self._crs_transform()

        rc_params = matplotlib.rcParams

        n_steps = self.n_steps

        x_ticks = self.xlocator.tick_values(x_lim[0], x_lim[1])
        y_ticks = self.ylocator.tick_values(y_lim[0], y_lim[1])

        # XXX this bit is cartopy specific. (for circular longitudes)
        # Purpose: omit plotting the last x line, as it may overlap the first.
        x_gridline_points = x_ticks[:]
        crs = self.crs
        if (isinstance(crs, Projection) and
                isinstance(crs, _RectangularProjection) and
                abs(np.diff(x_lim)) == abs(np.diff(crs.x_limits))):
            x_gridline_points = x_gridline_points[:-1]

        collection_kwargs = self.collection_kwargs
        if collection_kwargs is None:
            collection_kwargs = {}
        collection_kwargs = collection_kwargs.copy()
        collection_kwargs['transform'] = transform
        # XXX doesn't gracefully handle lw vs linewidth aliases...
        collection_kwargs.setdefault('color', rc_params['grid.color'])
        collection_kwargs.setdefault('linestyle', rc_params['grid.linestyle'])
        collection_kwargs.setdefault('linewidth', rc_params['grid.linewidth'])

        if self.xlines:
            lines = []
            for x in x_gridline_points:
                l = list(zip(np.zeros(n_steps) + x,
                         np.linspace(min(y_ticks), max(y_ticks), n_steps)))
                lines.append(l)

            x_lc = mcollections.LineCollection(lines, **collection_kwargs)
            self.xline_artists.append(x_lc)
            self.axes.add_collection(x_lc, autolim=False)

        if self.ylines:
            lines = []
            for y in y_ticks:
                l = list(zip(np.linspace(min(x_ticks), max(x_ticks), n_steps),
                             np.zeros(n_steps) + y))
                lines.append(l)

            y_lc = mcollections.LineCollection(lines, **collection_kwargs)
            self.yline_artists.append(y_lc)
            self.axes.add_collection(y_lc, autolim=False)

        #################
        # Label drawing #
        #################

        # Trim outside-area points from the label coords.
        # Tickers may round *up* the desired range to something tidy, not
        # all of which is necessarily visible.  We must be stricter with
        # our texts, as they are drawn *without clipping*.
        x_label_points = [x for x in x_ticks if x_lim[0] <= x <= x_lim[1]]
        y_label_points = [y for y in y_ticks if y_lim[0] <= y <= y_lim[1]]

        if self.xlabels_bottom or self.xlabels_top:
            self._assert_can_draw_ticks()
            self.xformatter.set_locs(x_label_points)
            for x in x_label_points:
                if self.xlabels_bottom:
                    self._add_gridline_label(x, axis='x', upper_end=False)
                if self.xlabels_top:
                    self._add_gridline_label(x, axis='x', upper_end=True)

        if self.ylabels_left or self.ylabels_right:
            self._assert_can_draw_ticks()
            self.yformatter.set_locs(y_label_points)
            for y in y_label_points:
                if self.ylabels_left:
                    self._add_gridline_label(y, axis='y', upper_end=False)
                if self.ylabels_right:
                    self._add_gridline_label(y, axis='y', upper_end=True)

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
        if not isinstance(self.axes.projection,
                          (cartopy.crs.PlateCarree, cartopy.crs.Mercator)):
            raise TypeError('Cannot label gridlines on a '
                            '{prj.__class__.__name__} plot.  Only PlateCarree'
                            ' and Mercator plots are currently '
                            'supported.'.format(prj=self.axes.projection))
        return True

    def _axes_domain(self, nx=None, ny=None, background_patch=None):
        """Returns x_range, y_range"""
        DEBUG = False

        transform = self._crs_transform()

        ax_transform = self.axes.transAxes
        desired_trans = ax_transform - transform

        nx = nx or 30
        ny = ny or 30
        x = np.linspace(1e-9, 1 - 1e-9, nx)
        y = np.linspace(1e-9, 1 - 1e-9, ny)
        x, y = np.meshgrid(x, y)

        coords = np.concatenate([x.flatten()[:, None],
                                 y.flatten()[:, None]],
                                1)

        in_data = desired_trans.transform(coords)

        ax_to_bkg_patch = self.axes.transAxes - \
            background_patch.get_transform()

        ok = np.zeros(in_data.shape[:-1], dtype=np.bool)
        # XXX Vectorise contains_point
        for i, val in enumerate(in_data):
            # convert the coordinates of the data to the background
            # patches coordinates
            background_coord = ax_to_bkg_patch.transform(coords[i:i + 1, :])
            bkg_patch_contains = background_patch.get_path().contains_point
            if bkg_patch_contains(background_coord[0, :]):
                color = 'r'
                ok[i] = True
            else:
                color = 'b'

            if DEBUG:
                import matplotlib.pyplot as plt
                plt.plot(coords[i, 0], coords[i, 1], 'o' + color,
                         clip_on=False, transform=ax_transform)
#                plt.text(coords[i, 0], coords[i, 1], str(val), clip_on=False,
#                         transform=ax_transform, rotation=23,
#                         horizontalalignment='right')

        inside = in_data[ok, :]

        # If there were no data points in the axes we just use the x and y
        # range of the projection.
        if inside.size == 0:
            x_range = self.crs.x_limits
            y_range = self.crs.y_limits
        else:
            x_range = np.nanmin(inside[:, 0]), np.nanmax(inside[:, 0])
            y_range = np.nanmin(inside[:, 1]), np.nanmax(inside[:, 1])

        # XXX Cartopy specific thing. Perhaps make this bit a specialisation
        # in a subclass...
        crs = self.crs
        if isinstance(crs, Projection):
            x_range = np.clip(x_range, *crs.x_limits)
            y_range = np.clip(y_range, *crs.y_limits)

            # if the limit is >90% of the full x limit, then just use the full
            # x limit (this makes circular handling better)
            prct = np.abs(np.diff(x_range) / np.diff(crs.x_limits))
            if prct > 0.9:
                x_range = crs.x_limits

        return x_range, y_range
