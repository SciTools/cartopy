# (C) British Crown Copyright 2011 - 2018, Met Office
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
import shapely.geometry as sgeom

import cartopy
from cartopy.crs import Projection, _RectangularProjection


degree_locator = mticker.MaxNLocator(nbins=9, steps=[1, 1.5, 1.8, 2, 3, 6, 10])

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


def _text_artists_overlaps_(artista, artistb, axis):
    """Check that two text artists don't overlap (approximately)"""
    xa, ya = artista.get_position()
    xb, yb = artistb.get_position()
    trans = artista.get_transform()
    xa, ya = trans.transform_point((xa, ya))
    xb, yb = trans.transform_point((xb, yb))
    size = artista.get_size()
    if axis == 'x':
        factor = 0.8 if 'monospace' not in artista.get_family() else 1
        dxa = size * len(artista.get_text()) / 2 * factor
        dxb = size * len(artistb.get_text()) / 2 * factor
        return not ((xa + dxa) < (xb - dxb) or (xb + dxb) < (xa - dxa))
    dya = size / 2
    dyb = size / 2
    return not ((ya + dya) < (yb - dyb) or (yb + dyb) < (ya - dya))


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
        collection_kwargs: optional
            Dictionary controlling line properties, passed to
            :class:`matplotlib.collections.Collection`. Defaults to None.

        """
        self.axes = axes

        #: The :class:`~matplotlib.ticker.Locator` to use for the x
        #: gridlines and labels.
        self.xlocator = xlocator or degree_locator

        #: The :class:`~matplotlib.ticker.Locator` to use for the y
        #: gridlines and labels.
        self.ylocator = ylocator or degree_locator

        #: The :class:`~matplotlib.ticker.Formatter` to use for the x labels.
        self.xformatter = LONGITUDE_FORMATTER

        #: The :class:`~matplotlib.ticker.Formatter` to use for the y labels.
        self.yformatter = LATITUDE_FORMATTER

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

        #: The padding from the map edge to the x labels in points.
        self.xpadding = 5

        #: The padding from the map edge to the y labels in points.
        self.ypadding = 5

        self.crs = crs

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

        Note
        ----
            The drawing transform depends on the transform of our 'axes', so
            it may change dynamically.

        """
        transform = self.crs
        if not isinstance(transform, mtrans.Transform):
            transform = transform._as_mpl_transform(self.axes)
        return transform

    def _add_gridline_label(self, loc, text, axis, edge):
        """
        Create a Text artist on our axes for a gridline label.

        Parameters
        ----------
        value
            Coordinate value of this gridline.  The text contains this
            value, and is positioned centred at that point.
        axis
            Which axis the label is on: 'x' or 'y'.
        edge: str
            If "top" or "right", place at the maximum of the "other"
            coordinate (Axes coordinate == 1.0).  Else 'lower' end
            (Axes coord = 0.0).

        """
        # Get label specs
        if axis == 'x':
            x = loc
            if edge == 'top':
                y = 1.0
                meth = self.axes.get_xaxis_text2_transform
            else:
                y = 0.0
                meth = self.axes.get_xaxis_text1_transform
            label_transform, v_align, h_align = meth(self.xpadding)
            user_label_style = self.xlabel_style
        elif axis == 'y':
            y = loc
            if edge == 'right':
                x = 1.0
                meth = self.axes.get_yaxis_text2_transform
            else:
                x = 0.0
                meth = self.axes.get_yaxis_text1_transform
            label_transform, v_align, h_align = meth(self.xpadding)
            if matplotlib.__version__ > '2.0':
                v_align = 'center_baseline'
            else:
                v_align = 'center'
            user_label_style = self.ylabel_style
        else:
            raise ValueError(
                "Unknown axis, {!r}, must be either 'x' or 'y'".format(axis))

        # Create and add a Text artist with these properties
        label_style = {'verticalalignment': v_align,
                       'horizontalalignment': h_align,
                       }
        label_style.update(user_label_style)
        text_artist = mtext.Text(x, y, text,
                                 clip_on=False,
                                 transform=label_transform, **label_style)

        # Check that this artist don't overlap this one
        artists = getattr(self, '{}label_{}_artists'.format(axis, edge))
        for ta in artists:
            if _text_artists_overlaps_(ta, text_artist, axis):
                return

        # Ok, register it
        artists.append(text_artist)
        self.axes.add_artist(text_artist)
        return text_artist

    def _draw_gridliner(self, nx=None, ny=None, background_patch=None):
        """Create Artists for all visible elements and add to our Axes."""
        x_lim, y_lim = self._axes_domain(nx=nx, ny=ny,
                                         background_patch=background_patch)

        transform = self._crs_transform()

        rc_params = matplotlib.rcParams

        n_steps = self.n_steps

        x_ticks = self.xlocator.tick_values(x_lim[0], x_lim[1])
        y_ticks = self.ylocator.tick_values(y_lim[0], y_lim[1])

        #####################
        # Gridlines drawing #
        #####################

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

        # Longitude lines
        lon_lines = []
        for x in x_gridline_points:
            ticks = list(zip(
                np.zeros(n_steps) + x,
                np.linspace(min(y_ticks), max(y_ticks), n_steps)))
            lon_lines.append(ticks)

        if self.xlines:
            x_lc = mcollections.LineCollection(lon_lines, **collection_kwargs)
            self.xline_artists.append(x_lc)
            self.axes.add_collection(x_lc, autolim=False)

        # Latitude lines
        lat_lines = []
        for y in y_ticks:
            ticks = list(zip(
                np.linspace(min(x_ticks), max(x_ticks), n_steps),
                np.zeros(n_steps) + y))
            lat_lines.append(ticks)
        if self.ylines:
            y_lc = mcollections.LineCollection(lat_lines, **collection_kwargs)
            self.yline_artists.append(y_lc)
            self.axes.add_collection(y_lc, autolim=False)

        #################
        # Label drawing #
        #################

        self.xlabel_bottom_artists = []
        self.xlabel_top_artists = []
        self.ylabel_left_artists = []
        self.ylabel_right_artists = []
        if not (self.ylabels_left or self.ylabels_right or
                self.xlabels_bottom or self.xlabels_top):
            return

        # Get the real map boundaries
        x0, x1 = self.axes.get_xlim()
        y0, y1 = self.axes.get_ylim()
        plot_boundary = sgeom.Polygon([[x0, y0], [x1, y0],
                                       [x1, y1], [x0, y1]])
        map_boundary_vertices = self.axes.outline_patch.get_path().vertices
        map_boundary = sgeom.Polygon(map_boundary_vertices)
        map_boundaries = plot_boundary.intersection(map_boundary)

        # Loop on longitude and latitude lines and collect what to draw
        to_draw = {}
        for lonlat, lines, line_ticks, formatter in (
                ('lon', lon_lines, x_ticks, self.xformatter),
                ('lat', lat_lines, y_ticks, self.yformatter)):

            to_draw[lonlat] = {'y': {'left': [], 'right': []},
                               'x': {'bottom': [], 'top': []}}

            for line, tick_value in zip(lines, line_ticks):

                ls = sgeom.LineString(line)
                lsp = self.axes.projection.project_geometry(ls)

                if lsp.intersects(map_boundaries):

                    lsb = lsp.intersection(map_boundaries)
                    if lsb.is_empty:
                        continue

                    lsb = lsb.boundary
                    for point in lsb:
                        x = point.x
                        y = point.y
                        text = formatter(tick_value)
                        if x == x0:
                            loc, axis, edge = y, 'y', 'left'
                        elif x == x1:
                            loc, axis, edge = y, 'y', 'right'
                        elif y == y0:
                            loc, axis, edge = x, 'x', 'bottom'
                        elif y == y1:
                            loc, axis, edge = x, 'x', 'top'
                        else:
                            continue

                        if getattr(self, '{}labels_{}'.format(axis, edge)):
                            to_draw[lonlat][axis][edge].append((loc, text))

        # Draw in the prefered order
        for lonlat, axis in (('lon', 'x'), ('lat', 'y'),
                             ('lon', 'y'), ('lat', 'x')):
            for edge, lts in to_draw[lonlat][axis].items():
                for loc, text in lts:
                    self._add_gridline_label(loc, text, axis, edge)

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
        """Return x_range, y_range"""
        DEBUG = False

        transform = self._crs_transform()

        ax_transform = self.axes.transAxes
        desired_trans = ax_transform - transform

        nx = nx or 30
        ny = ny or 30
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
