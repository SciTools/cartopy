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
                 ylocator=None, collection_kwargs=None,
                 xformatter=None, yformatter=None):
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
            longitude labels.
        yformatter: optional
            A :class:`matplotlib.ticker.Formatter` instance to format
            latitude labels.
        collection_kwargs: optional
            Dictionary controlling line properties, passed to
            :class:`matplotlib.collections.Collection`. Defaults to None.

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
            self.xlocator = LongitudeLocator()
        else:
            self.xlocator = classic_locator

        #: The :class:`~matplotlib.ticker.Locator` to use for the y
        #: gridlines and labels.
        if ylocator is not None:
            if not isinstance(ylocator, mticker.Locator):
                ylocator = mticker.FixedLocator(ylocator)
            self.ylocator = ylocator
        elif isinstance(crs, cartopy.crs.PlateCarree):
            self.ylocator = LatitudeLocator()
        else:
            self.ylocator = classic_locator

        #: The :class:`~matplotlib.ticker.Formatter` to use for the lon labels.
        self.xformatter = xformatter or LongitudeFormatter()

        #: The :class:`~matplotlib.ticker.Formatter` to use for the lat labels.
        self.yformatter = yformatter or LatitudeFormatter()

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

    def _add_gridline_label(self, loc, text, edge):
        """
        Create a Text artist on our axes for a gridline label.

        Parameters
        ----------
        loc
            Coordinate value of this gridline in projected units
            (not in degrees).
        text
            Text to display, positioned centred at loc.
        edge: str
            Edge of the plot on wich to draw the label.
            One of 'top', 'bottom', 'right' or 'left'.
        """
        # Get label specs
        if edge in ['top', 'bottom']:
            axis = 'x'
            x = loc
            if edge == 'top':
                y = 1.0
                meth = self.axes.get_xaxis_text2_transform
            else:
                y = 0.0
                meth = self.axes.get_xaxis_text1_transform
            label_transform, v_align, h_align = meth(self.xpadding)
            user_label_style = self.xlabel_style
        elif edge in ['right', 'left']:
            axis = 'y'
            y = loc
            if edge == 'right':
                x = 1.0
                meth = self.axes.get_yaxis_text2_transform
            else:
                x = 0.0
                meth = self.axes.get_yaxis_text1_transform
            label_transform, v_align, h_align = meth(self.ypadding)
            if matplotlib.__version__ > '2.0':
                v_align = 'center_baseline'
            else:
                v_align = 'center'
            user_label_style = self.ylabel_style
        else:
            raise ValueError(
                "Unknown edge, {!r}, must be either 'top', 'bottom', "
                "'left' or 'right'".format(axis))

        # Create and add a Text artist with these properties
        label_style = {'verticalalignment': v_align,
                       'horizontalalignment': h_align,
                       }
        label_style.update(user_label_style)
        text_artist = mtext.Text(x, y, text,
                                 clip_on=False,
                                 transform=label_transform, **label_style)

        # Check that this artist don't overlap this one
#        artists = getattr(self, '{}label_{}_artists'.format(axis, edge))
        artists = getattr(self, edge + '_label_artists')
        for ta in artists:
            if _text_artists_overlaps_(ta, text_artist, axis):
                return

        # Ok, register it
        artists.append(text_artist)
        self.axes.add_artist(text_artist)

    def _draw_gridliner(self, nx=None, ny=None, background_patch=None):
        """Create Artists for all visible elements and add to our Axes."""
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
                     if value >= crs.x_limits[0] and value <= crs.x_limits[1]]
        lat_ticks = [value for value in lat_ticks
                     if value >= crs.y_limits[0] and value <= crs.y_limits[1]]

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

        # Longitude lines
        lon_lines = []
        for x in lon_ticks:
            ticks = list(zip(
                np.zeros(n_steps) + x,
                np.linspace(min(lat_lim[0], lat_ticks[0]),
                            max(lat_lim[1], lat_ticks[-1]), n_steps)))
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

        # Latitude lines
        lat_lines = []
        for y in lat_ticks:
            ticks = list(zip(
                np.linspace(min(lon_lim[0], lon_ticks[0]),
                            max(lon_lim[1], lon_ticks[-1]), n_steps),
                np.zeros(n_steps) + y))
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
                ('lon', lon_lines, lon_ticks, self.xformatter),
                ('lat', lat_lines, lat_ticks, self.yformatter)):

            to_draw[lonlat] = {'y': {'left': [], 'right': []},
                               'x': {'bottom': [], 'top': []}}
            formatter.set_locs(line_ticks)

            for line, tick_value in zip(lines, line_ticks):

                ls = sgeom.LineString(line)
                lsp = self.axes.projection.project_geometry(ls)

                if lsp.intersects(map_boundaries):

                    lsb = lsp.intersection(map_boundaries)
                    if lsb.is_empty:
                        continue

                    lsb = lsb.boundary
                    for point in lsb:
                        x = round(point.x, 5)
                        y = round(point.y, 5)
                        text = formatter(tick_value)
                        checks = [(x, round(x0, 5), y, 'y', 'left'),
                                  (x, round(x1, 5), y, 'y', 'right'),
                                  (y, round(y0, 5), x, 'x', 'bottom'),
                                  (y, round(y1, 5), x, 'x', 'top')]
                        if lonlat == 'lon':
                            checks = checks[2:] + checks[:2]
                        for xy, xy01, loc, axis, edge in checks:
                            if xy == xy01:
                                break
                        else:
                            continue

                        if getattr(self, edge + '_labels'):
                            to_draw[lonlat][axis][edge].append((loc, text))

        # Draw in the prefered order
        for lonlat, axis in (('lon', 'x'), ('lat', 'y'),
                             ('lon', 'y'), ('lat', 'x')):
            for edge, lts in to_draw[lonlat][axis].items():
                for loc, text in lts:
                    self._add_gridline_label(loc, text, edge)

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
            lon_range = self.crs.x_limits
            lat_range = self.crs.y_limits
        else:
            lon_range = np.nanmin(inside[:, 0]), np.nanmax(inside[:, 0])
            lat_range = np.nanmin(inside[:, 1]), np.nanmax(inside[:, 1])

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
