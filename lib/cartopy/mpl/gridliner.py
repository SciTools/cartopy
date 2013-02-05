# (C) British Crown Copyright 2011 - 2012, Met Office
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
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib
import matplotlib.collections as mcollections
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.transforms as mtrans
import numpy

import cartopy
from cartopy.crs import Projection, _RectangularProjection


degree_locator = mticker.MaxNLocator(nbins=9, steps=[1, 2, 3, 6, 15, 18])


class Gridliner(object):
    # NOTE: In future, one of these objects will be add-able to a GeoAxes (and
    # maybe even a plain old mpl axes) and it will call the "do_gridlines"
    # method on draw. This will enable automatic gridline resolution
    # determination on zoom/pan.
    def __init__(self, axes, crs, draw_labels=False, collection_kwargs=None):
        """
        Object used by :function:`cartopy.mpl.geoaxes.Geoaxes.gridlines`
        to add gridlines to a map.

        Args:

        * axes
            The :class:`cartopy.mpl.geoaxes.GeoAxes` object to be drawn on.

        * crs
            The :class:`cartopy.crs.CRS` defining the coordinate system that
            the gridlines are drawn in.

        * draw_labels
            Label the gridlines at the map edges.

        * collection_kwargs
            Dictionary controlling line properties, passed to
            :class:`matplotlib.collections.Collection`.

        """
        self.axes = axes
        # might not be desirable for certain projections (osgb for instance)
        self.xlocator = degree_locator
        self.ylocator = degree_locator
        self.crs = crs
        if draw_labels:
            # Check labelling is supported, currently a limited set of options.
            if not isinstance(crs, cartopy.crs.PlateCarree):
                raise TypeError(
                    'Cannot label {} gridlines.  Only PlateCarree gridlines '
                    'are currently supported.'.format(crs.__class__.__name__))
            if not isinstance(axes.projection,
                              (cartopy.crs.PlateCarree, cartopy.crs.Mercator)):
                raise TypeError(
                    'Cannot label gridlines on a {} plot.  Only PlateCarree '
                    'and Mercator plots are currently supported.'.format(
                        axes.projection.__class__.__name__))
        self.add_labels = draw_labels
        self.collection_kwargs = collection_kwargs

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
        else:
            raise ValueError(
                "Unknown axis, {!r}, must be either 'x' or 'y'".format(axis))

        # Make a 'blended' transform for label text positioning.
        # One coord is geographic, and the other a plain Axes
        # coordinate with an appropriate offset.
        label_transform = mtrans.blended_transform_factory(
            x_transform=tr_x, y_transform=tr_y)

        # Create and add a Text artist with these properties
        text_artist = mtext.Text(
            x, y, '{:g}'.format(value),
            clip_on=False,
            verticalalignment=v_align,
            horizontalalignment=h_align,
            transform=label_transform)
        self.axes.add_artist(text_artist)

    def do_gridlines(self, nx=None, ny=None, background_patch=None):
        """Create Artists for all visible elements and add to our Axes."""
        x_lim, y_lim = self.get_domain(nx=nx, ny=ny,
                                       background_patch=background_patch)

        transform = self._crs_transform()

        rc_params = matplotlib.rcParams

        n_steps = 30

        lines = []

        x_ticks = self.xlocator.tick_values(x_lim[0], x_lim[1])
        y_ticks = self.ylocator.tick_values(y_lim[0], y_lim[1])

        # XXX this bit is cartopy specific. (for circular longitudes)
        # Purpose: omit plotting the last x line, as it may overlap the first.
        x_gridline_points = x_ticks[:]
        crs = self.crs
        if (isinstance(crs, Projection) and
                isinstance(crs, _RectangularProjection) and
                numpy.diff(x_lim) == 2 * crs._half_width):
            x_gridline_points = x_gridline_points[:-1]

        for x in x_gridline_points:
            l = zip(numpy.zeros(n_steps) + x,
                    numpy.linspace(min(y_ticks), max(y_ticks),
                                   n_steps, endpoint=True)
                    )
            lines.append(l)

        collection_kwargs = self.collection_kwargs
        if collection_kwargs is None:
            collection_kwargs = {}
        collection_kwargs = collection_kwargs.copy()
        collection_kwargs['transform'] = transform
        # XXX doesn't gracefully handle lw vs linewidth aliases...
        collection_kwargs.setdefault('color', rc_params['grid.color'])
        collection_kwargs.setdefault('linestyle', rc_params['grid.linestyle'])
        collection_kwargs.setdefault('linewidth', rc_params['grid.linewidth'])

        x_lc = mcollections.LineCollection(lines, **collection_kwargs)
        self.axes.add_collection(x_lc, autolim=False)

        lines = []

        for y in y_ticks:
            l = zip(numpy.linspace(min(x_ticks), max(x_ticks), n_steps),
                    numpy.zeros(n_steps) + y)
            lines.append(l)

        y_lc = mcollections.LineCollection(lines, **collection_kwargs)
        self.axes.add_collection(y_lc, autolim=False)

        if self.add_labels:
            # Trim outside-area points from the label coords.
            # Tickers may round *up* the desired range to something tidy, not
            # all of which is necessarily visible.  We must be stricter with
            # our texts, as they are drawn *without clipping*.
            x_label_points = [x for x in x_ticks if x_lim[0] <= x <= x_lim[1]]
            y_label_points = [y for y in y_ticks if y_lim[0] <= y <= y_lim[1]]

            # Add text labels at each end of every gridline.
            for upper_end in (False, True):
                for x in x_label_points:
                    self._add_gridline_label(x, axis='x', upper_end=upper_end)
                for y in y_label_points:
                    self._add_gridline_label(y, axis='y', upper_end=upper_end)

    def get_domain(self, nx=None, ny=None, background_patch=None):
        """Returns x_range, y_range"""
        DEBUG = False

        transform = self._crs_transform()

        ax_transform = self.axes.transAxes
        desired_trans = ax_transform - transform

        nx = nx or 30
        ny = ny or 30
        x = numpy.linspace(1e-9, 1 - 1e-9, nx)
        y = numpy.linspace(1e-9, 1 - 1e-9, ny)
        x, y = numpy.meshgrid(x, y)

        coords = numpy.concatenate([x.flatten()[:, None],
                                    y.flatten()[:, None]],
                                   1)

        in_data = desired_trans.transform(coords)

        ax_to_bkg_patch = self.axes.transAxes - \
            background_patch.get_transform()

        ok = numpy.zeros(in_data.shape[:-1], dtype=numpy.bool)
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
        x_range = numpy.nanmin(inside[:, 0]), numpy.nanmax(inside[:, 0])
        y_range = numpy.nanmin(inside[:, 1]), numpy.nanmax(inside[:, 1])

        # XXX Cartopy specific thing. Perhaps make this bit a specialisation
        # in a subclass...
        crs = self.crs
        if isinstance(crs, Projection):
            x_range = numpy.clip(x_range, *crs.x_limits)
            y_range = numpy.clip(y_range, *crs.y_limits)

            # if the limit is >90 of the full x limit, then just use the full
            # x limit (this makes circular handling better)
            prct = numpy.abs(numpy.diff(x_range) / numpy.diff(crs.x_limits))
            if prct > 0.9:
                x_range = crs.x_limits

        return x_range, y_range
