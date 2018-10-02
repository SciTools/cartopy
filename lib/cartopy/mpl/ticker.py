# (C) British Crown Copyright 2014 - 2018, Met Office
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
"""This module contains tools for handling tick marks in cartopy."""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from matplotlib.ticker import Formatter, MaxNLocator

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes


class _PlateCarreeFormatter(Formatter):
    """
    Base class for formatting ticks on geographical axes using a
    rectangular projection (e.g. Plate Carree, Mercator).

    """

    _target_projection = ccrs.PlateCarree()

    def __init__(self, degree_symbol=u'\u00B0', decimal=False,
                 minute_symbol=u"'", second_symbol=u"''",
                 number_format='g', seconds_number_format='g',
                 auto_hide=True, transform_precision=1e-8):
        """
        Base class for simpler implementation of specialised formatters
        for latitude and longitude axes.

        """
        self._degree_symbol = degree_symbol
        self._degrees_number_format = number_format
        self._transform_precision = transform_precision
        self._decimal = decimal
        self._minute_symbol = minute_symbol
        self._second_symbol = second_symbol
        self._seconds_num_format = seconds_number_format
        self._auto_hide = auto_hide
        self._auto_hide_degrees = False
        self._auto_hide_minutes = False
        self._precision = 5  # locator precision

    def __call__(self, value, pos=None):
        if self.axis is not None:

            if not isinstance(self.axis.axes, GeoAxes):
                raise TypeError("This formatter can only be "
                                "used with cartopy axes.")

            # We want to produce labels for values in the familiar Plate Carree
            # projection, so convert the tick values from their own projection
            # before formatting them.
            source = self.axis.axes.projection
            if not isinstance(source, (ccrs._RectangularProjection,
                                       ccrs.Mercator)):
                raise TypeError("This formatter cannot be used with "
                                "non-rectangular projections.")
            projected_value = self._apply_transform(value,
                                                    self._target_projection,
                                                    source)

            # Round the transformed value using a given precision for display
            # purposes. Transforms can introduce minor rounding errors that
            # make the tick values look bad, these need to be accounted for.
            f = 1. / self._transform_precision
            projected_value = round(f * projected_value) / f

        else:

            # There is no projection so we assume it is already PlateCarree
            projected_value = value

        # Return the formatted values, the formatter has both the re-projected
        # tick value and the original axis value available to it.
        return self._format_value(projected_value, value)

    def _format_value(self, value, original_value):
        hemisphere = self._hemisphere(value, original_value)

        if self._decimal:
            return (self._format_degrees(abs(value)) +
                    hemisphere)

        value, deg, mn, sec = self._get_dms(abs(value))

        # Format
        label = u''
        if sec:
            label = self._format_seconds(sec)

        if mn or (not self._auto_hide_minutes and label):
            label = self._format_minutes(mn) + label

        if not self._auto_hide_degrees or not label:
            label = self._format_degrees(deg) + hemisphere + label

        return label

    def _get_dms(self, x):
        x = np.asarray(x, 'd')
        degs = np.round(x, self._precision).astype('i')
        x = (x - degs) * 60
        mins = np.round(x, self._precision).astype('i')
        secs = np.round((x - mins) * 60, self._precision - 3)
        return x, degs, mins, secs

    def set_locs(self, locs):

        Formatter.set_locs(self, locs)
        if not self._auto_hide:
            return

        self.locs, degs, mins, secs = self._get_dms(self.locs)
        secs = np.round(secs, self._precision-3).astype('i')
        secs0 = secs == 0
        mins0 = mins == 0

        def auto_hide(valid, values):
            """Should I switch on auto_hide?"""
            if not valid.any():
                return False
            if valid.sum() == 1:
                return True
            return np.diff(degs.compress(valid)).max() == 1

        # Potentially hide minutes labels when pure minutes are all displayed
        self._auto_hide_minutes = auto_hide(secs0, mins)

        # Potentially hide degrees labels when pure degrees are all displayed
        self._auto_hide_degrees = auto_hide(secs0 & mins0, degs)

    def _format_degrees(self, deg):
        """Format degrees as an integer"""
        if not self._decimal:
            deg = int(deg)
            number_format = 'd'
        else:
            number_format = self._degrees_number_format
        return u'{value:{number_format}}{symbol}'.format(
                value=abs(deg),
                number_format=number_format,
                symbol=self._degree_symbol)

    def _format_minutes(self, mn):
        """Format minutes as an integer"""
        return u'{value:d}{symbol}'.format(
                value=int(mn),
                symbol=self._minute_symbol)

    def _format_seconds(self, sec):
        """Format seconds as an float"""
        return u'{value:{fmt}}{symbol}'.format(
                value=sec,
                fmt=self._seconds_num_format,
                symbol=self._second_symbol)

    def _apply_transform(self, value, target_proj, source_crs):
        """
        Given a single value, a target projection and a source CRS,
        transform the value from the source CRS to the target
        projection, returning a single value.

        """
        raise NotImplementedError("A subclass must implement this method.")

    def _hemisphere(self, value, value_source_crs):
        """
        Given both a tick value in the Plate Carree projection and the
        same value in the source CRS, return a string indicating the
        hemisphere that the value is in.

        Must be over-ridden by the derived class.

        """
        raise NotImplementedError("A subclass must implement this method.")


class LatitudeFormatter(_PlateCarreeFormatter):
    """Tick formatter for latitude axes."""
    def __init__(self,
                 degree_symbol=u'\u00B0',
                 decimal=False,
                 minute_symbol=u"'",
                 second_symbol=u"''",
                 number_format='g',
                 seconds_number_format='g',
                 auto_hide=True,
                 transform_precision=1e-8):
        """
        Tick formatter for a latitudes.

        When bounded to an axis, the axis must be part of an axes defined
        on a rectangular projection (e.g. Plate Carree, Mercator).


        Parameters
        ----------
        degree_symbol: optional
            The character(s) used to represent the degree symbol in the
            tick labels. Defaults to u'\u00B0' which is the unicode
            degree symbol. Can be an empty string if no degree symbol is
            desired.
        decimal: bool, optional
            Wether or not formatting as decimal degrees and not as
            degrees-minutes-seconds.
        minute_symbol: str, optional
            The character(s) used to represent the minute symbol.
        second_symbol: str, optional
            The character(s) used to represent the second symbol.
        auto_hide: bool, optional
            Auto-hide degrees or minutes when redondant.
        number_format: optional
            Format string to represent the degrees tick values.
            Defaults to 'g'.
        seconds_number_format: optional
            Format string to represent the seconds tick values.
            Defaults to 'g'.
        transform_precision: optional
            Sets the precision (in degrees) to which transformed tick
            values are rounded. The default is 1e-7, and should be
            suitable for most use cases. To control the appearance of
            tick labels use the *number_format* keyword.

        Note
        ----
            A formatter can only be used for one axis. A new formatter
            must be created for every axis that needs formatted labels.

        Examples
        --------
        Label latitudes from -90 to 90 on a Plate Carree projection::

            ax = plt.axes(projection=PlateCarree())
            ax.set_global()
            ax.set_yticks([-90, -60, -30, 0, 30, 60, 90],
                          crs=ccrs.PlateCarree())
            lat_formatter = LatitudeFormatter()
            ax.yaxis.set_major_formatter(lat_formatter)

        Label latitudes from -80 to 80 on a Mercator projection, this
        time omitting the degree symbol::

            ax = plt.axes(projection=Mercator())
            ax.set_global()
            ax.set_yticks([-90, -60, -30, 0, 30, 60, 90],
                          crs=ccrs.PlateCarree())
            lat_formatter = LatitudeFormatter(degree_symbol='')
            ax.yaxis.set_major_formatter(lat_formatter)

        When not bounded to an axis::

            lat_formatter = LatitudeFormatter()
            ticks = [-90, -60, -30, 0, 30, 60, 90]
            lat_formatter.set_locs(ticks)
            labels = [lat_formatter(value) for value in ticks]

        """
        super(LatitudeFormatter, self).__init__(
            degree_symbol=degree_symbol,
            number_format=number_format,
            transform_precision=transform_precision,
            decimal=decimal,
            minute_symbol=minute_symbol,
            second_symbol=second_symbol,
            seconds_number_format=seconds_number_format,
            auto_hide=auto_hide,
            )

    def _apply_transform(self, value, target_proj, source_crs):
        return target_proj.transform_point(0, value, source_crs)[1]

    def _hemisphere(self, value, value_source_crs):
        if value > 0:
            hemisphere = 'N'
        elif value < 0:
            hemisphere = 'S'
        else:
            hemisphere = ''
        return hemisphere


class LongitudeFormatter(_PlateCarreeFormatter):
    """Tick formatter for a longitude axis."""

    def __init__(self,
                 zero_direction_label=False,
                 dateline_direction_label=False,
                 degree_symbol=u'\u00B0',
                 decimal=False,
                 minute_symbol=u"'",
                 second_symbol=u"''",
                 number_format='g',
                 seconds_number_format='g',
                 auto_hide=True,
                 transform_precision=1e-8):
        """
        Create a formatter for longitudes.

        When bounded to an axis, the axis must be part of an axes defined
        on a rectangular projection (e.g. Plate Carree, Mercator).

        Parameters
        ----------
        zero_direction_label: optional
            If *True* a direction label (E or W) will be drawn next to
            longitude labels with the value 0. If *False* then these
            labels will not be drawn. Defaults to *False* (no direction
            labels).
        dateline_direction_label: optional
            If *True* a direction label (E or W) will be drawn next to
            longitude labels with the value 180. If *False* then these
            labels will not be drawn. Defaults to *False* (no direction
            labels).
        degree_symbol: optional
            The symbol used to represent degrees. Defaults to u'\u00B0'
            which is the unicode degree symbol.
        decimal: bool, optional
            Wether or not formatting as decimal degrees and not as
            degrees-minutes-seconds.
        auto_hide: bool, optional
            Auto-hide degrees or minutes when redondant.
        minute_symbol: str, optional
            The character(s) used to represent the minute symbol.
        second_symbol: str, optional
            The character(s) used to represent the second symbol.
        number_format: optional
            Format string to represent the degrees tick values.
            Defaults to 'g'.
        seconds_number_format: optional
            Format string to represent the seconds tick values.
            Defaults to 'g'.
        transform_precision: optional
            Sets the precision (in degrees) to which transformed tick
            values are rounded. The default is 1e-7, and should be
            suitable for most use cases. To control the appearance of
            tick labels use the *number_format* keyword.

        Note
        ----
            A formatter can only be used for one axis. A new formatter
            must be created for every axis that needs formatted labels.

        Examples
        --------
        Label longitudes from -180 to 180 on a Plate Carree projection
        with a central longitude of 0::

            ax = plt.axes(projection=PlateCarree())
            ax.set_global()
            ax.set_xticks([-180, -120, -60, 0, 60, 120, 180],
                          crs=ccrs.PlateCarree())
            lon_formatter = LongitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)

        Label longitudes from 0 to 360 on a Plate Carree projection
        with a central longitude of 180::

            ax = plt.axes(projection=PlateCarree(central_longitude=180))
            ax.set_global()
            ax.set_xticks([0, 60, 120, 180, 240, 300, 360],
                          crs=ccrs.PlateCarree())
            lon_formatter = LongitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)


        When not bounded to an axis::

            lon_formatter = LongitudeFormatter()
            ticks = [0, 60, 120, 180, 240, 300, 360]
            lon_formatter.set_locs(ticks)
            labels = [lon_formatter(value) for value in ticks]
        """
        super(LongitudeFormatter, self).__init__(
            degree_symbol=degree_symbol,
            number_format=number_format,
            transform_precision=transform_precision,
            decimal=decimal,
            minute_symbol=minute_symbol,
            second_symbol=second_symbol,
            seconds_number_format=seconds_number_format,
            auto_hide=auto_hide,
            )
        self._zero_direction_labels = zero_direction_label
        self._dateline_direction_labels = dateline_direction_label

    def _apply_transform(self, value, target_proj, source_crs):
        return target_proj.transform_point(value, 0, source_crs)[0]

    @classmethod
    def _fix_lons(cls, lons):
        if isinstance(lons, list):
            return map(cls._fix_lons, lons)
        lons = ((lons + 180) % 360) - 180
        if isinstance(lons, np.ndarray):
            valid = lons == -180
            if valid.any():
                lons[valid] = 180
        elif lons == -180:
            lons = 180
        return lons

    def set_locs(self, locs):
        _PlateCarreeFormatter.set_locs(self, self._fix_lons(locs))

    def _format_degrees(self, deg):
        return _PlateCarreeFormatter._format_degrees(self, self._fix_lons(deg))

    def _hemisphere(self, value, value_source_crs):
        value = self._fix_lons(value)
        # Perform basic hemisphere detection.
        if value < 0:
            hemisphere = 'W'
        elif value > 0:
            hemisphere = 'E'
        else:
            hemisphere = ''
        # Correct for user preferences:
        if value == 0 and self._zero_direction_labels:
            # Use the original tick value to determine the hemisphere.
            if value_source_crs < 0:
                hemisphere = 'E'
            else:
                hemisphere = 'W'
        if value in (-180, 180) and not self._dateline_direction_labels:
            hemisphere = ''
        return hemisphere


class LongitudeLocator(MaxNLocator):
    """A locator for longitude degrees that works even at very small scale
    on minutes and seconds

    Parameters
    ----------
    decimal: bool
        Force decimal degrees so the local does not stop specifically on
        fractions of minutes and seconds (False by default)
    """

    default_params = MaxNLocator.default_params.copy()
    default_params.update(nbins=8, decimal=False)

    def set_params(self, **kwargs):
        """Set parameters within this locator."""
        MaxNLocator.set_params(self, **kwargs)
        if 'decimal' in kwargs:
            self._decimal = kwargs['decimal']

    def _guess_steps(self, vmin, vmax):

        dv = abs(vmax - vmin)
        if dv > 180:
            dv -= 180

        if dv > 50.:

            steps = np.array([1, 2, 3, 6, 10])

        elif self._decimal or dv > 3.:

            steps = np.array([1, 1.5, 2, 2.5, 3, 5, 10])

        else:
            steps = np.array([1, 10/6., 15/6., 20/6., 30/6., 10])

        self.set_params(steps=np.array(steps))

    def _raw_ticks(self, vmin, vmax):
        self._guess_steps(vmin, vmax)
        return MaxNLocator._raw_ticks(self, vmin, vmax)


class LatitudeLocator(LongitudeLocator):
    """A locator for latitude degrees that works even at very small scale
    on minutes and seconds

    Parameters
    ----------
    decimal: bool
        Force decimal degrees so the local does not stop specifically on
        fractions of minutes and seconds (False by default)
    """

    def _raw_ticks(self, vmin, vmax):
        vmin = max(vmin, -90.)
        vmax = min(vmax, 90.)
        ticks = LongitudeLocator._raw_ticks(self, vmin, vmax)
        return [t for t in ticks if t >= -90 and t <= 90]
