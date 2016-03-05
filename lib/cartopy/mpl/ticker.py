# (C) British Crown Copyright 2014 - 2016, Met Office
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

from matplotlib.ticker import Formatter

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes


class _PlateCarreeFormatter(Formatter):
    """
    Base class for formatting ticks on geographical axes using a
    rectangular projection (e.g. Plate Carree, Mercator).

    """

    _target_projection = ccrs.PlateCarree()

    def __init__(self, degree_symbol=u'\u00B0', number_format='g',
                 transform_precision=1e-8):
        """
        Base class for simpler implementation of specialised formatters
        for latitude and longitude axes.

        """
        self._degree_symbol = degree_symbol
        self._number_format = number_format
        self._transform_precision = transform_precision

    def __call__(self, value, pos=None):
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
        projected_value = self._apply_transform(value, self._target_projection,
                                                source)
        # Round the transformed value using a given precision for display
        # purposes. Transforms can introduce minor rounding errors that make
        # the tick values look bad, these need to be accounted for.
        f = 1. / self._transform_precision
        projected_value = round(f * projected_value) / f
        # Return the formatted values, the formatter has both the re-projected
        # tick value and the original axis value available to it.
        return self._format_value(projected_value, value)

    def _format_value(self, value, original_value):
        hemisphere = self._hemisphere(value, original_value)
        fmt_string = u'{value:{number_format}}{degree}{hemisphere}'
        return fmt_string.format(value=abs(value),
                                 number_format=self._number_format,
                                 degree=self._degree_symbol,
                                 hemisphere=hemisphere)

    def _apply_transform(self, value, target_proj, source_crs):
        """
        Given a single value, a target projection and a source CRS,
        transforms the value from the source CRS to the target
        projection, returning a single value.

        """
        raise NotImplementedError("A subclass must implement this method.")

    def _hemisphere(self, value, value_source_crs):
        """
        Given both a tick value in the Plate Carree projection and the
        same value in the source CRS returns a string indicating the
        hemisphere that the value is in.

        Must be over-ridden by the derived class.

        """
        raise NotImplementedError("A subclass must implement this method.")


class LatitudeFormatter(_PlateCarreeFormatter):
    """Tick formatter for latitude axes."""
    def __init__(self, degree_symbol=u'\u00B0', number_format='g',
                 transform_precision=1e-8):
        """
        Tick formatter for a latitude axis.

        The axis must be part of an axes defined on a rectangular
        projection (e.g. Plate Carree, Mercator).

        .. note::

           A formatter can only be used for one axis. A new formatter
           must be created for every axis that needs formatted labels.

        Kwargs:

        * degree_symbol (string):
            The character(s) used to represent the degree symbol in the
            tick labels. Defaults to u'\u00B0' which is the unicode
            degree symbol. Can be an empty string if no degree symbol is
            desired.

        * number_format (string):
            Format string to represent the tick values. Defaults to 'g'.

        * transform_precision (float):
            Sets the precision (in degrees) to which transformed tick
            values are rounded. The default is 1e-7, and should be
            suitable for most use cases. To control the appearance of
            tick labels use the *number_format* keyword.

        Examples:

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

        """
        super(LatitudeFormatter, self).__init__(
            degree_symbol=degree_symbol,
            number_format=number_format,
            transform_precision=transform_precision)

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
                 number_format='g',
                 transform_precision=1e-8):
        """
        Create a formatter for longitude values.

        The axis must be part of an axes defined on a rectangular
        projection (e.g. Plate Carree, Mercator).

        .. note::

           A formatter can only be used for one axis. A new formatter
           must be created for every axis that needs formatted labels.

        Kwargs:

        * zero_direction_label (False | True):
            If *True* a direction label (E or W) will be drawn next to
            longitude labels with the value 0. If *False* then these
            labels will not be drawn. Defaults to *False* (no direction
            labels).

        * dateline_direction_label (False | True):
            If *True* a direction label (E or W) will be drawn next to
            longitude labels with the value 180. If *False* then these
            labels will not be drawn. Defaults to *False* (no direction
            labels).

        * degree_symbol (string):
            The symbol used to represent degrees. Defaults to u'\u00B0'
            which is the unicode degree symbol.

        * number_format (string):
            Format string to represent the longitude values. Defaults to
            'g'.

        * transform_precision (float):
            Sets the precision (in degrees) to which transformed tick
            values are rounded. The default is 1e-7, and should be
            suitable for most use cases. To control the appearance of
            tick labels use the *number_format* keyword.

        Examples:

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
            ont_formatter = LongitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)

        """
        super(LongitudeFormatter, self).__init__(
            degree_symbol=degree_symbol,
            number_format=number_format,
            transform_precision=transform_precision)
        self._zero_direction_labels = zero_direction_label
        self._dateline_direction_labels = dateline_direction_label

    def _apply_transform(self, value, target_proj, source_crs):
        return target_proj.transform_point(value, 0, source_crs)[0]

    def _hemisphere(self, value, value_source_crs):
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
