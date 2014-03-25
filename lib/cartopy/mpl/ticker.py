# (C) British Crown Copyright 2014, Met Office
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
"""This module contains tools for handling tick marks in cartopy."""
from matplotlib.ticker import Formatter

import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes


class GeoFormatter(Formatter):
    """Base class for formatting ticks on geographical axes."""

    def __init__(self, degree_symbol=u'\u00B0', number_format='g'):
        """
        Create a formatter for geographical axis values.

        Kwargs:

        * degree_symbol (string):
            The character(s) used to represent the degree symbol in the
            tick labels. Defaults to u'\u00B0' which is the unicode
            degree symbol. Can be an empty string if no degree symbol is
            desired.

        * number_format (string):
            Format string to represent the tick values. Defaults to 'g'.

        """
        self._degree_symbol = degree_symbol
        self._number_format = number_format

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
        target = ccrs.PlateCarree()
        projected_value = self.extract_transform_result(
            target.transform_point(*self.make_transform_args(value, source)))
        # Round the transformed values to the nearest 0.1 degree for display
        # purposes (transforms can introduce minor rounding errors that make
        # the tick values look bad).
        projected_value = round(10 * projected_value) / 10
        # Return the formatted values, the formatter has both the re-projected
        # tick value and the original axis value available to it.
        return self._format_value(projected_value, value)

    def _format_value(self, value, original_value):
        hemisphere = self.hemisphere(value, original_value)
        fmt_string = u'{value:{number_format}}{degree}{hemisphere}'
        return fmt_string.format(value=abs(value),
                                 number_format=self._number_format,
                                 degree=self._degree_symbol,
                                 hemisphere=hemisphere)

    def make_transform_args(self, value, source_crs):
        """
        Given a single coordinate value and a source `CRS` returns a
        3-tuple of arguments suitable for use by `CRS.transform_point`.

        Must be over-ridden by the derived class.

        """
        raise NotImplementedError("A subclass must implement this method.")

    def extract_transform_result(self, transform_result):
        """
        Given a 2-tuple returned from `CRS.transform_point` returns the
        required element.

        Must be over-ridden by the derived class.

        """
        raise NotImplementedError("A subclass must implement this method.")

    def hemisphere(self, value, value_source_crs):
        """
        Given both a tick value in the Plate Carree projection and the
        same value in the source CRS returns a string indicating the
        hemisphere that the value is in.

        Must be over-ridden by the derived class.

        """
        raise NotImplementedError("A subclass must implement this method.")


class LatitudeFormatter(GeoFormatter):
    """Tick formatter for latitude axes."""

    def make_transform_args(self, value, source_crs):
        return (0, value, source_crs)

    def extract_transform_result(self, transform_result):
        return transform_result[1]

    def hemisphere(self, value, value_source_crs):
        if value > 0:
            hemisphere = 'N'
        elif value < 0:
            hemisphere = 'S'
        else:
            hemisphere = ''
        return hemisphere


class LongitudeFormatter(GeoFormatter):
    """Tick formatter for longitude axes."""

    def __init__(self,
                 zero_direction_label=False,
                 dateline_direction_label=False,
                 degree_symbol=u'\u00B0',
                 number_format='g'):
        """
        Create a formatter for longitude values.

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

        """
        super(LongitudeFormatter, self).__init__(degree_symbol=degree_symbol,
                                                 number_format=number_format)
        self._zero_direction_labels = zero_direction_label
        self._dateline_direction_labels = dateline_direction_label

    def make_transform_args(self, value, source_crs):
        return (value, 0, source_crs)

    def extract_transform_result(self, transform_result):
        return transform_result[0]

    def hemisphere(self, value, value_source_crs):
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
