# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.


import matplotlib.ticker as mticker
import numpy as np
import re

_DEGREE_SYMBOL = '\u00B0'


class CardinalDirections():
    """
    This class holds all cardination symbols
    applied in the gridline.
    """

    def __init__(self,
                 top_label_symbol='N',
                 bottom_label_symbol='S',
                 left_label_symbol='W',
                 right_label_symbol='E'):

        self._west_hemisphere_symbol = left_label_symbol
        self._east_hemisphere_symbol = right_label_symbol
        self._north_hemisphere_symbol = top_label_symbol
        self._south_hemisphere_symbol = bottom_label_symbol

    @property
    def west_hemisphere_symbol(self):
        return self._west_hemisphere_symbol

    @west_hemisphere_symbol.setter
    def west_hemisphere_symbol(self, new_symbol):
        self._west_hemisphere_symbol = new_symbol

    @property
    def east_hemisphere_symbol(self):
        return self._east_hemisphere_symbol

    @east_hemisphere_symbol.setter
    def east_hemisphere_symbol(self, new_symbol):
        self._east_hemisphere_symbol = new_symbol

    @property
    def north_hemisphere_symbol(self):
        return self._north_hemisphere_symbol

    @north_hemisphere_symbol.setter
    def north_hemisphere_symbol(self, new_symbol):
        self._north_hemisphere_symbol = new_symbol

    @property
    def south_hemisphere_symbol(self):
        return self._south_hemisphere_symbol

    @south_hemisphere_symbol.setter
    def south_hemisphere_symbol(self, new_symbol):
        self._south_hemisphere_symbol = new_symbol


def _fix_lons(lons):
    """
    Fix the given longitudes into the range ``[-180, 180]``.

    """
    lons = np.array(lons, copy=False, ndmin=1)
    fixed_lons = ((lons + 180) % 360) - 180
    # Make the positive 180s positive again.
    fixed_lons[(fixed_lons == -180) & (lons > 0)] *= -1
    return fixed_lons


def _fix_lats(lons):
    """
    Fix the given longitudes into the range ``[-180, 180]``.

    """
    lons = np.array(lons, copy=False, ndmin=1)
    fixed_lons = ((lons + 90) % 180) - 90
    # Make the positive 180s positive again.
    fixed_lons[(fixed_lons == -90) & (lons > 0)] *= -1
    return fixed_lons


class BASE_CLASS_FORMATER(CardinalDirections):
    def __init__(self,
                 north_hemisphere_symbol='N',
                 south_hemisphere_symbol='S',
                 west_hemisphere_symbol='W',
                 east_hemisphere_symbol='E',
                 degree_symbol=None,
                 decimal_separator='.'):

        CardinalDirections.__init__(self,
                                    north_hemisphere_symbol,
                                    south_hemisphere_symbol,
                                    west_hemisphere_symbol,
                                    east_hemisphere_symbol)

        if degree_symbol is None:
            self._DEGREE_SYMBOL = _DEGREE_SYMBOL
        else:
            self._DEGREE_SYMBOL = degree_symbol

        self._decimal_separator = decimal_separator

    @property
    def decimal_separator(self):
        return self._decimal_separator

    @decimal_separator.setter
    def decimal_separator(self, newstring):
        self._decimal_separator = newstring

    def _reset(self):
        self.__init__(north_hemisphere_symbol=self.north_hemisphere_symbol,
                      south_hemisphere_symbol=self.south_hemisphere_symbol,
                      west_hemisphere_symbol=self.west_hemisphere_symbol,
                      east_hemisphere_symbol=self.east_hemisphere_symbol,
                      degree_symbol=self._DEGREE_SYMBOL)

    def _lat_hemisphere(self, latitude):
        """Return the hemisphere (N, S or '' for 0) for the given latitude."""

        latitude = _fix_lats(latitude)
        if latitude > 0:
            hemisphere = self.north_hemisphere_symbol
        elif latitude < 0:
            hemisphere = self.south_hemisphere_symbol

        elif np.isclose(latitude, 0, rtol=1e-7):
            hemisphere = self.north_hemisphere_symbol

        else:
            hemisphere = ''
        return hemisphere

    def _north_south_formatted(self, latitude, latitude_gl_num_format):

        lat_hemisphere = self._lat_hemisphere(latitude)
        if latitude_gl_num_format == 'g':
            fmt_string = u'{latitude:{num_format}}\
                           {degree}{hemisphere}'

            String = fmt_string.format(latitude=abs(latitude),
                                       num_format=latitude_gl_num_format,
                                       hemisphere=lat_hemisphere,
                                       degree=self._DEGREE_SYMBOL)

            return String

        else:
            coords = '{}'.format(latitude_gl_num_format).format(abs(latitude))

            coord_info = (u'{degree}{hemisphere}'
                          .format(degree=self._DEGREE_SYMBOL,
                                  hemisphere=lat_hemisphere)
                          )

            fmt_string = ''.join([coords, coord_info])

        return fmt_string.replace('.', self.decimal_separator)

    def _lon_hemisphere(self, longitude):
        """Return the hemisphere (E, W or '' for 0) for the given longitude."""
        longitude = _fix_lons(longitude)
        if longitude > 0:
            hemisphere = self.east_hemisphere_symbol
        elif longitude < 0:
            hemisphere = self.west_hemisphere_symbol

        elif np.isclose(longitude, 0, rtol=1e-7):
            hemisphere = self.east_hemisphere_symbol
        else:
            hemisphere = ''
        return hemisphere

    def _east_west_formatted(self, longitude, num_format):

        lon_hemisphere = self._lon_hemisphere(longitude)

        if num_format == 'g':
            fmt_string = u'{longitude:{num_format}}\
            {degree}{hemisphere}'

            fmt_string = fmt_string.format(longitude=abs(longitude),
                                           num_format=num_format,
                                           hemisphere=lon_hemisphere,
                                           degree=self._DEGREE_SYMBOL)
        else:

            coords = '{}'.format(num_format).format(abs(longitude))

            coord_info = (u'{degree}{hemisphere}'
                          .format(degree=self._DEGREE_SYMBOL,
                                  hemisphere=lon_hemisphere)
                          )

            fmt_string = ''.join([coords, coord_info])

        return fmt_string.replace('.', self.decimal_separator)


class LONGITUDE_FORMATTER_CLASS():

    def __init__(self,
                 degree_symbol=None,
                 gl_num_format='g'):

        self._longitude_gl_num_format = gl_num_format

        def _func(v, pos):
            return self._east_west_formatted(v,
                                             self._longitude_gl_num_format
                                             )

        base_xformatter = mticker.FuncFormatter(_func)

        self.base_xformatter = base_xformatter

    def __call__(self, x, pos=None):

        return self.base_xformatter(x, pos)

    def set_locs(self, line_ticks):
        self.base_xformatter.set_locs(line_ticks)

    @property
    def longitude_gl_num_format(self):
        return self._longitude_gl_num_format

    @longitude_gl_num_format.setter
    def longitude_gl_num_format(self, num_format):

        self._longitude_gl_num_format = num_format


class LATITUDE_FORMATTER_CLASS():

    def __init__(self,
                 gl_num_format='{0:3f}'):

        self._latitude_gl_num_format = gl_num_format

        def _func(v, pos):
            return self._north_south_formatted(v,
                                               self._latitude_gl_num_format
                                               )

        base_yformatter = mticker.FuncFormatter(_func)

        self.base_yformatter = base_yformatter

    @property
    def latitude_gl_num_format(self):
        return self._latitude_gl_num_format

    @latitude_gl_num_format.setter
    def latitude_gl_num_format(self, num_format):

        self._latitude_gl_num_format = num_format

    def __call__(self, x, pos=None):

        return self.base_yformatter(x, pos)

    def set_locs(self, line_ticks):

        self.base_yformatter.set_locs(line_ticks)


#########################################################


class Gridline_Base(BASE_CLASS_FORMATER, LONGITUDE_FORMATTER_CLASS,
                    LATITUDE_FORMATTER_CLASS):

    def __init__(self,
                 west_hemisphere_symbol='W',
                 east_hemisphere_symbol='E',
                 north_hemisphere_symbol='N',
                 south_hemisphere_symbol='S',
                 degree_symbol=_DEGREE_SYMBOL,
                 gl_num_format='{0:.2f}',
                 decimal_separator='.'):

        LONGITUDE_FORMATTER_CLASS.__init__(self,
                                           gl_num_format=gl_num_format)

        LATITUDE_FORMATTER_CLASS.__init__(self,
                                          gl_num_format=gl_num_format)

        BASE_CLASS_FORMATER.__init__(self,
                                     north_hemisphere_symbol,
                                     south_hemisphere_symbol,
                                     west_hemisphere_symbol,
                                     east_hemisphere_symbol,
                                     degree_symbol,
                                     decimal_separator)

    def change_gridline_tick_decimal_separator(
            self,
            gl_num_format='{0:.2f}',
            axis='both',
            decimal_separator=',',
            degree_symbol=_DEGREE_SYMBOL):
        """
        Description:
            This function changes how the gridline ticklabels are drawn.
            Through this function, one can change the n° of decimal cases,
            the decimal seperator symbol (i.e.: dots - international symbol,
            commas for latin countries, etc.), as well as the degree symbol
            (i.e.: the '°')

        Parameters:

            gl_num_format (string): it sets how the ticklabels are drawn,
                and how many decimal cases to use.
                Default: '{0:.2f}' -> float with 2 decimal cases

            axis (string): it sets in which gridline border
                           to apply the changes
                if   axis == 'both': all borders are changed
                elif axis == 'x': bottom and top borders are changed
                elif axis == 'y': right and left borders are changed

                Accepted inputs: 'both', 'x', 'y':

                Default: 'both' -> both y and x borders of the gridline are
                    altered.

            decimal_separator (string): the decimal symbol to use.
                Default: '.'

            degree_symbol (string): the degree symbol to use.
                Default: '°'
        Returns
            self
        """

        self.decimal_separator = decimal_separator
        self.decimal_separator = decimal_separator

        self.degree_symbol = degree_symbol
        self.degree_symbol = degree_symbol

        if axis.lower() == 'both':
            self.latitude_gl_num_format = gl_num_format
            self.longitude_gl_num_format = gl_num_format

        else:
            if axis.lower().startswith('x'):
                self.longitude_gl_num_format = gl_num_format
            else:
                self.latitude_gl_num_format = gl_num_format

        return self

    def set_cardinal_directions(self,
                                new_cardinal_directions
                                ):
        """
        Description:
            This function allows one to set the cardinal direction
            symbols (i.e.: 'E' for East, 'W' for West', 'N' for North,
            'S' for South')

        Parameters:
            new_cardinal_directions (string, list or dictionary):
                if string: the order of cardinal directions is assumed
                           as following: East-West-North-Ssouth

                if list: the same order as if string

                if dictionary:
                    a regex function is applied over the common patterns
                    which lookup for the cardinal directions:
                    Default example: {'left': 'W',
                                      'east':'E',
                                      'north':'N',
                                      'south':'S'}

        Return
            None

        """

        def is_in_searcher(searcher, string):

            return max(searcher.match(string).span()) > 0

        if isinstance(new_cardinal_directions, dict):

            west_searcher = re.compile('lef*| wes*|w*|o*|oes*')
            east_searcher = re.compile('righ*|eas*|e*|lest*')
            north_searcher = re.compile('north*|top*|up*|n*')
            south_searcher = re.compile('south*|bott*|low*|s*')

            for key, symbol in new_cardinal_directions.items():

                if is_in_searcher(east_searcher, key.lower()):

                    self.east_hemisphere_symbol = symbol

                elif is_in_searcher(west_searcher, key.lower()):

                    self.west_hemisphere_symbol = symbol

                elif is_in_searcher(north_searcher, key.lower()):

                    self.north_hemisphere_symbol = symbol

                elif is_in_searcher(south_searcher, key.lower()):

                    self.south_hemisphere_symbol = symbol
                else:
                    pass

        if isinstance(new_cardinal_directions, str):
            new_cardinal_directions = list(new_cardinal_directions)

        if isinstance(new_cardinal_directions, list):

            # we will assume that the cardinals are being passed
            # in the following order (east, west, north, south)

            self.east_hemisphere_symbol = new_cardinal_directions[0]
            self.west_hemisphere_symbol = new_cardinal_directions[1]

            self.north_hemisphere_symbol = new_cardinal_directions[2]
            self.south_hemisphere_symbol = new_cardinal_directions[3]
