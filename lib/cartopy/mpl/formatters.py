# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.


import matplotlib.ticker as mticker
import numpy as np

_DEGREE_SYMBOL = '\u00B0'


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


class BASE_CLASS_FORMATER():
    def __init__(self,
                 north_hemisphere_str='N',
                 south_hemisphere_str='S',
                 west_hemisphere_str='W',
                 east_hemisphere_str='E',
                 degree_symbol=None,
                 gridline_num_format='g',
                 decimal_separator='.'):

        self._west_hemisphere_str = west_hemisphere_str
        self._east_hemisphere_str = east_hemisphere_str
        self._north_hemisphere_str = north_hemisphere_str
        self._south_hemisphere_str = south_hemisphere_str

        if degree_symbol is None:
            self._DEGREE_SYMBOL = _DEGREE_SYMBOL
        else:
            self._DEGREE_SYMBOL = degree_symbol
            
        self._gridline_num_format = gridline_num_format
        
        self._decimal_separator = decimal_separator
    
    def _reset(self):
        self.__init__(north_hemisphere_str=self.north_hemisphere_str,
                      south_hemisphere_str=self.south_hemisphere_str,
                      west_hemisphere_str=self.west_hemisphere_str,
                      east_hemisphere_str=self.east_hemisphere_str,
                      degree_symbol=self._DEGREE_SYMBOL,
                      gridline_num_format=self.gridline_num_format)
    
    @property
    def west_hemisphere_str(self):
        return self._west_hemisphere_str
        
    @west_hemisphere_str.setter
    def west_hemisphere_str(self, string):

        self._west_hemisphere_str = string
        #self._reset()
        
        
    @property
    def decimal_separator(self):
        return self._decimal_separator
    
    @decimal_separator.setter
    def decimal_separator(self, decimal_separator):
        self._decimal_separator = decimal_separator    
        
    @property
    def east_hemisphere_str(self):
        return self._east_hemisphere_str
        
    @east_hemisphere_str.setter
    def east_hemisphere_str(self, string):
        
        self._east_hemisphere_str = string
        #self._reset()

    @property
    def north_hemisphere_str(self):
        return self._north_hemisphere_str
        
    @north_hemisphere_str.setter
    def north_hemisphere_str(self, string):

        self._north_hemisphere_str = string
        #self._reset()
        
    @property
    def south_hemisphere_str(self):
        return self._south_hemisphere_str
        
    @south_hemisphere_str.setter
    def south_hemisphere_str(self, string):
    
        self._south_hemisphere_str = string
        #self._reset()

    @property
    def gridline_num_format(self):
        return self._gridline_num_format
    
    @gridline_num_format.setter
    def gridline_num_format(self, num_format):

        self._gridline_num_format = num_format
        #self._reset()

    def _lon_hemisphere(self, longitude):
        """Return the hemisphere (E, W or '' for 0) for the given longitude."""
        longitude = _fix_lons(longitude)
        if longitude > 0:
            hemisphere = self.east_hemisphere_str
        elif longitude < 0:
            hemisphere = self.west_hemisphere_str

        elif np.isclose(longitude, 0, rtol=1e-7):
            hemisphere = self.east_hemisphere_str
        else:
            hemisphere = ''
        return hemisphere

    def _lat_hemisphere(self, latitude):
        """Return the hemisphere (N, S or '' for 0) for the given latitude."""

        latitude = _fix_lats(latitude)
        if latitude > 0:
            hemisphere = self.north_hemisphere_str
        elif latitude < 0:
            hemisphere = self.south_hemisphere_str

        elif np.isclose(latitude, 0, rtol=1e-7):
            hemisphere = self.north_hemisphere_str

        else:
            hemisphere = ''
        return hemisphere

    def _north_south_formatted(self, latitude):

        if self.gridline_num_format == 'g':
            fmt_string = u'{latitude:{gridline_num_format}}{degree}{hemisphere}'
            return fmt_string.format(latitude=abs(latitude), gridline_num_format=self.gridline_num_format,
                                     hemisphere=self._lat_hemisphere(latitude),
                                     degree=self._DEGREE_SYMBOL)

        else:
            coords = '{}'.format(self.gridline_num_format).format(latitude)

            coord_info = u'{degree}{hemisphere}'.format(degree=self._DEGREE_SYMBOL,
                                                        hemisphere=self._lat_hemisphere(latitude),
                                                        )

            fmt_string = ''.join([coords , coord_info])

        return  fmt_string.replace('.', self.decimal_separator)

    def _east_west_formatted(self, longitude):
        
        if self.gridline_num_format == 'g':
            fmt_string = u'{longitude:{gridline_num_format}}{degree}{hemisphere}'
        
            fmt_string = fmt_string.format(longitude=abs(longitude), gridline_num_format=self.gridline_num_format,
                                     hemisphere=self._lon_hemisphere(longitude),
                                     degree=self._DEGREE_SYMBOL)
        else:

            coords = '{}'.format(self.gridline_num_format).format(longitude)

            coord_info = u'{degree}{hemisphere}'.format(degree=self._DEGREE_SYMBOL,
                                                        hemisphere=self._lon_hemisphere(longitude)
                                                        )

            fmt_string = ''.join([coords , coord_info])

        return fmt_string.replace('.', self.decimal_separator)


class LONGITUDE_FORMATTER_CLASS(BASE_CLASS_FORMATER):

    def __init__(self,
                 west_hemisphere_str='W',
                 east_hemisphere_str='E',
                 degree_symbol=None,
                 gridline_num_format='g',
                 decimal_separator='.'):

        BASE_CLASS_FORMATER.__init__(
            self,
            west_hemisphere_str=west_hemisphere_str,
            east_hemisphere_str=east_hemisphere_str,
            degree_symbol=degree_symbol,
            gridline_num_format=gridline_num_format,
            decimal_separator = decimal_separator)

        xformatter = mticker.FuncFormatter(lambda v, pos:
                                           self._east_west_formatted(v)
                                           )
        
        self.xformatter = xformatter


    def __call__(self, x, pos=None):

        return self.xformatter(x, pos)

    def set_locs(self, line_ticks):
        self.xformatter.set_locs(line_ticks)


class LATITUDE_FORMATTER_CLASS(BASE_CLASS_FORMATER):

    def __init__(self,
                 north_hemisphere_str='N',
                 south_hemisphere_str='S',
                 degree_symbol=None,
                 gridline_num_format='{0:3f}',
                 decimal_separator='.'):
        
        BASE_CLASS_FORMATER.__init__(
            self,
            north_hemisphere_str=north_hemisphere_str,
            south_hemisphere_str=south_hemisphere_str,
            degree_symbol=degree_symbol,
            gridline_num_format=gridline_num_format,
            decimal_separator = decimal_separator)


        yformatter = mticker.FuncFormatter(lambda v, pos:
                                                self._north_south_formatted(v)
                                           )

        self.yformatter = yformatter


    def __call__(self, x, pos=None):
        
        return self.yformatter(x, pos)

    def set_locs(self, line_ticks):
        
        self.yformatter.set_locs(line_ticks)

#########################################################


class Gridline_Base():
    
    def __init__(self,
                 west_hemisphere_str='W',
                 east_hemisphere_str='E',
                 north_hemisphere_str='N',
                 south_hemisphere_str='S',
                 degree_symbol=_DEGREE_SYMBOL,
                 gridline_num_format='{0:.2f}'):
        
        self.LONGITUDE_FORMATTER = LONGITUDE_FORMATTER_CLASS(
                            degree_symbol=degree_symbol,
                            west_hemisphere_str=west_hemisphere_str,
                            east_hemisphere_str=east_hemisphere_str,
                            gridline_num_format=gridline_num_format)
        
        
        self.LATITUDE_FORMATTER = LATITUDE_FORMATTER_CLASS(
                            degree_symbol=degree_symbol,
                            north_hemisphere_str=north_hemisphere_str,
                            south_hemisphere_str=south_hemisphere_str,
                            gridline_num_format=gridline_num_format)
        
        self.xformatter = self.LONGITUDE_FORMATTER
        self.yformatter = self.LATITUDE_FORMATTER

    def set_number_of_ticks(self, nbins=4, locator='xlocator'):
        '''
            Description:
                This is a helper function for setting (inplace) \
                the maximum number of ticks in the gridliner

            Parameters:
                nbins (int): number of bins to use in the given locator (axis)

                locator(str): the gridliner locator to set the number of bins
                    Standard value: 'xlocator'

            Return
                self (gridliner instance)

        '''
        locator = getattr(self, locator)

        locator = mticker.MaxNLocator(nbins)

        return self

    def change_gridline_tick_decimal_separator(
            self,
            gridline_num_format='{0:.2f}',
            axis='both',
            decimal_separator=',',
            degree_symbol=_DEGREE_SYMBOL):
        
        self.LONGITUDE_FORMATTER.decimal_separator = decimal_separator
        self.LATITUDE_FORMATTER.decimal_separator = decimal_separator
        
        self.LONGITUDE_FORMATTER.degree_symbol = degree_symbol
        self.LATITUDE_FORMATTER.degree_symbol = degree_symbol

        if axis.lower() == 'both':
            self.LONGITUDE_FORMATTER.gridline_num_format = gridline_num_format
            self.LATITUDE_FORMATTER.gridline_num_format = gridline_num_format
            
           

        else:
            if axis.lower().startswith('x'):
                self.LONGITUDE_FORMATTER.gridline_num_format = gridline_num_format
            else:
                self.LATITUDE_FORMATTER.gridline_num_format = gridline_num_format
        
        self.xformatter = self.LONGITUDE_FORMATTER
        self.yformatter = self.LATITUDE_FORMATTER

        return self

    def set_longitude_hemisphere_str(self,
                                     west_hemisphere_str='W',
                                     east_hemisphere_str='E'
                                     ):
        
        self.xformatter.west_hemisphere_str = west_hemisphere_str
        self.xformatter.east_hemisphere_str = east_hemisphere_str
        
        
    def set_latitude_hemisphere_str(self,
                                    north_hemisphere_str='N',
                                    south_hemisphere_str='S'
                                    ):
        
        self.yformatter.north_hemisphere_str = north_hemisphere_str
        self.yformatter.south_hemisphere_str = south_hemisphere_str
