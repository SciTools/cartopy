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
                 north_hemisphere_str=' - N',
                 south_hemisphere_str=' - S',
                 west_hemisphere_str=' - O', 
                 east_hemisphere_str=' - L',
                 degree_symbol=None):
        
        self.west_hemisphere_str = west_hemisphere_str
        self.east_hemisphere_str = east_hemisphere_str
        self.north_hemisphere_str = north_hemisphere_str
        self.south_hemisphere_str = south_hemisphere_str
        
        if degree_symbol is None:
            self._DEGREE_SYMBOL = _DEGREE_SYMBOL
        else:
            self._DEGREE_SYMBOL = degree_symbol
    
        
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
    
    
    def _east_west_formatted(self, longitude, num_format='g'):
        fmt_string = u'{longitude:{num_format}}{degree}{hemisphere}'.format(longitude=abs(longitude), 
        num_format=num_format,hemisphere=self._lon_hemisphere(longitude),
        degree=self._DEGREE_SYMBOL)
        
        
        return fmt_string
    
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
    
    def _north_south_formatted(self, latitude, num_format='g'):
        fmt_string = u'{latitude:{num_format}}{degree}{hemisphere}'
        return fmt_string.format(latitude=abs(latitude), num_format=num_format,
                                 hemisphere=self._lat_hemisphere(latitude),
                                 degree=self._DEGREE_SYMBOL)
    

class LONGITUDE_FORMATTER_CLASS(BASE_CLASS_FORMATER):
    
    def __init__(self, 
                 west_hemisphere_str=' - O', 
                 east_hemisphere_str=' - L',
                 degree_symbol=None):
    
        
        BASE_CLASS_FORMATER.__init__(self,west_hemisphere_str=west_hemisphere_str,
                                         east_hemisphere_str=east_hemisphere_str,
                                         degree_symbol=degree_symbol)

     
        
        self.XFormatter = mticker.FuncFormatter(lambda v, pos:
                                           self._east_west_formatted(v)
                                           )
        
    def __call__(self, lon, pos=None):
        
        return self.XFormatter(lon, pos)
    
    
    def set_locs(self, line_ticks):
        self.XFormatter.set_locs(line_ticks)

    
    
class LATITUDE_FORMATTER_CLASS(BASE_CLASS_FORMATER):
    
    def __init__(self, 
                 north_hemisphere_str=' - Norte',
                 south_hemisphere_str=' - SUl',
                 degree_symbol=None):
    
        
        BASE_CLASS_FORMATER.__init__(self, north_hemisphere_str=north_hemisphere_str,
                            south_hemisphere_str=south_hemisphere_str,
                            degree_symbol=degree_symbol)
        
        
        
        self.YFormatter = mticker.FuncFormatter(lambda v, pos:
                                           self._north_south_formatted(v)
                                           )
        
    def __call__(self, lat, lpos=None):
        
        return self.YFormatter(lat, lpos)
    
    
    def set_locs(self, line_ticks):
        self.YFormatter.set_locs(line_ticks)

#########################################################





class Gridline_Base():

    def __init__(self):
        pass
    
        
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
       

    def change_gridline_tick_decimal_separator(self, 
                                               gridline_tick_formating='{0:.2f}', 
                                               axis='yaxis', 
                                               decimal_separator=',', 
                                               geographical_symbol=_DEGREE_SYMBOL):  
         
        def custom_tick_func_formatter(x, y):
            ticklabel = '{0}'.format(gridline_tick_formating).format(x) 
            
            return ticklabel.replace('.', decimal_separator) + geographical_symbol
         
        if axis.lower() == 'both':
     
            self.xformatter = mticker.FuncFormatter(lambda x, y : custom_tick_func_formatter(x, y))
            self.yformatter = mticker.FuncFormatter(lambda x, y : custom_tick_func_formatter(x, y))
                 
        else:
            if axis.lower().startswith('x'):
                self.xformatter = mticker.FuncFormatter(lambda x, y : custom_tick_func_formatter(x, y))
            else:
                self.yformatter = mticker.FuncFormatter(lambda x, y : custom_tick_func_formatter(x, y))
        
        
        return self



    def set_longitude_hemisphere_str(self,
                                     west_hemisphere_str='W',
                                     east_hemisphere_str='E',
                                     degree_symbol=_DEGREE_SYMBOL
                                     ):
        
            
        lon_formatter = LONGITUDE_FORMATTER_CLASS(degree_symbol=degree_symbol,
                                           west_hemisphere_str=west_hemisphere_str,
                                           east_hemisphere_str=east_hemisphere_str)
                                           
        self.xformatter = lon_formatter



    
    def set_latitude_hemisphere_str(self,
                                    north_hemisphere_str='N',
                                    south_hemisphere_str='S',
                                    degree_symbol=_DEGREE_SYMBOL
                                    ):
        
        lat_formatter = LATITUDE_FORMATTER_CLASS(degree_symbol=degree_symbol,
                                          north_hemisphere_str=north_hemisphere_str,
                                          south_hemisphere_str=south_hemisphere_str)
            
        self.yformatter = lat_formatter