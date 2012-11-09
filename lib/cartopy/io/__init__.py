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
"""
Provides a collection of sub-packages for loading, saving and retrieving various data formats.

"""


def fh_getter(fh, mode='r', needs_filename=False):
    """
    Convenience function for opening files.
    
    Args:
    
        * fh - File handle, filename or (file handle, filename) tuple
    
    Kwargs:
    
        * mode - Open mode. Defaults to "r".
        
    Returns:

        * (file handle, filename), opened in the given mode.

    """
    if mode != 'r':
        raise ValueError('Only mode "r" currently supported.')

    if isinstance(fh, basestring):
        filename = fh
        fh = open(fh, mode)
    elif isinstance(fh, tuple):
        fh, filename = fh

    if filename is None:
        try:
            filename = fh.name
        except AttributeError: # does this occur?
            if needs_filename:
                raise ValueError('filename cannot be determined')
            else:
                filename = ''

    return fh, filename


import os
import urllib2

from cartopy import config


#config = {
#          'downloads': {}
#          }


import string
class MissingKeyFormatter(string.Formatter):
    def get_value(self, key, args, kwargs):
        if isinstance(key, (int, long)):
            raise ValueError('Index based format not applicable.')
        else:
            # XXX raise a warning?
            return kwargs.get(key, '{' + key + '}')
    

class DownloadableItem(object):
    """
    XXX 
    """
    def __init__(self, url_template, target_path_template, 
                 pre_downloaded_path_template=''):
        """
        url_template - the template of the downloadable item - 
                        formatted later on when given a full specification
                        
        target_path_template - the target location of the downloaded item - 
                           formatted later on when given a full specification
                           
        pre_downloaded_path_template - a location where a pre-downloaded version of
                                  this file might be found - formatted later on
                                   when given a full specification
        
        Ultimately you will want to use self.path to get the path of the resource.
        
        """
        self.url_template = url_template
        self.target_path_template = target_path_template
        self.pre_downloaded_path_template = pre_downloaded_path_template
        # provide a caching mechanism for this DownloadableItem to store the
        # path of the file that it references. Will be None until self.path()
        # is called.
        self._cached_path = None
        
    def url(self, format_dict):
        return MissingKeyFormatter().format(self.url_template, **format_dict) 
    
    def target_path(self, format_dict):
        """no guarantee it exists."""
        return MissingKeyFormatter().format(self.target_path_template,
                                            **format_dict)
        
    def pre_downloaded_path(self, format_dict):
        """no guarantee it exists."""
        return MissingKeyFormatter().format(self.pre_downloaded_path_template,
                                            **format_dict)
        
    def path(self, format_dict):
        """
        Return the path of this DownloadableItem, downloading the file if necessary.
        """ 
        if self._cached_path is None:
            pre_downloaded_path = self.pre_downloaded_path(format_dict)
            target_path = self.target_path(format_dict)
            if pre_downloaded_path is not None and os.path.exists(pre_downloaded_path):
                self._cached_path = pre_downloaded_path
            elif os.path.exists(target_path):
                self._cached_path = target_path
            else:
                # we need to download the file
                result = self._aquire_resource(target_path, format_dict)
                self._cached_path = result
        
        # XXX it may be that we want the result of _aquire_resource to be optional... (SRTM3)
        if isinstance(self._cached_path, Exception):
            raise self._cached_path 
        else:
            return self._cached_path
    
    def _aquire_resource(self, target_path, format_dict):
        """
        Downloads the resource that this DownloadableItem represents.
        """
        target_dir = os.path.dirname(target_path)
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        
        url = self.url(format_dict)
        
        # try getting the resource (no exception handling, just let it raise)
        response = urllib2.urlopen(url)
            
        with open(target_path, 'wb') as fh:
            fh.write(response.read())
            
        return target_path
    
    @staticmethod
    def from_config(specification):
        """
        Returns the appropriate DownloadableItem for the given specification.
        Looks at all levels of the specification, if there isn't one at the top
        level to start with, one will be created (and returned). 
        
        """
        spec_depth = len(specification)
        downloaders = config['downloads']
        
        for i in range(spec_depth, 0, -1):
            lookup = specification[:i]
#            print 'l:', lookup, lookup in downloaders, downloaders.keys()
            downloadable_item = downloaders.get(lookup, None)
            if downloadable_item is not None:
                # put the DownloadableItem at the top level 
                # (so that it is quicker to find next time).
                if i < spec_depth:
                    import copy
                    downloaders[specification] = copy.copy(downloadable_item)
                return downloadable_item
        else:
            # should never really happen, but could if the user does
            # some strange things...
            raise ValueError('No generic downloadable item in the '
                             '{} for {}'.format("config['downloads']",
                                                specification))
