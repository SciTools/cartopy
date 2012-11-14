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

import os
import string
import urllib2
import warnings

from cartopy import config


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


class _MissingKeyFormatter(string.Formatter):
    """
    Provides a formatter class which can handle missing keys in a 
    format string.
    
    """
    def get_value(self, key, args, kwargs):
        if isinstance(key, (int, long)):
            raise ValueError('Index based format not applicable.')
        else:
            return kwargs.get(key, '{' + key + '}')
    

class DownloadableItem(object):
    """
    Represents a resource, that can be configured easily, which knows
    how to acquire itself (perhaps via HTTP).
    
    The key interface methods is :meth:`path` - typically *all* external calls
    will be made to that method.
    
    Args:
    
        ``url_template`` - a string which can be formatted using the 
                           ``.format`` method with the ``format_dict``
                           dictionary passed to :meth:`path` to create
                           the full URL of representing this resource.
                           
        ``target_path_template`` - a string which can be formatted using the
                                   ``.format`` method with the ``format_dict``
                                   dictionary passed to :meth:`path` to create
                                   the path of the target file representing
                                   this resource. If this file doesn't exist
                                   it will be downloaded via HTTP (subclasses 
                                   may implement alternative protocols).
                                   
    Kwargs:
    
        ``pre_downloaded_path_template`` - a string which can be formatted 
                                           using the ``.format`` method with 
                                           the ``format_dict`` dictionary 
                                           passed to :meth:`path` to create
                                           the path of a file which represents
                                           this resource. If the file does not
                                           exist, it will not necessarily be
                                           downloaded - instead target_path
                                           will be used to determine whether
                                           the resource already exists.
     
    """
    FORMAT_KEYS = ('config', )
    """
    The minimum keys which should be provided in the ``format_dict``
    argument for the ``path``, ``url``, ``target_path`` and 
    ``pre_downloaded_path`` methods.
    
    """
    
    def __init__(self, url_template, target_path_template, 
                 pre_downloaded_path_template=''):
        self.url_template = url_template
        self.target_path_template = target_path_template
        self.pre_downloaded_path_template = pre_downloaded_path_template
        # provide a caching mechanism for this DownloadableItem to store the
        # path of the file that it references. Will be None until self.path()
        # is called.
        self._cached_path = None
        """
        False if path caching is allowed for this DownloadableItem, None if
        it is. This will be updated by :meth:`path` if not False.
        
        """
        
    def url(self, format_dict):
        """
        The full URL that this resource represents.
        
        Args:
        
            ``format_dict`` - The dictionary which is used to replace
                              certain template variables. Subclasses should
                              document which keys are expected as a minimum
                              in their ``FORMAT_KEYS`` class attribute.
        
        """
        return _MissingKeyFormatter().format(self.url_template, **format_dict) 
    
    def target_path(self, format_dict):
        """
        The path on disk of the file that this resource represents, must
        either exist, or be writable by the current user. This method
        does not check either of these conditions.        
        
        Args:
        
            ``format_dict`` - The dictionary which is used to replace
                              certain template variables. Subclasses should
                              document which keys are expected as a minimum
                              in their ``FORMAT_KEYS`` class attribute.
        
        """
        return _MissingKeyFormatter().format(self.target_path_template,
                                             **format_dict)
        
    def pre_downloaded_path(self, format_dict):
        """
        The path on disk of the file that this resource represents, if it does
        not exist, then no further action will be taken with this path, and all
        further processing will be done using :meth:`target_path` instead.
        
        Args:
        
            ``format_dict`` - The dictionary which is used to replace
                              certain template variables. Subclasses should
                              document which keys are expected as a minimum
                              in their ``FORMAT_KEYS`` class attribute.
                              
        """
        return _MissingKeyFormatter().format(self.pre_downloaded_path_template,
                                             **format_dict)
        
    def path(self, format_dict):
        """
        Returns the path to a file on disk that this resource represents.
        
        If the file doesn't exist in :meth:`pre_downloaded_path` then it
        will check whether it exists in :meth:`target_path`, otherwise
        the resource will be downloaded via :meth:`acquire_resouce` from
        :meth:`url` to :meth:`target_path`.
        
        Typically, this is the method that most applications will call, 
        allowing implementors of new DownloadableItems to specialise
        :meth:`acquire_resource`.
        
        Args:
        
            ``format_dict`` - The dictionary which is used to replace
                              certain template variables. Subclasses should
                              document which keys are expected as a minimum
                              in their ``FORMAT_KEYS`` class attribute.
                              
        .. note::
        
            After the first call, unless ``_cached_path`` is False, the result will
            be cached, meaning the same path is always returned.
        
        """
        result_path = self._cached_path
        if result_path is None or result_path is False:
            pre_downloaded_path = self.pre_downloaded_path(format_dict)
            target_path = self.target_path(format_dict)
            if pre_downloaded_path is not None and os.path.exists(pre_downloaded_path):
                result_path = pre_downloaded_path
            elif os.path.exists(target_path):
                result_path = target_path
            else:
                # we need to download the file
                result_path = self.acquire_resource(target_path, format_dict)
             
            if self._cached_path is not False:    
                self._cached_path = result_path
        
        return self._cached_path
    
    def acquire_resource(self, target_path, format_dict):
        """
        Downloads, via HTTP, the file that this resource represents.
        Subclasses will typically override this method.
        
        Args:
        
            ``format_dict`` - The dictionary which is used to replace
                              certain template variables. Subclasses should
                              document which keys are expected as a minimum
                              in their ``FORMAT_KEYS`` class attribute.
        
        """
        target_dir = os.path.dirname(target_path)
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
        
        url = self.url(format_dict)
        
        # try getting the resource (no exception handling, just let it raise)
        response = self._urlopen(url)
            
        with open(target_path, 'wb') as fh:
            fh.write(response.read())
            
        return target_path
    
    def _urlopen(self, url):
        """
        Return a file handle to the given HTTP resource URL.

        Caller should close the file handle when finished with it.
        
        """ 
        warnings.warn('Downloading: {}'.format(url))
        return urllib2.urlopen(url)
    
    @staticmethod
    def from_config(specification):
        """
        The ``from_config`` static method implements the logic for acquiring a
        DownloadableItem (sub)class instance from ``cartopy.config``.
        
        The given specification should be iterable, as it will be traversed
        in reverse order before it finds the appropriate   
        
        Returns the appropriate DownloadableItem for the given specification.
        Looks at all levels of the specification, if there isn't one at the top
        level to start with, one will be created (and returned). 
        
        """
        spec_depth = len(specification)
        downloaders = config['downloads']
        result_downloader = None        
        
        for i in range(spec_depth, 0, -1):
            lookup = specification[:i]
            downloadable_item = downloaders.get(lookup, None)
            if downloadable_item is not None:
                # put the DownloadableItem at the top level 
                # (so that it is quicker to find next time).
                if i < spec_depth:
                    import copy
                    result_downloader = copy.copy(downloadable_item)
                    downloaders[specification] = result_downloader
                else:
                    result_downloader = downloadable_item
                return result_downloader
        else:
            # should never really happen, but could if the user does
            # some strange things...
            raise ValueError('No generic downloadable item in the '
                             '{} for {}'.format("config['downloads']",
                                                specification))
