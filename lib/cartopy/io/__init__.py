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
    Given a file handle, filename or (file handle, filename) tuple return
    (file handle, filename) in the given read/write mode.

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
