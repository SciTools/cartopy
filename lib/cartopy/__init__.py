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

import os.path

# convenient name for creating projections
import cartopy.crs as prj

# Enable shapely performance enhancements
import shapely.speedups
if shapely.speedups.available:
    shapely.speedups.enable()

# Commonly used sub-modules. Imported here to provide end-user
# convenience.
import cartopy.crs
import cartopy.feature


__version__ = '0.5.x'


config = {'data_dir': os.path.join(os.path.dirname(__file__), 'data'),
          'downloaders': {},
          }
"""
The config dictionary stores global configuration values for cartopy.

In the first instance, the config is defined in ``cartopy/__init__``. From
there, it is possible to provide site wide customisations by including a
``siteconfig.py`` file, along with the cartopy source code, which contains
a function ``update_config`` which takes the config dictionary instance as its
first and only argument (from where it is possible to update the dictionary
howsoever desired).

For users without write permission to the cartopy source directory, a package
called ``cartopy_userconfig`` should be made importable (consider putting it
in ``site.getusersitepackages()``) and should expose a
function called ``update_config`` which takes the config dictionary as its
first and only argument.


Keys in the config dictionary:

 * ``data_dir`` - the absolute path to a directory where standard data (such
                  as that from NaturalEarth) can be found. If it is not found
                  and the item is downloadable cartopy will download the
                  appropriate file(s) to a subdirectory of this directory,
                  therefore ``data_dir`` should be writable by the user.

 * ``downloaders`` - a dictionary mapping standard "specifications" to the
                     appropriate :class:`~cartopy.io.Downloader`. For further
                     documentation and an example see
                     :func:`cartopy.io.Downloader.from_config`.

"""


# try importing a siteconfig file which exposes an update_config function,
# otherwise, fail gracefully.
try:
    from cartopy.siteconfig import update_config as _update_config
    _update_config(config)
except ImportError:
    pass


# try importing a cartopy_userconfig file which exposes an update_config
# function, otherwise, fail gracefully.
try:
    from cartopy_userconfig import update_config as _update_config
    _update_config(config)
except ImportError:
    pass
