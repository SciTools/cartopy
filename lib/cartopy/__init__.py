# (C) British Crown Copyright 2011 - 2017, Met Office
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

from __future__ import (absolute_import, division, print_function)

__version__ = '0.15.0'
__document_these__ = ['config']

# Enable shapely performance enhancements
import shapely.speedups
if shapely.speedups.available:
    shapely.speedups.enable()


# Configuration
import os.path

# for the writable data directory (i.e. the one where new data goes), follow
# the XDG guidelines found at
# https://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html
_writable_dir = os.path.join(os.path.expanduser('~'), '.local', 'share')
_data_dir = os.path.join(os.environ.get("XDG_DATA_HOME", _writable_dir),
                         'cartopy')

config = {'pre_existing_data_dir': '',
          'data_dir': _data_dir,
          'repo_data_dir': os.path.join(os.path.dirname(__file__), 'data'),
          'downloaders': {},
          }
"""
The config dictionary stores global configuration values for cartopy.

In the first instance, the config is defined in ``cartopy/__init__.py``. It
is possible to provide site wide customisations by including a
``siteconfig.py`` file along with the cartopy source code. ``siteconfig.py``
should contain a function called ``update_config`` which takes the config
dictionary instance as its first and only argument (from where it is
possible to update the dictionary howsoever desired).

For users without write permission to the cartopy source directory, a package
called ``cartopy_userconfig`` should be made importable (consider putting it
in ``site.getusersitepackages()``) and should expose a
function called ``update_config`` which takes the config dictionary as its
first and only argument.


Keys in the config dictionary:

 * ``pre_existing_data_dir`` - the absolute path to a directory where standard
                               data (such as that from NaturalEarth) can be
                               found. If it is not found in this location
                               the ``data_dir`` config item will be used.

 * ``data_dir`` - the absolute path to a directory where standard data (such
                  as that from NaturalEarth) can be found. If it is not found
                  and the item is downloadable cartopy will download the
                  appropriate file(s) to a subdirectory of this directory,
                  therefore ``data_dir`` should be writable by the user.

 * ``repo_data_dir`` - the absolute path to the directory where the data
                       delivered with the cartopy repository is stored.
                       Typically this will only be set by OS packagers and
                       system administrators for site wide deployments.

 * ``downloaders`` - a dictionary mapping standard "specifications" to the
                     appropriate :class:`~cartopy.io.Downloader`. For further
                     documentation and an example see
                     :func:`cartopy.io.Downloader.from_config`.

"""  # n.b. docstring changes should be propagated to docs/source/cartopy.rst

del _data_dir
del _writable_dir


# Try importing a siteconfig file which exposes an update_config function,
# otherwise, fail gracefully.
try:
    from cartopy.siteconfig import update_config as _update_config
    _update_config(config)
except ImportError:
    pass


# Try importing a cartopy_userconfig file which exposes an update_config
# function, otherwise, fail gracefully.
try:
    from cartopy_userconfig import update_config as _update_config
    _update_config(config)
except ImportError:
    pass


# Commonly used sub-modules. Imported here to provide end-user
# convenience.
import cartopy.crs
import cartopy.feature
