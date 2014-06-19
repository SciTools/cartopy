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
Distribution definition for Cartopy.

"""

from distutils.core import setup, Command, Extension
from distutils.sysconfig import get_config_var
from distutils.util import convert_path
import fnmatch
import os
import sys

from Cython.Distutils import build_ext
import numpy as np


if sys.platform.startswith('win'):
    def get_config_var(name):
        return '.'
    geos = 'geos'
    extra_extension_args = {}
else:
    geos = 'geos_c'
    extra_extension_args = dict(runtime_library_dirs=[get_config_var('LIBDIR')])


def file_walk_relative(top, remove=''):
    """
    Returns a generator of files from the top of the tree, removing
    the given prefix from the root/file result.

    """
    for root, dirs, files in os.walk(top):
        for file in files:
            yield os.path.join(root, file).replace(remove, '')


def find_package_tree(root_path, root_package):
    """
    Returns the package and all its sub-packages.

    Automated package discovery - extracted/modified from Distutils Cookbook:
    http://wiki.python.org/moin/Distutils/Cookbook/AutoPackageDiscovery

    """
    packages = [root_package]
    root_count = len(root_path.split('/'))
    for (dir_path, dir_names, _) in os.walk(convert_path(root_path)):
        # Prune dir_names *in-place* to prevent unwanted directory recursion
        for dir_name in list(dir_names):
            if not os.path.isfile(os.path.join(dir_path, dir_name, '__init__.py')):
                dir_names.remove(dir_name)
        if dir_names:
            prefix = dir_path.split('/')[root_count:]
            packages.extend(['.'.join([root_package] + prefix + [dir_name]) for dir_name in dir_names])
    return packages


class MissingHeaderError(Exception):
    """
    Raised when one or more files do not have the required copyright
    and licence header.

    """
    pass


class HeaderCheck(Command):
    """
    Checks that all the necessary files have the copyright and licence
    header.

    """

    description = "check for copyright/licence headers"
    user_options = []

    exclude_patterns = ('./setup.py',
                        './build/*',
                        './docs/build/*',
                        './dist/*',
                        './lib/cartopy/examples/*.py')

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        check_paths = []
        for root, dirs, files in os.walk('.'):
            for file in files:
                if file.endswith('.py') or file.endswith('.c'):
                    path = os.path.join(root, file)
                    check_paths.append(path)

        for pattern in self.exclude_patterns:
            exclude = lambda path: not fnmatch.fnmatch(path, pattern)
            check_paths = list(filter(exclude, check_paths))

        bad_paths = list(filter(self._header_bad, check_paths))
        if bad_paths:
            raise MissingHeaderError(bad_paths)

    def _header_bad(self, path):
        target = '(C) British Crown Copyright 2011 - 2012, Met Office'
        with open(path, 'rt') as text_file:
            # Check for the header on the first line.
            line = text_file.readline().rstrip()
            bad = target not in line

            # Check if it was an executable script, with the header
            # starting on the second line.
            if bad and line == '#!/usr/bin/env python':
                line = text_file.readline().rstrip()
                bad = target not in line
        return bad


setup(
    name='Cartopy',
    version='0.11.0',
    url='http://github.com/SciTools/cartopy',
    author='UK Met Office',

    packages=find_package_tree('lib/cartopy', 'cartopy'),
    package_dir={'': 'lib'},
    package_data={'cartopy': list(file_walk_relative('lib/cartopy/tests/mpl/baseline_images/',
                                                     remove='lib/cartopy/')) +\
                             list(file_walk_relative('lib/cartopy/data/raster',
                                                     remove='lib/cartopy/')) +\
                             list(file_walk_relative('lib/cartopy/data/netcdf',
                                                     remove='lib/cartopy/')) +\
                             list(file_walk_relative('lib/cartopy/data/wmts',
                                                     remove='lib/cartopy/')) +\
                             list(file_walk_relative('lib/cartopy/data/shapefiles/natural_earth',
                                                     remove='lib/cartopy/')) +\
                             list(file_walk_relative('lib/cartopy/data/shapefiles/gshhs',
                                                     remove='lib/cartopy/')) +\
                             ['io/srtm.json']
                 },


    # requires proj4 headers
    ext_modules=[
        Extension('cartopy.trace', ['lib/cartopy/trace.pyx', 'lib/cartopy/_trace.cpp'],
                  include_dirs=[get_config_var('INCLUDEDIR'), './lib/cartopy'],
                  libraries=[geos, 'proj'],
                  library_dirs=[get_config_var('LIBDIR')],
                  language='c++',
                  **extra_extension_args
                  ),
        Extension('cartopy._crs', ['lib/cartopy/_crs.pyx'],
                  include_dirs=[get_config_var('INCLUDEDIR'), np.get_include()],
                  libraries=['proj'],
                  library_dirs=[get_config_var('LIBDIR')],
                  **extra_extension_args
                  ),
    ],

    cmdclass={'build_ext': build_ext, 'header_check': HeaderCheck},
)
