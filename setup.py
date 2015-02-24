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

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension
from distutils.core import Command
from distutils.sysconfig import get_config_var
from distutils.util import convert_path
import fnmatch
import os
import subprocess
import sys
import warnings

from Cython.Distutils import build_ext
import numpy as np


PY3 = (sys.version_info[0] == 3)

# Please keep in sync with INSTALL file.
GEOS_MIN_VERSION = (3, 3, 3)

HERE = os.path.dirname(__file__)


def file_walk_relative(top, remove=''):
    """
    Returns a generator of files from the top of the tree, removing
    the given prefix from the root/file result.

    """
    top = top.replace('/', os.path.sep)
    remove = remove.replace('/', os.path.sep)
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
    # Accept a root_path with Linux path separators.
    root_path = root_path.replace('/', os.path.sep)
    root_count = len(root_path.split(os.path.sep))
    for (dir_path, dir_names, _) in os.walk(convert_path(root_path)):
        # Prune dir_names *in-place* to prevent unwanted directory recursion
        for dir_name in list(dir_names):
            if not os.path.isfile(os.path.join(dir_path, dir_name, '__init__.py')):
                dir_names.remove(dir_name)
        if dir_names:
            prefix = dir_path.split(os.path.sep)[root_count:]
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


# Dependency checks
# =================

# GEOS
try:
    geos_version = subprocess.check_output(['geos-config', '--version'])
    geos_version = tuple(int(v) for v in geos_version.split(b'.'))
    geos_includes = subprocess.check_output(['geos-config', '--includes'])
    geos_clibs = subprocess.check_output(['geos-config', '--clibs'])
except (OSError, ValueError, subprocess.CalledProcessError):
    warnings.warn(
        'Unable to determine GEOS version. Ensure you have %s or later '
        'installed, or installation may fail.' % (
            '.'.join(str(v) for v in GEOS_MIN_VERSION), ))

    geos_includes = []
    geos_library_dirs = []
    if sys.platform.startswith('win'):
        geos_libraries = ['geos']
    else:
        geos_libraries = ['geos_c']
else:
    if geos_version < GEOS_MIN_VERSION:
        print('GEOS version %s is installed, but cartopy requires at least '
              'version %s.' % ('.'.join(str(v) for v in geos_version),
                               '.'.join(str(v) for v in GEOS_MIN_VERSION)))
        exit(1)

    if PY3:
        geos_includes = geos_includes.decode()
        geos_clibs = geos_clibs.decode()

    geos_includes = geos_includes.split()
    geos_libraries = []
    geos_library_dirs = []
    for entry in geos_clibs.split():
        if entry.startswith('-L'):
            geos_library_dirs.append(entry[2:])
        elif entry.startswith('-l'):
            geos_libraries.append(entry[2:])

if sys.platform.startswith('win'):
    def get_config_var(name):
        return '.'
    extra_extension_args = {}
else:
    extra_extension_args = dict(
        runtime_library_dirs=[get_config_var('LIBDIR')])

# Description
# ===========

with open(os.path.join(HERE, 'README.rst'), 'r') as fh:
    description = ''.join(fh.readlines())

# Main setup
# ==========
setup(
    name='Cartopy',
    version='0.12.x',
    url='http://scitools.org.uk/cartopy/docs/latest/',
    download_url='https://github.com/SciTools/cartopy',
    author='UK Met Office',
    description='A cartographic python library with matplotlib support for visualisation',
    long_description=description,
    license = "LGPLv3",
    keywords = "cartography map transform projection proj.4 geos shapely shapefile",

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
                  include_dirs=[get_config_var('INCLUDEDIR'), './lib/cartopy'] + geos_includes,
                  libraries=['proj'] + geos_libraries,
                  library_dirs=[get_config_var('LIBDIR')] + geos_library_dirs,
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
    classifiers=[
            'Development Status :: 4 - Beta',
            'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Operating System :: POSIX :: AIX',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: C++',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: GIS',
            'Topic :: Scientific/Engineering :: Visualization',
          ],
)
