# (C) British Crown Copyright 2011 - 2015, Met Office
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


here = os.path.dirname(__file__)
with open(os.path.join(here, 'README.rst'), 'r') as fh:
    description = ''.join(fh.readlines())

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

    cmdclass={'build_ext': build_ext},
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
