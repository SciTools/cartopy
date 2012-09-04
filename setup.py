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

from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var
from distutils.util import convert_path
import os

from Cython.Distutils import build_ext
import numpy


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


setup(
    name='Cartopy',
    version='0.1',
    url='http://github.com/SciTools/cartopy',
    author='Philip Elson',
    author_email='pelson.pub@gmail.com',
      
    packages=find_package_tree('lib/cartopy', 'cartopy'),
    package_dir={'': 'lib'},
    #package_data={'cartopy': ['data/*']},
    # requires proj4 headers
    ext_modules=[
        Extension('cartopy.trace', ['lib/cartopy/trace.pyx', 'lib/cartopy/_trace.cpp'],
                  include_dirs=[get_config_var('INCLUDEDIR'), './lib/cartopy'],
                  libraries=['geos_c', 'proj'],
                  library_dirs=[get_config_var('LIBDIR')],
                  runtime_library_dirs=[get_config_var('LIBDIR')],
                  language='c++',
                  ),
        Extension('cartopy._crs', ['lib/cartopy/_crs.pyx'],
                  include_dirs=[get_config_var('INCLUDEDIR'), './lib/cartopy', numpy.get_include()],
                  libraries=['proj'],
                  library_dirs=[get_config_var('LIBDIR')],
                  runtime_library_dirs=[get_config_var('LIBDIR')],
                  #language='c++',
                  ),
    ],

    cmdclass={'build_ext': build_ext},
)
