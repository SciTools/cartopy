# (C) British Crown Copyright 2011 - 2019, Met Office
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
from __future__ import print_function

import fnmatch
import os
import subprocess
import sys
import warnings
from collections import defaultdict
from distutils.spawn import find_executable
from distutils.sysconfig import get_config_var

from setuptools import Command, Extension, convert_path, setup

import versioneer

"""
Distribution definition for Cartopy.

"""




# The existence of a PKG-INFO directory is enough to tell us whether this is a
# source installation or not (sdist).
HERE = os.path.dirname(__file__)
IS_SDIST = os.path.exists(os.path.join(HERE, 'PKG-INFO'))

if not IS_SDIST:
    import Cython
    if Cython.__version__ < '0.28':
        raise ImportError(
            "Cython 0.28+ is required to install cartopy from source.")

    from Cython.Distutils import build_ext as cy_build_ext


try:
    import numpy as np
except ImportError:
    raise ImportError('NumPy 1.10+ is required to install cartopy.')


PY3 = (sys.version_info[0] == 3)

# Please keep in sync with INSTALL file.
GEOS_MIN_VERSION = (3, 3, 3)
PROJ_MIN_VERSION = (4, 9, 0)


def file_walk_relative(top, remove=''):
    """
    Return a generator of files from the top of the tree, removing
    the given prefix from the root/file result.

    """
    top = top.replace('/', os.path.sep)
    remove = remove.replace('/', os.path.sep)
    for root, dirs, files in os.walk(top):
        for file in files:
            yield os.path.join(root, file).replace(remove, '')


def find_package_tree(root_path, root_package):
    """
    Return the package and all its sub-packages.

    Automated package discovery - extracted/modified from Distutils Cookbook:
    https://wiki.python.org/moin/Distutils/Cookbook/AutoPackageDiscovery

    """
    packages = [root_package]
    # Accept a root_path with Linux path separators.
    root_path = root_path.replace('/', os.path.sep)
    root_count = len(root_path.split(os.path.sep))
    for (dir_path, dir_names, _) in os.walk(convert_path(root_path)):
        # Prune dir_names *in-place* to prevent unwanted directory recursion
        for dir_name in list(dir_names):
            if not os.path.isfile(os.path.join(dir_path, dir_name,
                                               '__init__.py')):
                dir_names.remove(dir_name)
        if dir_names:
            prefix = dir_path.split(os.path.sep)[root_count:]
            packages.extend(['.'.join([root_package] + prefix + [dir_name])
                             for dir_name in dir_names])
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
    geos_libraries = ['geos_c']
else:
    if geos_version < GEOS_MIN_VERSION:
        print('GEOS version %s is installed, but cartopy requires at least '
              'version %s.' % ('.'.join(str(v) for v in geos_version),
                               '.'.join(str(v) for v in GEOS_MIN_VERSION)),
              file=sys.stderr)
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


# Proj
def find_proj_version_by_program(conda=None):
    proj = find_executable('proj')
    if proj is None:
        print(
            'Proj %s must be installed.' % (
                '.'.join(str(v) for v in PROJ_MIN_VERSION), ),
            file=sys.stderr)
        exit(1)

    if conda is not None and conda not in proj:
        print(
            'Proj %s must be installed in Conda environment "%s".' % (
                '.'.join(str(v) for v in PROJ_MIN_VERSION), conda),
            file=sys.stderr)
        exit(1)

    try:
        proj_version = subprocess.check_output([proj],
                                               stderr=subprocess.STDOUT)
        proj_version = proj_version.split()[1].split(b'.')
        proj_version = tuple(int(v.strip(b',')) for v in proj_version)
    except (OSError, IndexError, ValueError, subprocess.CalledProcessError):
        warnings.warn(
            'Unable to determine Proj version. Ensure you have %s or later '
            'installed, or installation may fail.' % (
                '.'.join(str(v) for v in PROJ_MIN_VERSION), ))
        proj_version = (0, 0, 0)

    return proj_version


def get_proj_libraries():
    """
    This function gets the PROJ libraries to cythonize with
    """
    proj_libraries = ["proj"]
    if os.name == "nt" and proj_version >= (6, 0, 0):
        proj_libraries = [
            "proj_{}_{}".format(proj_version[0], proj_version[1])
        ]
    return proj_libraries


conda = os.getenv('CONDA_DEFAULT_ENV')
if conda is not None and conda in sys.prefix:
    # Conda does not provide pkg-config compatibility, but the search paths
    # should be set up so that nothing extra is required. We'll still check
    # the version, though.
    proj_version = find_proj_version_by_program(conda)
    if proj_version < PROJ_MIN_VERSION:
        print(
            'Proj version %s is installed, but cartopy requires at least '
            'version %s.' % ('.'.join(str(v) for v in proj_version),
                             '.'.join(str(v) for v in PROJ_MIN_VERSION)),
            file=sys.stderr)
        exit(1)

    proj_includes = []
    proj_libraries = get_proj_libraries()
    proj_library_dirs = []

else:
    try:
        proj_version = subprocess.check_output(['pkg-config', '--modversion',
                                                'proj'],
                                               stderr=subprocess.STDOUT)
        proj_version = tuple(int(v) for v in proj_version.split(b'.'))
        proj_includes = subprocess.check_output(['pkg-config', '--cflags',
                                                 'proj'])
        proj_clibs = subprocess.check_output(['pkg-config', '--libs', 'proj'])
    except (OSError, ValueError, subprocess.CalledProcessError):
        proj_version = find_proj_version_by_program()
        if proj_version < PROJ_MIN_VERSION:
            print(
                'Proj version %s is installed, but cartopy requires at least '
                'version %s.' % ('.'.join(str(v) for v in proj_version),
                                 '.'.join(str(v) for v in PROJ_MIN_VERSION)),
                file=sys.stderr)
            exit(1)

        proj_includes = []
        proj_libraries = get_proj_libraries()
        proj_library_dirs = []
    else:
        if proj_version < PROJ_MIN_VERSION:
            print(
                'Proj version %s is installed, but cartopy requires at least '
                'version %s.' % ('.'.join(str(v) for v in proj_version),
                                 '.'.join(str(v) for v in PROJ_MIN_VERSION)),
                file=sys.stderr)
            exit(1)

        if PY3:
            proj_includes = proj_includes.decode()
            proj_clibs = proj_clibs.decode()

        proj_includes = [
            proj_include[2:] if proj_include.startswith('-I') else
            proj_include for proj_include in proj_includes.split()]

        proj_libraries = []
        proj_library_dirs = []
        for entry in proj_clibs.split():
            if entry.startswith('-L'):
                proj_library_dirs.append(entry[2:])
            elif entry.startswith('-l'):
                proj_libraries.append(entry[2:])

# Python dependencies
extras_require = {}
for name in os.listdir(os.path.join(HERE, 'requirements')):
    with open(os.path.join(HERE, 'requirements', name), 'r') as fh:
        section, ext = os.path.splitext(name)
        extras_require[section] = []
        for line in fh:
            if line.startswith('#'):
                pass
            elif line.startswith('-'):
                pass
            else:
                extras_require[section].append(line.strip())
install_requires = extras_require.pop('default')
tests_require = extras_require.pop('tests', [])

# General extension paths
if sys.platform.startswith('win'):
    def get_config_var(name):
        return '.'
include_dir = get_config_var('INCLUDEDIR')
library_dir = get_config_var('LIBDIR')
extra_extension_args = defaultdict(list)
if not sys.platform.startswith('win'):
    extra_extension_args["runtime_library_dirs"].append(
        get_config_var('LIBDIR')
    )

# Description
# ===========
with open(os.path.join(HERE, 'README.md'), 'r') as fh:
    description = ''.join(fh.readlines())


cython_coverage_enabled = os.environ.get('CYTHON_COVERAGE', None)
if proj_version >= (6, 0, 0):
    extra_extension_args["define_macros"].append(
        ('ACCEPT_USE_OF_DEPRECATED_PROJ_API_H', '1')
    )
if cython_coverage_enabled:
    extra_extension_args["define_macros"].append(
        ('CYTHON_TRACE_NOGIL', '1')
    )

extensions = [
    Extension(
        'cartopy.trace',
        ['lib/cartopy/trace.pyx'],
        include_dirs=([include_dir, './lib/cartopy', np.get_include()] +
                      proj_includes + geos_includes),
        libraries=proj_libraries + geos_libraries,
        library_dirs=[library_dir] + proj_library_dirs + geos_library_dirs,
        language='c++',
        **extra_extension_args),
    Extension(
        'cartopy._crs',
        ['lib/cartopy/_crs.pyx'],
        include_dirs=[include_dir, np.get_include()] + proj_includes,
        libraries=proj_libraries,
        library_dirs=[library_dir] + proj_library_dirs,
        **extra_extension_args),
    # Requires proj v4.9
    Extension(
        'cartopy.geodesic._geodesic',
        ['lib/cartopy/geodesic/_geodesic.pyx'],
        include_dirs=[include_dir, np.get_include()] + proj_includes,
        libraries=proj_libraries,
        library_dirs=[library_dir] + proj_library_dirs,
        **extra_extension_args),
]


if cython_coverage_enabled:
    # We need to explicitly cythonize the extension in order
    # to control the Cython compiler_directives.
    from Cython.Build import cythonize

    directives = {'linetrace': True,
                  'binding': True}
    extensions = cythonize(extensions, compiler_directives=directives)


def decythonize(extensions, **_ignore):
    # Remove pyx sources from extensions.
    # Note: even if there are changes to the pyx files, they will be ignored.
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in ('.pyx',):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


cmdclass = versioneer.get_cmdclass()

if IS_SDIST:
    extensions = decythonize(extensions)
else:
    cmdclass.update({'build_ext': cy_build_ext})


# Main setup
# ==========
setup(
    name='Cartopy',
    version=versioneer.get_version(),
    url='http://scitools.org.uk/cartopy/docs/latest/',
    download_url='https://github.com/SciTools/cartopy',
    author='UK Met Office',
    description='A cartographic python library with Matplotlib support for '
                'visualisation',
    long_description=description,
    long_description_content_type='text/markdown',
    license="LGPLv3",
    keywords="cartography map transform projection proj proj.4 geos shapely "
             "shapefile",

    install_requires=install_requires,
    extras_require=extras_require,
    tests_require=tests_require,

    packages=find_package_tree('lib/cartopy', 'cartopy'),
    package_dir={'': 'lib'},
    package_data={'cartopy': list(file_walk_relative('lib/cartopy/tests/'
                                                     'mpl/baseline_images/',
                                                     remove='lib/cartopy/')) +
                  list(file_walk_relative('lib/cartopy/data/raster',
                                          remove='lib/cartopy/')) +
                  list(file_walk_relative('lib/cartopy/data/netcdf',
                                          remove='lib/cartopy/')) +
                  list(file_walk_relative('lib/cartopy/data/'
                                          'shapefiles/gshhs',
                                          remove='lib/cartopy/')) +
                  list(file_walk_relative('lib/cartopy/tests/lakes_shapefile',
                                          remove='lib/cartopy/')) +
                  ['io/srtm.npz']},


    # requires proj headers
    ext_modules=extensions,
    cmdclass=cmdclass,
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*',
    classifiers=[
            'Development Status :: 4 - Beta',
            'License :: OSI Approved :: GNU Lesser General Public License v3 '
            'or later (LGPLv3+)',
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
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: GIS',
            'Topic :: Scientific/Engineering :: Visualization',
          ],
)
