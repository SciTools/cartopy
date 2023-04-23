# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

# NOTE: This file must remain Python 2 compatible for the foreseeable future,
# to ensure that we error out properly for people with outdated setuptools
# and/or pip.
import sys


PYTHON_MIN_VERSION = (3, 8)

if sys.version_info < PYTHON_MIN_VERSION:
    error = """
Beginning with Cartopy 0.21, Python {} or above is required.
You are using Python {}.

This may be due to an out of date pip.

Make sure you have pip >= 9.0.1.
""".format('.'.join(str(n) for n in PYTHON_MIN_VERSION),
           '.'.join(str(n) for n in sys.version_info[:3]))
    sys.exit(error)


from collections import defaultdict
import os
from pathlib import Path
import subprocess
from sysconfig import get_config_var
import warnings

from setuptools import Extension, find_packages, setup


# The existence of a PKG-INFO directory is enough to tell us whether this is a
# source installation or not (sdist).
HERE = Path(__file__).parent
IS_SDIST = (HERE / 'PKG-INFO').exists()
FORCE_CYTHON = os.environ.get('FORCE_CYTHON', False)

if not IS_SDIST or FORCE_CYTHON:
    import Cython
    if Cython.__version__ < '0.29':
        raise ImportError(
            "Cython 0.29+ is required to install cartopy from source.")

    from Cython.Distutils import build_ext as cy_build_ext


try:
    import numpy as np
except ImportError:
    raise ImportError('NumPy 1.19+ is required to install cartopy.')


# Please keep in sync with INSTALL file.
GEOS_MIN_VERSION = (3, 7, 2)


def file_walk_relative(root):
    """
    Return a list of files from the top of the tree, removing
    the lib/cartopy prefix from the resulting paths.

    """
    return [str(p.relative_to(Path('lib') / 'cartopy'))
            for p in root.glob("**/*")]


# Dependency checks
# =================

# GEOS
try:
    geos_version = subprocess.check_output(['geos-config', '--version'])
    geos_version = tuple(int(v) for v in geos_version.split(b'.')
                         if 'dev' not in str(v))
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

    geos_includes = geos_includes.decode().split()
    geos_libraries = []
    geos_library_dirs = []
    for entry in geos_clibs.decode().split():
        if entry.startswith('-L'):
            geos_library_dirs.append(entry[2:])
        elif entry.startswith('-l'):
            geos_libraries.append(entry[2:])


# Python dependencies
extras_require = {}
for name in (HERE / 'requirements').iterdir():
    with open(name) as fh:
        section, ext = name.stem, name.suffix
        extras_require[section] = []
        for line in fh:
            if line.startswith('#'):
                pass
            elif line.startswith('-'):
                pass
            else:
                extras_require[section].append(line.strip())
install_requires = extras_require.pop('default')
tests_require = extras_require.get('tests', [])

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
with open(HERE / 'README.md') as fh:
    description = ''.join(fh.readlines())


cython_coverage_enabled = os.environ.get('CYTHON_COVERAGE', None)
if cython_coverage_enabled:
    extra_extension_args["define_macros"].append(
        ('CYTHON_TRACE_NOGIL', '1')
    )

extensions = [
    Extension(
        'cartopy.trace',
        ['lib/cartopy/trace.pyx'],
        include_dirs=([include_dir, './lib/cartopy', np.get_include()] +
                      geos_includes),
        libraries=geos_libraries,
        library_dirs=[library_dir] + geos_library_dirs,
        language='c++',
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
            spath = Path(sfile)
            ext = spath.suffix
            if ext in ('.pyx',):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = str(spath.with_suffix(ext))
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


if IS_SDIST and not FORCE_CYTHON:
    extensions = decythonize(extensions)
    cmdclass = {}
else:
    cmdclass = {'build_ext': cy_build_ext}

base_path = Path('lib') / 'cartopy'
package_data = (
    file_walk_relative(base_path / 'tests' / 'mpl' / 'baseline_images')
    + file_walk_relative(base_path / 'data')
    + file_walk_relative(base_path / 'tests' / 'lakes_shapefile')
    + ['io/srtm.npz'])

# Main setup
# ==========
setup(
    name='Cartopy',
    url='https://scitools.org.uk/cartopy/docs/latest/',
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

    use_scm_version={
        'write_to': 'lib/cartopy/_version.py',
    },

    packages=find_packages("lib"),
    package_dir={'': 'lib'},
    package_data={'cartopy': package_data},

    scripts=['tools/cartopy_feature_download.py'],
    ext_modules=extensions,
    cmdclass=cmdclass,
    python_requires='>=' + '.'.join(str(n) for n in PYTHON_MIN_VERSION),
    classifiers=[
            'Development Status :: 4 - Beta',
            'Framework :: Matplotlib',
            'License :: OSI Approved :: GNU Lesser General Public License v3 '
            'or later (LGPLv3+)',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Operating System :: POSIX :: AIX',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: C++',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3 :: Only',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: GIS',
            'Topic :: Scientific/Engineering :: Visualization',
          ],
)
