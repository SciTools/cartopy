# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD license.
# See LICENSE in the root of the repository for full licensing details.

import os
from pathlib import Path

import numpy as np
from setuptools import Extension, setup


# The existence of a PKG-INFO directory is enough to tell us whether this is a
# source installation or not (sdist).
HERE = Path(__file__).parent
IS_SDIST = (HERE / 'PKG-INFO').exists()
FORCE_CYTHON = os.environ.get('FORCE_CYTHON', False)

USE_CYTHON = not IS_SDIST or FORCE_CYTHON
if USE_CYTHON:
    import Cython
    if Cython.__version__ < '0.29':
        raise ImportError(
            "Cython 0.29+ is required to install cartopy from source.")
    ext = '.pyx'
else:
    ext = '.cpp'

# Macros to enable Cython coverage
define_macros = []
if os.environ.get('CYTHON_COVERAGE'):
    define_macros.append(('CYTHON_TRACE_NOGIL', '1'))

extensions = [
    Extension(
        'cartopy.trace',
        [f'lib/cartopy/trace{ext}'],
        include_dirs=[np.get_include()],
        language='c++',
        define_macros=define_macros),
]

if USE_CYTHON:
    # We need to explicitly cythonize the extension in order
    # to control the Cython compiler_directives.
    from Cython.Build import cythonize
    compiler_directives = {"profile": True, "linetrace": True}
    extensions = cythonize(extensions, compiler_directives=compiler_directives)


# Main setup
# ==========
setup(
    ext_modules=extensions,
)
