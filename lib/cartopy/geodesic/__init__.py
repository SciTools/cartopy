# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

# NOTE: noqa for flake8 = unused import
from cartopy.geodesic._geodesic import Geodesic  # noqa: F401

# Need to be explicit because it's from a 'different' module.
__document_these__ = ['Geodesic']
