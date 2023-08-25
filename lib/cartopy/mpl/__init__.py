# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib
import packaging.version


_MPL_VERSION = packaging.version.parse(matplotlib.__version__)
_MPL_34 = _MPL_VERSION.release[:2] >= (3, 4)
_MPL_35 = _MPL_VERSION.release[:2] >= (3, 5)
_MPL_38 = _MPL_VERSION.release[:2] >= (3, 8)

assert _MPL_34, 'Cartopy is only supported with Matplotlib 3.4 or greater.'
