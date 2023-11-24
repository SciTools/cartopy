# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

import matplotlib
import packaging.version


_MPL_VERSION = packaging.version.parse(matplotlib.__version__)
_MPL_36 = _MPL_VERSION.release[:2] >= (3, 6)
_MPL_38 = _MPL_VERSION.release[:2] >= (3, 8)
