# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

import os
from pathlib import Path
import sys


# asv runs the benchmarks without activating the conda environment, so the PROJ
# data path that conda's activation script sets is missing, and GDAL (used by
# fiona) fails to open the PROJ database. Set it the same way.
_proj_data = Path(sys.prefix) / 'share' / 'proj'
if _proj_data.is_dir():
    os.environ.setdefault('PROJ_DATA', str(_proj_data))
