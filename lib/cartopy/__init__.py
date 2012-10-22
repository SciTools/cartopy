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

# convenient name for creating projections
import cartopy.crs as prj

# Enable shapely performance enhancements
import shapely.speedups
if shapely.speedups.available:
    shapely.speedups.enable()

# Commonly used sub-modules. Imported here to provide end-user
# convenience.
import cartopy.crs
import cartopy.feature


__version__ = '0.5.x'
