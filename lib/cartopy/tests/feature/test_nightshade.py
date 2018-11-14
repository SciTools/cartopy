# (C) British Crown Copyright 2018, Met Office
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

from __future__ import (absolute_import, division, print_function)

from datetime import datetime

import pytest

from cartopy.feature.nightshade import _julian_day, _solar_position


def test_julian_day():
    # Using Vallado 2007 "Fundamentals of Astrodynamics and Applications"
    # Example 3.4
    dt = datetime(1996, 10, 26, 14, 20)
    jd = _julian_day(dt)
    assert pytest.approx(jd) == 2450383.09722222


# Using timeanddate website. Not very many digits of accuracy,
# but it should do OK for simple testing
# 1) Early Sept 2018
#    https://www.timeanddate.com/worldclock/sunearth.html
#    ?month=9&day=29&year=2018&hour=00&min=00&sec=0&n=&ntxt=&earth=0
# 2) Later Sept 2018
#    https://www.timeanddate.com/worldclock/sunearth.html
#    ?month=9&day=29&year=2018&hour=14&min=00&sec=0&n=&ntxt=&earth=0
# 3) Different year/month (Feb 1992)
#    https://www.timeanddate.com/worldclock/sunearth.html
#    ?month=2&day=14&year=1992&hour=00&min=0&sec=0&n=&ntxt=&earth=0
# 4) Test in the future near summer solstice (June 2030)
#    https://www.timeanddate.com/worldclock/sunearth.html
#    ?month=6&day=21&year=2030&hour=0&min=0&sec=0&n=&ntxt=&earth=0

@pytest.mark.parametrize('dt, true_lat, true_lon', [
    (datetime(2018, 9, 29, 0, 0), -(2 + 18/60), (177 + 37/60)),
    (datetime(2018, 9, 29, 14, 0), -(2 + 32/60), -(32 + 25/60)),
    (datetime(1992, 2, 14, 0, 0), -(13 + 20/60), -(176 + 26/60)),
    (datetime(2030, 6, 21, 0, 0), (23 + 26/60), -(179 + 34/60))
    ])
def test_solar_position(dt, true_lat, true_lon):
    lat, lon = _solar_position(dt)
    assert pytest.approx(true_lat, 0.1) == lat
    assert pytest.approx(true_lon, 0.1) == lon
