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
import numpy as np


def julian_day(date):
    """
    Calculate the Julian day from an input datetime.

    Parameters
    ----------
    date
        A UTC datetime object.

    Note
    ----
    Algorithm implemented following equations from Chapter 3 (Algorithm 14):
    Vallado, David 'Fundamentals of Astrodynamics and Applications', (2007)

    Julian day epoch is: noon on January 1, 4713 BC (proleptic Julian)
                         noon on November 24, 4714 BC (proleptic Gregorian)

    """
    year = date.year
    month = date.month
    day = date.day
    hour = date.hour
    minute = date.minute
    second = date.second

    # January/February correspond to months 13/14 respectively
    # for the constants to work out properly
    if month < 3:
        month += 12
        year -= 1

    B = 2 - int(year/100) + int(int(year/100)/4)
    C = ((second/60 + minute)/60 + hour)/24

    JD = (int(365.25*(year + 4716)) + int(30.6001*(month+1)) +
          day + B - 1524.5 + C)
    return JD


def solar_position(date):
    """
    Calculate the latitude and longitude point where the sun is
    directly overhead for the given date.

    Parameters
    ----------
    date
        A UTC datetime object.

    Returns
    -------
    (latitude, longitude) in degrees

    Note
    ----
    Algorithm implemented following equations from Chapter 5 (Algorithm 29):
    Vallado, David 'Fundamentals of Astrodynamics and Applications', (2007)

    """
    # NOTE: Constants are in degrees in the textbook,
    #       so we need to convert the values from deg2rad when taking sin/cos

    # Centuries from J2000
    T_UT1 = (julian_day(date) - 2451545.0)/36525

    # solar longitude (deg)
    lambda_M_sun = (280.460 + 36000.771*T_UT1) % 360

    # solar anomaly (deg)
    M_sun = (357.5277233 + 35999.05034*T_UT1) % 360

    # ecliptic longitude
    lambda_ecliptic = (lambda_M_sun + 1.914666471*np.sin(np.deg2rad(M_sun)) +
                       0.019994643*np.sin(np.deg2rad(2*M_sun)))

    # obliquity of the ecliptic (epsilon in Vallado's notation)
    epsilon = 23.439291 - 0.0130042*T_UT1

    # declination of the sun
    delta_sun = np.rad2deg(np.arcsin(np.sin(np.deg2rad(epsilon)) *
                                     np.sin(np.deg2rad(lambda_ecliptic))))

    # Greenwich mean sidereal time (seconds)
    theta_GMST = (67310.54841 +
                  (876600*3600 + 8640184.812866)*T_UT1 +
                  0.093104*T_UT1**2 -
                  6.2e-6*T_UT1**3)
    # Convert to degrees
    theta_GMST = (theta_GMST % 86400)/240

    # Right ascension calculations
    numerator = (np.cos(np.deg2rad(epsilon)) *
                 np.sin(np.deg2rad(lambda_ecliptic)) /
                 np.cos(np.deg2rad(delta_sun)))
    denominator = (np.cos(np.deg2rad(lambda_ecliptic)) /
                   np.cos(np.deg2rad(delta_sun)))

    alpha_sun = np.rad2deg(np.arctan2(numerator, denominator))

    # longitude is opposite of Greenwich Hour Angle (GHA)
    # GHA == theta_GMST - alpha_sun
    lon = -(theta_GMST-alpha_sun)
    if lon < -180:
        lon += 360

    return (delta_sun, lon)
