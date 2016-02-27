"""
The Peirce Quincuncial projection

Original code in Tcl by Kevin B. Kenny, 2007
Code translated from Tcl to Python by Jonathan Feinberg, 2016
"""

import numpy as np
from scipy.special import ellipk, ellipkinc, ellipj

__author__ = "Jonathan Feinberg"
__email__ = "jonathan@feinberg.no"
__credits__ = ["Jonathan Feinberg", "Kevin B. Kenny"]
__all__ = ["forward_pq", "inverse_pq"]


def forward_ssc(lon, lat):
    """
Calculate the stereographic projection and Swartz-Cristoffel transformation.

Args:
    lat (np.ndarray) : Latitude of the point in radians.
    lon (np.ndarray) : Longtetude of the point in radians.

Returns:
    (np.ndarray, np.ndarray) : Coordinates on the [-1,1]^2 cube.
    """

    # Compute the auxiliary quantities 'm' and 'n'. Set 'm' to match
    # the sign of 'lon' and 'n' to be positive if |lon| > pi/2
    cos_phiosqrt2 = np.sqrt(2) / 2 * np.cos(lat)

    cos_a = cos_phiosqrt2 * (np.sin(lon) + np.cos(lon))
    cos_b = cos_phiosqrt2 * (np.sin(lon) - np.cos(lon))
    sin_a = np.sqrt(1.0 - cos_a * cos_a)
    sin_b = np.sqrt(1.0 - cos_b * cos_b)
    sin2_m = 1.0 + cos_a*cos_b - sin_a*sin_b
    sin2_n = 1.0 - cos_a*cos_b - sin_a*sin_b

    sin2_m = np.where(sin2_m < 0, 0, sin2_m)
    sin_m = np.sqrt(sin2_m)
    sin2_m = np.where(sin2_m > 1, 1, sin2_m)
    cos_m = np.sqrt(1.0 - sin2_m)

    sin_m *= np.sign(np.sin(lon))
    sin2_n *= sin2_n > 0

    sin_n = np.sqrt(sin2_n)
    sin2_n = np.where(sin2_n > 1, 1, sin2_n)
    cos_n = np.sqrt(1.0 - sin2_n)

    sin_n *= -np.sign(np.cos(lon))

    # Compute elliptic integrals to map the disc to the square
    xcorr = ellipkinc(np.arctan(sin_m/cos_m), 0.5)
    ycorr = ellipkinc(np.arctan(sin_n/cos_n), 0.5)

    return xcorr, ycorr


def inverse_ssc(xcorr, ycorr):
    """
Calculate the inverse stereographic projection and Swartz-Cristoffel
transformation.

Args:
    xcorr (array_like): x co-ordinates of a point on the map
    ycorr (array_like): y co-ordinates of a point on the map
    """
    # Compute the elliptic functions to map the plane onto the sphere
    cnx = ellipj(xcorr, 0.5)[1]
    cny = ellipj(ycorr, 0.5)[1]

    # Undo the mapping to latitude and longitude
    arc1 = np.arccos(-cnx**2)
    arc2 = np.arccos(cny**2)

    cos_a = np.cos(0.5*(arc1-arc2))
    cos_b = -np.cos(0.5*(arc1+arc2))

    lon = np.pi/4 - np.arctan2(cos_b, cos_a)
    lat = np.arccos(np.hypot(cos_b, cos_a))

    return lon, lat


def forward_pq(lon, lat, lon0=20.):
    """
Converts geodetic co-ordinates to the Peirce Quincuncial projection.

Args:
    lon (array_like) : Longitude of the point to be projected in degrees.
    lat (array_like) : Latitude of the point to be projected in degrees.

Kwargs:
    lon0 (float) : Longitude of the central meridian.  (Conventionally 20.)

Returns:
    (np.ndarray, np.ndarray) : Coordinates in Peirce Quincuncial projection

Examples:
    >>> lon = [-90, -90, 90, 90]
    >>> lat = [-45, 45, -45, 45]
    >>> x,y = forward_pq(lon, lat)
    >>> print(x)
    [-2.26997292 -0.74939046  2.26997292  0.74939046]
    >>> print(y)
    [-1.87266709 -0.35208463  1.87266709  0.35208463]
    """

    # numpy safety check
    lon = np.asfarray(lon)
    lat = np.asfarray(lat)

    if len(lon.shape) < len(lat.shape):
        lon = lon*np.ones(lat.shape)
    elif len(lon.shape) > len(lat.shape):
        lat = lat*np.ones(lon.shape)

    # rotate to angles of interest
    lon = lon - lon0 + 180
    lon = np.where((lon < 0)+(lon > 360),
                   lon - 360 * np.floor(lon/360), lon)

    # from degree to radians
    lon = (lon - 180) * np.pi / 180
    lat = lat * np.pi / 180

    # perform transformation
    xcorr, ycorr = forward_ssc(lon, lat)

    # Reflect the Southern Hemisphere outward
    neglat = lat < 0
    scale = ellipk(0.5)*2

    sector0 = neglat * (lon < -0.75*np.pi)
    negs = ~sector0
    sector1 = neglat * negs * (lon < -0.25*np.pi)
    negs *= ~sector1
    sector2 = neglat * negs * (lon < 0.25*np.pi)
    negs *= ~sector2
    sector3 = neglat * negs * (lon < 0.75*np.pi)
    negs *= ~sector3
    sector4 = neglat * negs * (lon >= 0.75*np.pi)

    ycorr = np.where(sector0, scale-ycorr, ycorr)
    xcorr = np.where(sector1, -scale-xcorr, xcorr)
    ycorr = np.where(sector2, -scale-ycorr, ycorr)
    xcorr = np.where(sector3, scale-xcorr, xcorr)
    ycorr = np.where(sector4, scale-ycorr, ycorr)

    # Rotate the square by 45 degrees to fit the screen better
    xcorr, ycorr = (xcorr - ycorr) * np.sqrt(2) / 2,\
        (xcorr + ycorr) * np.sqrt(2) / 2

    return xcorr, ycorr


def inverse_pq(xcorr, ycorr, lon0=20.):
    """
Converts Peirce Quincuncial map co-ordinates to latitude and longitude.

Args:
    xcorr (array_like): normalized x co-ordinates of a point on the map
    ycorr (array_like): normalized y co-ordinates of a point on the map

Kwargs:
    lon0 (array_like): Longitude of the center of projection

Returns:
    (np.ndarray, np.ndarray): Returns a list consisting of the longitude and
    latitude in degrees.

Examples:
    >>> xcorr = [-1, -1, 1, 1]
    >>> ycorr = [-1, 1, -1, 1]
    >>> lon, lat = inverse_pq(xcorr, ycorr)
    >>> print(lon)
    [ -70. -160.   20.  110.]
    >>> print(lat)
    [ 18.10370722  18.10370722  18.10370722  18.10370722]
    """
    # numpy safety check
    xcorr = np.asfarray(xcorr)
    ycorr = np.asfarray(ycorr)

    if len(xcorr.shape) < len(ycorr.shape):
        xcorr = xcorr*np.ones(ycorr.shape)
    elif len(xcorr.shape) > len(ycorr.shape):
        ycorr = ycorr*np.ones(xcorr.shape)

    # Rotate xcorr and ycorr 45 degrees
    xcorr, ycorr = (xcorr + ycorr) * np.sqrt(2) / 2,\
        (ycorr - xcorr) * np.sqrt(2) / 2

    # Reflect Southern Hemisphere into the Northern
    limit = ellipk(0.5)
    scale = ellipk(0.5)*2

    sector0 = xcorr < -limit
    sector1 = ~sector0 * (xcorr > limit)
    sector2 = ~sector1 * (ycorr < -limit)
    sector3 = ~sector2 * (ycorr > limit)

    xcorr = sector0 * (-scale - xcorr) + sector1 * (scale - xcorr) +\
        ~sector0 * ~sector1 * xcorr
    ycorr = sector2 * (-scale - ycorr) + sector3 * (scale - ycorr) +\
        ~sector2 * ~sector3 * ycorr
    southern = sector0 + sector1 + sector2 + sector3

    # Now we know that latitude will be positive.  If xcorr is negative, then
    # longitude will be negative; reflect the Western Hemisphere into the
    # Eastern.
    western = xcorr < 0
    xcorr *= (-1)**western

    # If ycorr is positive, the point is in the back hemisphere.  Reflect
    # it to the front.
    back = ycorr > 0
    ycorr *= (-1)**back

    # Finally, constrain longitude to be less than pi/4, by reflecting across
    # the 45 degree meridian.
    complement = xcorr > -ycorr
    xcorr, ycorr = np.where(complement, ycorr, xcorr),\
        np.where(complement, -xcorr, ycorr)

    # Transform
    lon, lat = inverse_ssc(xcorr, ycorr)

    # Undo the reflections that were done above, to get correct latitude
    # and longitude
    lon = np.where(complement, np.pi/2 - lon, lon)
    lon = np.where(back, np.pi - lon, lon)
    lon *= (-1)**western
    lat *= (-1)**southern

    # Convert latitude and longitude to degrees
    lon = lon * 180/np.pi + 180 + lon0
    lon = np.where((lon < 0)+(lon > 360), lon - 360*np.floor(lon/360.), lon)

    lon = lon - 180
    lat = lat * 180 / np.pi

    return lon, lat


if __name__ == "__main__":

    import doctest
    doctest.testmod()
