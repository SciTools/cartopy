"""
The Peirce Quincuncial projection

Original code in Tcl by Kevin B. Kenny, 2007
Code translated from Tcl to Python by Jonathan Feinberg, 2016
"""

__author__ = "Jonathan Feinberg"
__email__ = "jonathan@feinberg.no"
__credits__ = ["Jonathan Feinberg", "Kevin B. Kenny"]

import numpy as np
from scipy.special import ellipk, ellipkinc, ellipj


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
    [ 4.49472464 -0.74939046  2.26997292  0.74939046]
    >>> print(y)
    [ 4.89203048 -0.35208463  1.87266709  0.35208463]
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
    lon = np.where((lon < 0)+(lon>360),
            lon - 360 * np.floor(lon/360), lon)

    # from degree to radians
    lon = (lon - 180) * np.pi / 180
    lat = lat * np.pi / 180

    # Compute the auxiliary quantities 'm' and 'n'. Set 'm' to match
    # the sign of 'lon' and 'n' to be positive if |lon| > pi/2
    cos_phiosqrt2 = np.sqrt(2) / 2 * np.cos(lat)
    cos_lon = np.cos(lon)
    sin_lon = np.sin(lon)

    cos_a = cos_phiosqrt2 * (sin_lon + cos_lon)
    cos_b = cos_phiosqrt2 * (sin_lon - cos_lon)
    sin_a = np.sqrt(1.0 - cos_a * cos_a)
    sin_b = np.sqrt(1.0 - cos_b * cos_b)
    sin2_m = 1.0 + cos_a*cos_b - sin_a*sin_b
    sin2_n = 1.0 - cos_a*cos_b - sin_a*sin_b

    sin2_m = np.where(sin2_m < 0, 0, sin2_m)
    sin_m = np.sqrt(sin2_m)
    sin2_m = np.where(sin2_m > 1, 1, sin2_m)
    cos_m = np.sqrt(1.0 - sin2_m)

    sin_m *= np.sign(sin_lon)
    sin2_n *= sin2_n > 0

    sin_n = np.sqrt(sin2_n)
    sin2_n = np.where(sin2_n > 1, 1, sin2_n)
    cos_n = np.sqrt(1.0 - sin2_n)

    sin_n *= -np.sign(cos_lon)

    # Compute elliptic integrals to map the disc to the square
    x = ellipkinc(np.arctan(sin_m/cos_m), 0.5)
    y = ellipkinc(np.arctan(sin_n/cos_n), 0.5)


    # Reflect the Southern Hemisphere outward
    neglat = lat < 0
    scale = ellipk(0.5)*2

    s0 = neglat * (lon < -0.75*np.pi)
    s1 = neglat * ~s0 * (lon < -0.25*np.pi)
    s2 = neglat * ~s1 * (lon < 0.25*np.pi)
    s3 = neglat * ~s2 * (lon < 0.75*np.pi)
    s4 = neglat * (lon >= 0.75*np.pi)

    y = np.where(s0, scale-y, y)
    x = np.where(s1, -scale-x, x)
    y = np.where(s2, -scale-y, y)
    x = np.where(s3, scale-x, x)
    y = np.where(s4, scale-y, y)

    # Rotate the square by 45 degrees to fit the screen better
    X = (x - y) * np.sqrt(2) / 2
    Y = (x + y) * np.sqrt(2) / 2

    return X, Y



def inverse_pq(x, y, lon0=20.):
    """
Converts Peirce Quincuncial map co-ordinates to latitude and longitude.

Args:
    x (array_like): normalized x co-ordinates of a point on the map
    y (array_like): normalized y co-ordinates of a point on the map

Kwargs:
    lon0 (array_like): Longitude of the center of projection

Returns:
    (np.ndarray, np.ndarray): Returns a list consisting of the longitude and latitude in degrees.

Examples:
    >>> x = [-1, -1, 1, 1]
    >>> y = [-1, 1, -1, 1]
    >>> lon, lat = inverse_pq(x, y)
    >>> print(lon)
    [ -70. -160.   20.  110.]
    >>> print(lat)
    [ 18.10370722  18.10370722  18.10370722  18.10370722]
    """
    # numpy safety check
    x = np.asfarray(x)
    y = np.asfarray(y)

    if len(x.shape) < len(y.shape):
        x = x*np.ones(y.shape)
    elif len(x.shape) > len(y.shape):
        y = y*np.ones(x.shape)

    # Rotate x and y 45 degrees
    X = (x + y) * np.sqrt(2) / 2
    Y = (y - x) * np.sqrt(2) / 2

    # Reflect Southern Hemisphere into the Northern
    limit = ellipk(0.5)
    scale = ellipk(0.5)*2

    s0 = X < -limit
    s1 = ~s0 * (X > limit)
    s2 = ~s1 * (Y < -limit)
    s3 = ~s2 * (Y > limit)

    X = s0 * (-scale - X) + s1 * (scale - X) + ~s0 * ~s1 * X
    Y = s2 * (-scale - Y) + s3 * (scale - Y) + ~s2 * ~s3 * Y
    southern = s0 + s1 + s2 + s3

    # Now we know that latitude will be positive.  If X is negative, then
    # longitude will be negative; reflect the Western Hemisphere into the
    # Eastern.
    western = X < 0
    X *= (-1)**western

    # If Y is positive, the point is in the back hemisphere.  Reflect
    # it to the front.
    back = Y > 0
    Y *= (-1)**back

    # Finally, constrain longitude to be less than pi/4, by reflecting across
    # the 45 degree meridian.
    complement = X > -Y
    X, Y = np.where(complement, Y, X), np.where(complement, -X, Y)

    # Compute the elliptic functions to map the plane onto the sphere
    cnx = ellipj(X, 0.5)[1]
    cny = ellipj(Y, 0.5)[1]

    # Undo the mapping to latitude and longitude
    a1 = np.arccos(-cnx**2)
    a2 = np.arccos(cny**2)
    b = 0.5*(a1+a2)
    a = 0.5*(a1-a2)

    cos_a = np.cos(a)
    cos_b = -np.cos(b)

    lon = np.pi/4 - np.arctan2(cos_b, cos_a)
    lat = np.arccos(np.hypot(cos_b, cos_a))

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

    # import matplotlib.pyplot as plt
    #
    # r = np.linspace(-180, 180, 200)[1:-1]
    # s = np.linspace(-90, 90, 200)[1:-1]
    # r,s = np.meshgrid(r, s)
    # z = r+s
    #
    # plt.contourf(r, s, z, 50)
    # plt.show()
    # plt.clf()
    #
    # R = forward_pq(r, s)
    #
    # plt.contourf(R[0], R[1], z, 50)
    # plt.show()
    # plt.clf()
    #
    # S = inverse_pq(R[0], R[1])
    #
    # plt.contourf(S[0], S[1], z, 50)
    # plt.show()
    #
