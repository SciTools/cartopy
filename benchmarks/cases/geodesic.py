# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import numpy as np

from cartopy.geodesic import Geodesic


class DirectSuite:
    params = [
        ('scalar', 'array'),
        ('scalar', 'array'),
        ('zero', 'scalar', 'array'),
    ]
    param_names = ['points', 'azimuths', 'distances']

    def setup(self, points, azimuths, distances):
        self.geod = Geodesic()

        self.points = {
            'scalar': (0, 0),
            'array': np.arange(200).reshape(-1, 2),
        }[points]

        self.azimuths = {
            'scalar': 0,
            'array': np.arange(100) * 2,
        }[azimuths]

        self.distances = {
            'zero': 0,
            'scalar': 100,
            'array': np.logspace(2, 5, 100),
        }[distances]

    def time_direct(self, points, azimuths, distances):
        self.geod.direct(self.points, self.azimuths, self.distances)


class InverseSuite:
    params = [
        ('scalar', 'array'),
        ('scalar', 'array'),
    ]
    param_names = ['points', 'endpoints']

    def setup(self, points, endpoints):
        self.geod = Geodesic()

        possible_points = {
            'scalar': (0, 0),
            'array': np.arange(200).reshape(-1, 2),
        }

        self.points = possible_points[points]
        self.endpoints = possible_points[endpoints]

    def time_inverse(self, points, endpoints):
        self.geod.inverse(self.points, self.endpoints)


class CircleSuite:
    params = [100, 1000, 10000, 100e3, 1000e3, 10000e3]  # in metres
    param_names = ['radius']

    def time_circle(self, radius):
        geod = Geodesic()
        geod.circle(0, 0, radius)
