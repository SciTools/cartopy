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

import itertools
import os

import numpy as np

import cartopy.crs as ccrs


def find_projections():
    for obj_name, o in vars(ccrs).copy().items():
#        o = getattr(ccrs, obj_name)
        if (isinstance(o, type) and issubclass(o, ccrs.Projection) and
            not obj_name.startswith('_') and obj_name not in ['Projection']):

            # yield the projections
            yield o

def projection_rst(projection_cls):
    name = projection_cls.__name__
    print(name)


SPECIAL_CASES = {ccrs.PlateCarree: ['PlateCarree()', 'PlateCarree(central_longitude=180)'],
                 ccrs.RotatedPole: ['RotatedPole(pole_longitude=177.5, pole_latitude=37.5)'],
                 }


COASTLINE_RESOLUTION = {ccrs.OSNI: '10m',
                        ccrs.OSGB: '50m',
                        ccrs.EuroPP: '50m'}

PRJ_SORT_ORDER = {'PlateCarree': 1, 'Mercator': 2, 'Mollweide': 2, 'Robinson': 2,
                  'TransverseMercator': 2, 'LambertCylindrical': 2,
                  'LambertConformal': 2, 'Stereographic': 2, 'Miller': 2,
                  'Orthographic': 2, 'InterruptedGoodeHomolosine': 3,
                  'RotatedPole': 3, 'OSGB': 4}


groups = [('cylindrical', [ccrs.PlateCarree, ccrs.Mercator, ccrs.TransverseMercator,
                           ccrs.OSGB, ccrs.LambertCylindrical, ccrs.Miller, ccrs.RotatedPole]),
          ('pseudo-cylindrical', [ccrs.Mollweide, ccrs.Robinson]),
#          ('conic', [ccrs.aed]),
          ('azimuthal', [ccrs.Stereographic, ccrs.NorthPolarStereo,
                         ccrs.SouthPolarStereo, ccrs.Gnomonic, ccrs.Orthographic
                         ]),
          ('misc', [ccrs.InterruptedGoodeHomolosine]),
          ]


all_projections_in_groups = list(itertools.chain.from_iterable([g[1] for g in groups]))


if __name__ == '__main__':
    fname = os.path.join(os.path.dirname(__file__), 'source',
                         'crs', 'projections.rst')
    table = open(fname, 'w')

    table.write('.. _cartopy_projections:\n\n')
    table.write('Cartopy projection list\n')
    table.write('=======================\n\n\n')

    prj_class_sorter = lambda cls: (PRJ_SORT_ORDER.get(cls.__name__, []), cls.__name__)
    for prj in sorted(find_projections(), key=prj_class_sorter):
        name = prj.__name__
#        print prj in SPECIAL_CASES, prj in all_projections_in_groups, prj

        # put the class documentation on the left, and a sidebar on the right.

        aspect = (np.diff(prj().x_limits) / np.diff(prj().y_limits))[0]
        width = 3 * aspect
        if width == int(width):
            width = int(width)

        table.write(name + '\n')
        table.write('-' * len(name) + '\n\n')

        table.write('.. autoclass:: cartopy.crs.%s\n' % name)

#        table.write('Ipsum lorum....')

#        table.write("""\n\n
#
#.. sidebar:: Example
#""")

        for instance_creation_code in SPECIAL_CASES.get(prj, ['%s()' % name]):
            code = """
.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=({width}, 3))
    ax = plt.axes(projection=ccrs.{proj_constructor})
    ax.coastlines(resolution={coastline_resolution!r})
    ax.gridlines()

\n""".format(width=width, proj_constructor=instance_creation_code,
             coastline_resolution=COASTLINE_RESOLUTION.get(prj, '110m'))

            table.write(code)
