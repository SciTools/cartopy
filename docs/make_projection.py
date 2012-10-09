import itertools
import os

import numpy as np

import cartopy.crs as ccrs


def find_projections():
    for obj_name, o in vars(ccrs).iteritems():
#        o = getattr(ccrs, obj_name)
        if (isinstance(o, type) and issubclass(o, ccrs.Projection) and
            not obj_name.startswith('_') and obj_name not in ['Projection']):

            # yield the projections
            yield o

def projection_rst(projection_cls):
    name = projection_cls.__name__
    print name


SPECIAL_CASES = {ccrs.PlateCarree: ['PlateCarree()', 'PlateCarree(central_longitude=180)'],
                 ccrs.RotatedPole: ['RotatedPole(pole_longitude=177.5, pole_latitude=37.5)'],
                 }


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
                         'projections', 'table.rst')
    table = open(fname, 'w')

    table.write('Cartopy projection list\n')
    table.write('=======================\n\n\n')


    for prj in find_projections():
        name = prj.__name__
#        print prj in SPECIAL_CASES, prj in all_projections_in_groups, prj

        # put the class documentation on the left, and a sidebar on the right.

        aspect = (np.diff(prj().x_limits) / np.diff(prj().y_limits))[0]
        width = 3 * aspect


        table.write(name + '\n')
        table.write('-' * len(name) + '\n\n')

        table.write(':class:`~cartopy.crs.%s`\n' % name)

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

    plt.figure(figsize=(%s, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.%s)
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()

\n""" % (width, instance_creation_code)

            table.write(code)
