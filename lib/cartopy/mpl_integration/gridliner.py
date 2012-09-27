import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.collections as mcollections
import matplotlib.transforms as mtrans


from cartopy.crs import Projection, _RectangularProjection

import numpy

import matplotlib.ticker as mticker


degree_locator = mticker.MaxNLocator(nbins=9, steps=[1, 2, 3, 6, 15, 18])


class Gridliner(object):
    # NOTE: In future, one of these objects will be add-able to a GeoAxes (and maybe even a plain old mpl axes)
    # and it will call the "do_gridlines" method on draw. This will enable automatic gridline resolution determination
    # on zoom/pan.
    def __init__(self, axes):
        self.axes = axes
        self.gridlines = {}
        # might not be desirable for certain projections (osgb for instance)
        self.xlocator = degree_locator 
        self.ylocator = degree_locator

    def do_gridlines(self, ax, crs, nx=None, ny=None, background_patch=None, collection_kwargs=None):
        x_lim, y_lim = self.get_domain(ax, crs, nx=nx, ny=ny, background_patch=background_patch)

        if not isinstance(crs, mtrans.Transform):
            transform = crs._as_mpl_transform(ax)
        else:
            transform = crs
        
        rc_params = matplotlib.rcParams

        n_steps = 30

        lines = []
        
        # XXX this bit is cartopy specific. (for circular longitudes)
        if isinstance(crs, Projection) and isinstance(crs, _RectangularProjection) and numpy.diff(x_lim) == 2 * crs._half_width:
            endpoint = False
        else:
            endpoint = True
        
        x_ticks = self.xlocator.tick_values(x_lim[0], x_lim[1])
        y_ticks = self.ylocator.tick_values(y_lim[0], y_lim[1])
        
        for x in numpy.linspace(min(x_ticks), max(x_ticks), len(y_ticks), endpoint=endpoint):
            l = zip(numpy.zeros(n_steps) + x, 
                    numpy.linspace(min(y_ticks), max(y_ticks), n_steps, endpoint=True)
                    )
            lines.append(l)

        if collection_kwargs is None:
            collection_kwargs = {}
        collection_kwargs = collection_kwargs.copy()
        collection_kwargs['transform'] = transform
        # XXX doesn't gracefully handle lw vs linewidth aliases...
        collection_kwargs.setdefault('color', rc_params['grid.color'])
        collection_kwargs.setdefault('linestyle', rc_params['grid.linestyle'])
        collection_kwargs.setdefault('linewidth', rc_params['grid.linewidth'])
        
        x_lc = mcollections.LineCollection(lines, **collection_kwargs)
        self.gridlines['x'] = x_lc
        self.axes.add_collection(x_lc, autolim=False)


        lines = []
        
        for y in numpy.linspace(min(y_ticks), max(y_ticks), len(x_ticks), endpoint=True):
#            l = zip(numpy.linspace(x_lim[0], x_lim[1], n_steps, endpoint=True),
#                              numpy.zeros(n_steps) + y)
#            l = zip(x_ticks, numpy.zeros(len(x_ticks)) + y)
            l = zip(
                    numpy.linspace(min(x_ticks), max(x_ticks), n_steps, endpoint=True),
                    numpy.zeros(n_steps) + y,
                    )
            lines.append(l)

        y_lc = mcollections.LineCollection(lines, **collection_kwargs)
        self.gridlines['y'] = y_lc # <--- not sure about this....
        self.axes.add_collection(y_lc, autolim=False)
        return x_lc, y_lc

    def get_domain(self, ax, crs, nx=None, ny=None, background_patch=None):
        """Returns x_range, y_range"""
        
        DEBUG = False
        
        if not isinstance(crs, mtrans.Transform):
            transform = crs._as_mpl_transform(ax)
        else:
            transform = crs
            
        ax_transform = ax.transAxes
        
        desired_trans = ax.transAxes - transform

        nx = nx or 30
        ny = ny or 30
        x = numpy.linspace(1e-9, 1-1e-9, nx)
        y = numpy.linspace(1e-9, 1-1e-9, ny)
        x, y = numpy.meshgrid(x, y)

        coords = numpy.concatenate([x.flatten()[:, None], y.flatten()[:, None]], 1)

        in_data = desired_trans.transform(coords)

        ax_to_bkg_patch = self.axes.transAxes - background_patch.get_transform()

        ok = numpy.zeros(in_data.shape[:-1], dtype=numpy.bool)
        # XXX Vectorise contains_point
        for i, val in enumerate(in_data):
            # convert the coordinates of the data to the background patches coordinates
            background_coord = ax_to_bkg_patch.transform(coords[i:i+1, :])
            if background_patch.get_path().contains_point(background_coord[0, :]):
                color = 'r'
                ok[i] = True
            else:
                color = 'b'

            if DEBUG:
                import matplotlib.pyplot as plt
                plt.plot(coords[i, 0], coords[i, 1], 'o' + color, clip_on=False, transform=ax_transform)
#                plt.text(coords[i, 0], coords[i, 1], str(val), clip_on=False,
#                         transform=ax_transform, rotation=23,
#                         horizontalalignment='right')

        inside = in_data[ok, :]
        x_range = numpy.nanmin(inside[:, 0]), numpy.nanmax(inside[:, 0])
        y_range = numpy.nanmin(inside[:, 1]), numpy.nanmax(inside[:, 1])
        
        # XXX Cartopy specific thing. Perhaps make this bit a specialisation in a subclass...
        if isinstance(crs, Projection):
            x_range = numpy.clip(x_range, *crs.x_limits)
            y_range = numpy.clip(y_range, *crs.y_limits)
            
            # if the limit is >90 of the full x limit, then just use the full x limit (this makes circular handling better)
            prct = numpy.abs(numpy.diff(x_range) / numpy.diff(crs.x_limits))
            if prct > 0.9:
                x_range = crs.x_limits
                
        return x_range, y_range