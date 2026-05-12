# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

from matplotlib.contour import QuadContourSet


class GeoContourSet(QuadContourSet):
    """
    A contourset designed to handle things like contour labels.

    """
    # nb. No __init__ method here - most of the time a GeoContourSet will
    # come from GeoAxes.contour[f]. These methods morph a ContourSet by
    # fiddling with instance.__class__.

    def clabel(self, *args, **kwargs):
        # Where contour paths exist at the edge of the globe, sometimes a
        # complete path in data space will become multiple paths when
        # transformed into axes or screen space.  Matplotlib's contour
        # labelling does not account for this so we need to give it the
        # pre-transformed paths to work with.

        # Define the transform that will take us from collection
        # coordinates through to axes projection coordinates.
        data_t = self.axes.transData
        col_to_data = self.get_transform() - data_t

        # Now that we have the transform, project all of this
        # collection's paths.
        paths = self.get_paths()
        new_paths = [col_to_data.transform_path(path) for path in paths]
        self.set_paths(new_paths)

        # The collection will now be referenced in axes projection
        # coordinates.
        self.set_transform(data_t)

        # Now that we have prepared the collection paths, call on
        # through to the underlying implementation.
        return super().clabel(*args, **kwargs)
