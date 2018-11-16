# (C) British Crown Copyright 2014 - 2018, Met Office
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
"""
Define the SlippyImageArtist class, which interfaces with
:class:`cartopy.io.RasterSource` instances at draw time, for interactive
dragging and zooming of raster data.

"""

from __future__ import (absolute_import, division, print_function)

from matplotlib.image import AxesImage
import matplotlib.artist


class SlippyImageArtist(AxesImage):

    """
    A subclass of :class:`~matplotlib.image.AxesImage` which provides an
    interface for getting a raster from the given object with interactive
    slippy map type functionality.

    Kwargs are passed to the AxesImage constructor.

    """
    def __init__(self, ax, raster_source, **kwargs):
        self.raster_source = raster_source
        super(SlippyImageArtist, self).__init__(ax, **kwargs)
        self.set_clip_path(ax.outline_patch)
        self.cache = []

        ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        ax.figure.canvas.mpl_connect('button_release_event', self.on_release)

        self.on_release()

    def on_press(self, event=None):
        self.user_is_interacting = True

    def on_release(self, event=None):
        self.user_is_interacting = False
        self.stale = True

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer, *args, **kwargs):
        if not self.get_visible():
            return

        ax = self.axes
        window_extent = ax.get_window_extent()
        [x1, y1], [x2, y2] = ax.viewLim.get_points()
        if not self.user_is_interacting:
            located_images = self.raster_source.fetch_raster(
                ax.projection, extent=[x1, x2, y1, y2],
                target_resolution=(window_extent.width, window_extent.height))
            self.cache = located_images

        for img, extent in self.cache:
            self.set_array(img)
            with ax.hold_limits():
                self.set_extent(extent)
            super(SlippyImageArtist, self).draw(renderer, *args, **kwargs)
