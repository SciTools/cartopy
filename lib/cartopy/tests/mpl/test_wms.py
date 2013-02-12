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

from nose.tools import assert_equal
from matplotlib.testing.decorators import image_comparison as mpl_image_comparison
import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
from cartopy.io.wms import WMS
from cartopy.tests.mpl import ImageTesting

import os
import os.path
import warnings


class MockLineString(object):
    """
    Mimicks the target_domain sent to WmsFactory.image_for_domain.
    
    """
    geom_type = "LineString"
    def __init__(self, bounds):
        self.bounds = bounds


def create_wms_plot():
    """
    Construct a plot which requests a WMS image but do not draw.
    
    Specifies a publicly hosted server but doesn't trigger the WMS request.
    
    """         
    server = "http://129.206.228.72/cached/osm?"
    layers = "osm_auto:all"
    origin = 'lower'
    size=(256, 256)
    plt.figure(figsize=(12,10))
    ax = plt.subplot(1,1,1, projection=ccrs.PlateCarree())
    ax.add_image(WMS(server, layers, origin=origin, size=size))


def test_wms_retrieve():
    # This tests the retrieval of a wms image, not the image itself.
    # Display of a mock wms image is tested in tests.mpl.test_wms_display.
    # TODO: Can we mock a wms server? This is currently an intergration test.

    # Construct a WmsFactory.
    create_wms_plot()
    
    # Find the factory we just created.
    for item in plt.gca().img_factories:
        if isinstance(item[0], WMS):
            wms_factory = item[0]
            break
    
    # Simulate draw-time, when the actual WMS retrieval happens.
    target_domain = MockLineString([-180, -90, 180, 90])
    img, bbox, origin = wms_factory.image_for_domain(target_domain, None, srs=ccrs.PlateCarree())
    
    # We've retrieved an image. 
    assert_equal(img.format, "JPEG")
    assert_equal(img.size, (256, 256))
    
    plt.close()
    
    
def image_comparison(folder, filename):
    # Based on the ImageTesting decorator. 
    root_image_results = os.path.dirname(__file__) 
    expected_fname = os.path.join(root_image_results,
                                  'baseline_images', 'mpl', folder, filename)
    result_fname = os.path.join(root_image_results,
                                'output', folder, filename)

    if not os.path.isdir(os.path.dirname(expected_fname)):
        os.makedirs(os.path.dirname(expected_fname))
    if not os.path.isdir(os.path.dirname(result_fname)):
        os.makedirs(os.path.dirname(result_fname))

    if not os.path.exists(expected_fname):
        warnings.warn("Creating baseline image: {}".format(expected_fname))
        plt.savefig(expected_fname)

    plt.savefig(result_fname)
    
    import matplotlib.testing.compare as mcompare
    err = mcompare.compare_images(expected_fname, result_fname, tol=1e-3)
    
    if err is not None:
        raise Exception(err)
    os.remove(result_fname)


def striped_image(size):
    """Create a striped image of the given size."""
    from PIL import Image
    import numpy as np
    pels = np.zeros((size[1], size[0], 3), dtype=np.int8)
    height = size[1]
    quarter_w = size[0]/4
    # 4 columns of different colour.
    for i in range(height):
        v = (i*256)/height
        pels[i, 0:quarter_w, 0] = v  # 1st column red
        pels[i, quarter_w:quarter_w*2, 1] = v  # 2nd column green
        pels[i, quarter_w*2:quarter_w*3, 2] = v  # 3rd column blue
    return Image.frombuffer("RGB", size, pels, 'raw', "RGB", 0, 1)


def test_wms_integration():
    # Test user-style WMS plotting.
    # Uses a mock urllib to elliminate dependency on an external server.

    class MockUrllib(object):
        """A mock for urllib that just returns a striped jpeg."""
        @classmethod
        def urlopen(cls, *args, **kwargs):
            """Return a file-like object containing jpeg bytes"""
            import StringIO
            s = StringIO.StringIO()
            img = striped_image((256, 256))
            img.save(s, "JPEG")
            s.seek(0)
            return s

    # Replace urllib with our mock
    import cartopy.io.wms
    original_urllib = cartopy.io.wms.urllib
    cartopy.io.wms.urllib = MockUrllib
    
    server = "http://monty-python.com/wms?"
    layers = "camelot"
    origin = 'upper'

    ax = plt.subplot(2,2,1, projection=ccrs.PlateCarree())
    ax.add_image(WMS(server, layers, origin=origin))

    ax = plt.subplot(2,2,2, projection=ccrs.PlateCarree(), xlim=[-20, 20], ylim=[30, 70])
    ax.add_image(WMS(server, layers, origin=origin))

    ax = plt.subplot(2,2,3, projection=ccrs.PlateCarree(), xlim=[-13, 3], ylim=[50, 60])
    ax.add_image(WMS(server, layers, origin=origin))

    ax = plt.subplot(2,2,4, projection=ccrs.PlateCarree(), xlim=[-7, -3], ylim=[49, 53])
    ax.add_image(WMS(server, layers, origin=origin))


    # We can't use the ImageTesting decorator because we need to clean up after image testing.
    image_comparison("test_wms", "wms_integration.png")

    # Restore normal urllib operations.
    cartopy.io.wms.urllib = original_urllib


if __name__ == '__main__':
    
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)

    

