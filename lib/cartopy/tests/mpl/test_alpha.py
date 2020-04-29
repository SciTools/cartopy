from __future__ import (absolute_import, division, print_function)

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs


def test_alpha():

    # tests that both image and alpha arrays are warped

    plt_crs = ccrs.Geostationary(central_longitude=-155.)

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1, projection=plt_crs)

    latlon_crs = ccrs.PlateCarree()

    coords = [-162., -148., 17.5, 23.]

    ax.set_extent(coords, crs=latlon_crs)

    fake_data = np.zeros([100, 100])

    fake_alphas = np.zeros(fake_data.shape)
    
    image = ax.imshow(fake_data, extent=coords, transform=latlon_crs,
                      alpha=fake_alphas)

    plt.close()

    image_data = image.get_array()
    image_alpha = image.get_alpha()

    assert image_data.shape == image_alpha.shape
