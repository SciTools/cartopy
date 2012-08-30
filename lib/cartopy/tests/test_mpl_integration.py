import numpy as np

from matplotlib.testing.decorators import image_comparison as mpl_image_comparison
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def image_comparison(baseline_images=None, extensions=('png', ), tol=1e-3):
    # changes the mpl default to only use PNGs
    return mpl_image_comparison(baseline_images, extensions, tol)


@image_comparison(baseline_images=['global_map'])
def test_global_map():
    ax = plt.axes(projection=ccrs.Robinson())
#    ax.coastlines()
#    ax.gridlines(5)

    plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
    
    plt.plot([-0.08, 132], [51.53, 43.17], color='red', 
             transform=ccrs.PlateCarree())
    
    plt.plot([-0.08, 132], [51.53, 43.17], color='blue', 
             transform=ccrs.Geodetic()) 
    

if __name__=='__main__':
    import nose
    nose.runmodule(argv=['-s','--with-doctest'], exit=False)
