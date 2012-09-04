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
    
    
@image_comparison(baseline_images=['multiple_projections1'])
def test_multiple_projections():
    
    projections = [ccrs.PlateCarree(), 
                   ccrs.Robinson(), 
                   ccrs.RotatedPole(pole_latitude=45, pole_longitude=180),
                   ccrs.OSGB(),
                   ccrs.TransverseMercator(),
                   ccrs.Mercator(),
                   ccrs.LambertCylindrical(),
                   ccrs.Miller(),
                   ccrs.Gnomonic(),
                   ccrs.Stereographic(),
                   ccrs.NorthPolarStereo(),
                   ccrs.SouthPolarStereo(),
                   ccrs.Orthographic(),
                   ccrs.Mollweide(),
#                   ccrs.InterruptedGoodeHomolosine(),
                   ]
    
    fig = plt.figure(figsize=(10, 10))
    for i, prj in enumerate(projections, 1):
        ax = fig.add_subplot(5, 5, i, projection=prj)
    
        ax.set_global()
        
        ax.set_title(prj.__class__.__name__)
        
        ax.coastlines()
            
        plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
        
        plt.plot([-0.08, 132], [51.53, 43.17], color='red', 
                 transform=ccrs.PlateCarree())
        
        plt.plot([-0.08, 132], [51.53, 43.17], color='blue', 
                 transform=ccrs.Geodetic())

#
#@image_comparison(baseline_images=['image_transform'])
#def test_image_transforms():
#    plt.subplot(131, projection=ccrs.PlateCarree())
#    

if __name__=='__main__':
    import nose
    nose.runmodule(argv=['-s','--with-doctest'], exit=False)
