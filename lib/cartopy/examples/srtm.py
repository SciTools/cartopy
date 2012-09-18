"""
"""
import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
import cartopy.io.srtm as io_srtm


def main():
    ax = plt.axes(projection=ccrs.PlateCarree())

    elev, crs, extent = io_srtm.srtm(-4, 50)

    elev = np.ma.masked_less_equal(elev, 0, copy=False)

    ax.imshow(np.ma.log(elev),
                extent=extent,
                transform=crs,
                cmap='Greens',
                )

#    ax.coastlines()
    plt.show()


def set_shade(a, intensity=None, cmap='jet', scale=10.0, azdeg=165.0, altdeg=45.0):
    """
    from http://rnovitsky.blogspot.co.uk/2010/04/using-hillshade-image-as-intensity.html

    sets shading for data array based on intensity layer
    or the data's value itself.
    inputs:
    a - a 2-d array or masked array
    intensity - a 2-d array of same size as a (no chack on that)
                      representing the intensity layer. if none is given
                      the data itself is used after getting the hillshade values
                      see hillshade for more details.
    cmap - a colormap (e.g matplotlib.colors.LinearSegmentedColormap
                instance)
    scale,azdeg,altdeg - parameters for hilshade function see there for
                more details
    output:
    rgb - an rgb set of the Pegtop soft light composition of the data and
             intensity can be used as input for imshow()
    based on ImageMagick's Pegtop_light:
    http://www.imagemagick.org/Usage/compose/#pegtoplight
    """
    if intensity is None:
        # hilshading the data
        intensity = hillshade(a,scale=10.0,azdeg=165.0,altdeg=45.0)
    else:
        # or normalize the intensity
        intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())

    # get rgb of normalized data based on cmap
    rgb = cmap((a-a.min())/float(a.max()-a.min()))[:,:,:3]
    # form an rgb eqvivalent of intensity
    d = intensity.repeat(3).reshape(rgb.shape)
    # simulate illumination based on pegtop algorithm.
    rgb = 2*d*rgb+(rgb**2)*(1-2*d)
    return rgb


def hillshade(data, scale=10.0, azdeg=165.0, altdeg=45.0):
    ''' convert data to hillshade based on matplotlib.colors.LightSource class.
    input:
         data - a 2-d array of data
         scale - scaling value of the data. higher number = lower gradient
         azdeg - where the light comes from: 0 south ; 90 east ; 180 north ;
                      270 west
         altdeg - where the light comes from: 0 horison ; 90 zenith
    output: a 2-d array of normalized hilshade
    '''
    # convert alt, az to radians
    az = azdeg * np.pi/180.0
    alt = altdeg * np.pi/180.0
    # gradient in x and y directions
    dx, dy = np.gradient(data / float(scale))
    slope = 0.5 * np.pi - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(dx, dy)
    intensity = np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * np.cos(-az - aspect - 0.5 * np.pi)
    intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
    return intensity


def main3():
    ax = plt.axes(projection=ccrs.PlateCarree())

    elev, crs, extent = io_srtm.srtm_composite(-5, 52, 2, 2)

    elev = np.ma.masked_less_equal(elev, 0, copy=False)

    use_mpl_light_source = False
    
    if use_mpl_light_source:
        from matplotlib.colors import LightSource
    
        ls = LightSource(azdeg=90, altdeg=80,
                         hsv_min_val=0.6, hsv_min_sat=0.9,
                         hsv_max_val=0.8, hsv_max_sat=1,
                         )
    
        rgb = ls.shade(elev, plt.get_cmap('Greens', 3))
    else:
        import matplotlib.colors as mcolors
        rgb = set_shade(elev,
                        cmap=mcolors.ListedColormap([plt.get_cmap('Greens', 3)(0.5)])
                    )


    ax.imshow(rgb,
                extent=extent,
                transform=crs
                )

    x = np.linspace(extent[0], extent[1], elev.shape[0])
    y = np.linspace(extent[2], extent[3], elev.shape[1])
#
#    ax.contour(x, y, elev, 100,
#               linestyles='-',
#               colors='blue',
#               linewidths=0.3,
#               alpha=0.4,
#               transform=crs,
#               )

    plt.show()


def main2():
    ax = plt.axes(projection=ccrs.PlateCarree())

    elev, crs, extent = io_srtm.srtm_composite(-6, 50, 4, 3)

    elev = np.ma.masked_less_equal(elev, 0, copy=False)

    ax.imshow(np.ma.log(elev**2),
              extent=extent,
              transform=crs,
              cmap='Greens',
              )

    plt.show()


if __name__ == '__main__':
    main3()
