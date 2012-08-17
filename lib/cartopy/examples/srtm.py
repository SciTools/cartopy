"""
"""
import matplotlib.pyplot as plt
import numpy

import cartopy.crs as ccrs
import cartopy.io.srtm as io_srtm


def main():
    ax = plt.axes(projection=ccrs.PlateCarree())

    elev, crs, extent = io_srtm.srtm(-4, 50)

    elev = numpy.ma.masked_less_equal(elev, 0, copy=False)

    ax.imshow(numpy.ma.log(elev),
                extent=extent,
                transform=crs,
                cmap='Greens',
                )

#    ax.coastlines()
    plt.show()


from pylab import *
def set_shade(a,intensity=None,cmap=cm.jet,scale=10.0,azdeg=165.0,altdeg=45.0):
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

def hillshade(data,scale=10.0,azdeg=165.0,altdeg=45.0):
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
    az = azdeg*pi/180.0
    alt = altdeg*pi/180.0
    # gradient in x and y directions
    dx, dy = gradient(data/float(scale))
    slope = 0.5*pi - arctan(hypot(dx, dy))
    aspect = arctan2(dx, dy)
    intensity = sin(alt)*sin(slope) + cos(alt)*cos(slope)*cos(-az - aspect - 0.5*pi)
    intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
    return intensity




def main3():
    ax = plt.axes(projection=ccrs.PlateCarree())

    elev, crs, extent = io_srtm.srtm(-4, 50)

    elev = numpy.ma.masked_less_equal(elev, 0, copy=False)

    from matplotlib.colors import LightSource

    ls = LightSource(azdeg=90, altdeg=80,
                     hsv_min_val=0.6, hsv_min_sat=0.9,
                     hsv_max_val=0.8, hsv_max_sat=1,
                     )

#    rgb = ls.shade(elev, plt.get_cmap('Greens', 3))
    import matplotlib.colors as mcolors
    rgb = set_shade(elev,
                    cmap=mcolors.ListedColormap([plt.get_cmap('Greens', 3)(0.5)])
                    )


    ax.imshow(rgb,
                extent=extent,
                transform=crs
                )

    x = numpy.linspace(extent[0], extent[1], elev.shape[0])
    y = numpy.linspace(extent[2], extent[3], elev.shape[1])

    ax.contour(x, y, elev, 100,
               linestyles='-',
               colors='blue',
               linewidths=0.3,
               alpha=0.4,
               transform=crs,
               )

    plt.show()


def main2():
    ax = plt.axes(projection=ccrs.PlateCarree())

    elev, crs, extent = io_srtm.srtm_composite(-6, 50, 4, 3)

    elev = numpy.ma.masked_less_equal(elev, 0, copy=False)

    ax.imshow(numpy.ma.log(elev**2),
              extent=extent,
              transform=crs,
              cmap='Greens',
              )

#    ax.coastlines()

    plt.show()

if __name__ == '__main__':
    main2()
