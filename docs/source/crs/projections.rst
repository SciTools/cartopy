.. _cartopy_projections:

Cartopy projection list
=======================


PlateCarree
-----------

.. autoclass:: cartopy.crs.PlateCarree

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='110m')
    ax.gridlines()



.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.coastlines(resolution='110m')
    ax.gridlines()


AlbersEqualArea
---------------

.. autoclass:: cartopy.crs.AlbersEqualArea

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(5.129856429268199, 3))
    ax = plt.axes(projection=ccrs.AlbersEqualArea())
    ax.coastlines(resolution='110m')
    ax.gridlines()


AzimuthalEquidistant
--------------------

.. autoclass:: cartopy.crs.AzimuthalEquidistant

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3, 3))
    ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_latitude=90))
    ax.coastlines(resolution='110m')
    ax.gridlines()


LambertConformal
----------------

.. autoclass:: cartopy.crs.LambertConformal

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(4.28969332204859, 3))
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.coastlines(resolution='110m')
    ax.gridlines()


LambertCylindrical
------------------

.. autoclass:: cartopy.crs.LambertCylindrical

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(9.42477796076938, 3))
    ax = plt.axes(projection=ccrs.LambertCylindrical())
    ax.coastlines(resolution='110m')
    ax.gridlines()


Mercator
--------

.. autoclass:: cartopy.crs.Mercator

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3.5090701847348473, 3))
    ax = plt.axes(projection=ccrs.Mercator())
    ax.coastlines(resolution='110m')
    ax.gridlines()


Miller
------

.. autoclass:: cartopy.crs.Miller

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(4.091529019548417, 3))
    ax = plt.axes(projection=ccrs.Miller())
    ax.coastlines(resolution='110m')
    ax.gridlines()


Mollweide
---------

.. autoclass:: cartopy.crs.Mollweide

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines(resolution='110m')
    ax.gridlines()


Orthographic
------------

.. autoclass:: cartopy.crs.Orthographic

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3, 3))
    ax = plt.axes(projection=ccrs.Orthographic())
    ax.coastlines(resolution='110m')
    ax.gridlines()


Robinson
--------

.. autoclass:: cartopy.crs.Robinson

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(5.914966076674721, 3))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.coastlines(resolution='110m')
    ax.gridlines()


Sinusoidal
----------

.. autoclass:: cartopy.crs.Sinusoidal

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6.0100710855457855, 3))
    ax = plt.axes(projection=ccrs.Sinusoidal())
    ax.coastlines(resolution='110m')
    ax.gridlines()


Stereographic
-------------

.. autoclass:: cartopy.crs.Stereographic

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3, 3))
    ax = plt.axes(projection=ccrs.Stereographic())
    ax.coastlines(resolution='110m')
    ax.gridlines()


TransverseMercator
------------------

.. autoclass:: cartopy.crs.TransverseMercator

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.TransverseMercator())
    ax.coastlines(resolution='110m')
    ax.gridlines()


UTM
---

.. autoclass:: cartopy.crs.UTM

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(0.12857142857142856, 3))
    ax = plt.axes(projection=ccrs.UTM(zone=30))
    ax.coastlines(resolution='110m')
    ax.gridlines()


InterruptedGoodeHomolosine
--------------------------

.. autoclass:: cartopy.crs.InterruptedGoodeHomolosine

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6.922806295266478, 3))
    ax = plt.axes(projection=ccrs.InterruptedGoodeHomolosine())
    ax.coastlines(resolution='110m')
    ax.gridlines()


RotatedPole
-----------

.. autoclass:: cartopy.crs.RotatedPole

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5))
    ax.coastlines(resolution='110m')
    ax.gridlines()


OSGB
----

.. autoclass:: cartopy.crs.OSGB

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(1.6153846153846154, 3))
    ax = plt.axes(projection=ccrs.OSGB())
    ax.coastlines(resolution='50m')
    ax.gridlines()


EuroPP
------

.. autoclass:: cartopy.crs.EuroPP

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(2.6153846153846154, 3))
    ax = plt.axes(projection=ccrs.EuroPP())
    ax.coastlines(resolution='50m')
    ax.gridlines()


Geostationary
-------------

.. autoclass:: cartopy.crs.Geostationary

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3, 3))
    ax = plt.axes(projection=ccrs.Geostationary())
    ax.coastlines(resolution='110m')
    ax.gridlines()


Gnomonic
--------

.. autoclass:: cartopy.crs.Gnomonic

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3, 3))
    ax = plt.axes(projection=ccrs.Gnomonic())
    ax.coastlines(resolution='110m')
    ax.gridlines()


LambertAzimuthalEqualArea
-------------------------

.. autoclass:: cartopy.crs.LambertAzimuthalEqualArea

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3.0065680601446605, 3))
    ax = plt.axes(projection=ccrs.LambertAzimuthalEqualArea())
    ax.coastlines(resolution='110m')
    ax.gridlines()


NorthPolarStereo
----------------

.. autoclass:: cartopy.crs.NorthPolarStereo

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3, 3))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.coastlines(resolution='110m')
    ax.gridlines()


OSNI
----

.. autoclass:: cartopy.crs.OSNI

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(2.4323374137753486, 3))
    ax = plt.axes(projection=ccrs.OSNI())
    ax.coastlines(resolution='10m')
    ax.gridlines()


SouthPolarStereo
----------------

.. autoclass:: cartopy.crs.SouthPolarStereo

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3, 3))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax.coastlines(resolution='110m')
    ax.gridlines()


