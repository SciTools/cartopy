Cartopy projection list
=======================


RotatedPole
-----------

:class:`~cartopy.crs.RotatedPole`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.RotatedPole(pole_longitude=177.5, pole_latitude=37.5))
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


Mercator
--------

:class:`~cartopy.crs.Mercator`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.Mercator())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


LambertCylindrical
------------------

:class:`~cartopy.crs.LambertCylindrical`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(9.42477796077, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.LambertCylindrical())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


OSGB
----

:class:`~cartopy.crs.OSGB`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(1.61538461538, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.OSGB())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


Stereographic
-------------

:class:`~cartopy.crs.Stereographic`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3.0, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.Stereographic())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


Gnomonic
--------

:class:`~cartopy.crs.Gnomonic`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3.0, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.Gnomonic())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


PlateCarree
-----------

:class:`~cartopy.crs.PlateCarree`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.PlateCarree())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()



.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.PlateCarree(central_longitude=180))
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


Mollweide
---------

:class:`~cartopy.crs.Mollweide`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6.0, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.Mollweide())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


InterruptedGoodeHomolosine
--------------------------

:class:`~cartopy.crs.InterruptedGoodeHomolosine`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6.92280629527, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.InterruptedGoodeHomolosine())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


Miller
------

:class:`~cartopy.crs.Miller`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(4.09152901955, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.Miller())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


SouthPolarStereo
----------------

:class:`~cartopy.crs.SouthPolarStereo`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3.0, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.SouthPolarStereo())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


Orthographic
------------

:class:`~cartopy.crs.Orthographic`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3.0, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.Orthographic())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


TransverseMercator
------------------

:class:`~cartopy.crs.TransverseMercator`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(6, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.TransverseMercator())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


NorthPolarStereo
----------------

:class:`~cartopy.crs.NorthPolarStereo`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(3.0, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.NorthPolarStereo())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()


Robinson
--------

:class:`~cartopy.crs.Robinson`

.. plot::

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    plt.figure(figsize=(5.91496652704, 3))
    delta = 0.125
    ax = plt.axes([0+delta, 0+delta, 1-delta, 1-delta], projection=ccrs.Robinson())
    #ax.set_global()
    ax.coastlines()
    ax.gridlines()
