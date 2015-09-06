.. _map_styling:

Cartography and Map Styling
===========================

Copyright note
--------------

Since basemap data is taken from external sources, it is advised to cite and credit the data provider accordingly.

Some data providers request or at least recommend to include a corresponding note along with the map.
Depending on the policy, this may either be only the data source and licese or terms of use.

The corresponding information for the data providers included in cartopy is listed in the :ref:`data_copyright_table`:

.. todo::

    Turn typical licenses into http://docutils.sourceforge.net/docs/ref/rst/directives.html#replacement-text
    
.. todo::

    maybe even add the typical annotation text as an string property to each function / class

.. _data_copyright_table:

.. csv-table:: List of Data Providers and Form of Citation
   :file: ../_static/copyright_license.csv
   :delim: ;
   :header-rows: 1
   :widths: 10, 10, 10, 10, 10, 10

.. |copy| unicode:: 0xA9 .. copyright sign
.. |TM| unicode:: U+2122
   .. with trademark sign
.. |---| unicode:: U+02014 .. em dash
   :trim:

Here is an :ref:`example <data_copyright_plot>` showing such annotation for Natural Earth data:

.. _data_copyright_plot:

.. plot::
    :include-source:

    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnchoredText
    # adapt to your data source
    data_source = 'Natural Earth'
    # adapt to your data source
    data_license = 'public domain'

    # Create basic map for illustration.
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)

    # Add a text annotation in the bottom right corner.
    ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
    text = AnchoredText(r'$\mathcircled{c}$ ' + data_source + '; license: ' + data_license, loc=4,
                        prop={'size': 10}, frameon=True)
    ax.add_artist(text)

    plt.show()
