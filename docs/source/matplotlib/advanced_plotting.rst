More advanced mapping with cartopy and matplotlib
=================================================

From the outset, cartopy's purpose has been to simplify and improve the quality of
mapping visualisations available for scientific data. Thanks to the simplicity of the cartopy
interface, in many cases the hardest part of producing such visualisations is getting
hold of the data in the first place. To address this, a Python package,
`Iris <http://scitools.org.uk/iris/>`_, has been created to make loading and saving data from a
variety of gridded datasets easier. Some of the following examples make use of the Iris
loading capabilities, while others use the netCDF4 Python package so as to show a range
of different approaches to data loading.


Contour plots
-------------

A very common way of displaying scalar data over a map is to draw contour
lines, or filled contours with a colormap.  This will display continuous
coloured regions, smoothing over any gridcell edges even if they are large.

See the :ref:`contour plot example <examples-contour-plot>`

.. figure:: ../gallery/scalar_data/images/sphx_glr_contour_plot_001.png
   :target: ../gallery/scalar_data/contour_plot.html
   :align: center
   :scale: 50


Block plots
-----------

To see the actual data grid, another approach is to display individual
gridcells on a plot, in colours mapped from the data values.  The individual
gridcells can then be seen on the plot, if they are large enough.

See the :ref:`blockplot example <examples-blockplot>`

.. figure:: ../gallery/scalar_data/images/sphx_glr_block_plot_001.png
   :target: ../gallery/scalar_data/block_plot.html
   :align: center
   :scale: 50


Images
------

An image can be overlaid on a map, if you provide the actual location of its
corners.

See the :ref:`image example <examples-image-display>`

.. figure:: ../gallery/scalar_data/images/sphx_glr_image_display_001.png
   :target: ../gallery/scalar_data/image_display.html
   :align: center
   :scale: 50


.. _vector_plotting:

Vector plotting
---------------

Cartopy comes with powerful vector field plotting functionality. There are 3 distinct options for
visualising vector fields:
:meth:`quivers <cartopy.mpl.geoaxes.GeoAxes.quiver>` (:ref:`example <examples-arrows>`),
:meth:`barbs <cartopy.mpl.geoaxes.GeoAxes.barbs>` (:ref:`example <examples-barbs>`) and
:meth:`streamplots <cartopy.mpl.geoaxes.GeoAxes.streamplot>` (:ref:`example <examples-streamplot>`)
each with their own benefits for displaying certain vector field forms.

.. figure:: ../gallery/vector_data/images/sphx_glr_arrows_001.png
   :target: ../gallery/vector_data/arrows.html
   :align: center
   :scale: 50

Since both :meth:`~cartopy.mpl.geoaxes.GeoAxes.quiver` and :meth:`~cartopy.mpl.geoaxes.GeoAxes.barbs`
are visualisations which draw every vector supplied, there is an additional option to "regrid" the
vector field into a regular grid on the target projection (done via
:func:`cartopy.vector_transform.vector_scalar_to_grid`). This is enabled with the ``regrid_shape``
keyword and can have a massive impact on the effectiveness of the visualisation:

.. figure:: ../gallery/vector_data/images/sphx_glr_regridding_arrows_001.png
   :target: ../gallery/vector_data/regridding_arrows.html
   :align: center
   :scale: 50
