Cartopy map gridlines and tick labels
=====================================

The :class:`~cartopy.mpl.gridliner.Gridliner` instance, often created by calling the
:meth:`cartopy.mpl.geoaxes.GeoAxes.gridlines` method on a
:class:`cartopy.mpl.geoaxes.GeoAxes` instance, has a variety of attributes which can be
used to determine draw time behaviour of the gridlines and labels.

.. important::

    The current :class:`~cartopy.mpl.gridliner.Gridliner` interface is likely to undergo
    a significant change in the versions following v0.6 in order to fix some of the underying
    limitations of the current implementation.

.. autoclass:: cartopy.mpl.gridliner.Gridliner
    :members:
    :undoc-members:


Sample Usage
------------
See :ref:`the gridliner features example <examples-gridliner>` for sample code
demonstrating some of the more advanced features.
