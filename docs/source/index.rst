.. toctree::
    :hidden:

    whats_new.rst
    building_from_source/index.rst
    citation.rst
    cartopy_outline.rst


Introduction
============

Cartopy is a Python package designed to make drawing maps for data analysis and visualisation as easy as possible.

Cartopy makes use of the powerful PROJ.4, numpy and shapely libraries and has a simple and intuitive
drawing interface to matplotlib for creating publication quality maps.

Some of the key features of cartopy are:

 * object oriented projection definitions
 * point, line, polygon and image transformations between projections
 * integration to expose advanced mapping in matplotlib with a simple and intuitive interface
 * powerful vector data handling by integrating shapefile reading with Shapely capabilities 

Cartopy is licensed under `GNU Lesser General Public License <http://www.gnu.org/licenses/lgpl.html>`_ 
and is at version |version|. You can find the source code for cartopy on 
`our github page <http://github.com/SciTools/cartopy>`_.


Installation
============

Installation of cartopy can currently only be done from source. 
Build instructions for specific operating systems can be found in the 
:ref:`building from source <building_from_source>` section.


Getting started
===============

The cartopy documentation is arranged in a userguide form with reference documentation available inline.

.. toctree::
    :maxdepth: 1

    crs/index.rst
    crs/projections.rst
    tutorials/using_the_shapereader.rst
    matplotlib/intro.rst
    matplotlib/feature_interface.rst

    developer_interfaces.rst


The :doc:`outline <cartopy_outline>` link found above the cartopy logo on all pages can be used to
quickly find the reference documentation for known classes or functions.

For those updating from an older version of cartopy, the :doc:`what's new <whats_new>` page
outlines recent changes, new features, and future development plans.


Getting involved
================

Cartopy was originally developed at the UK Met Office to allow scientists to visualise 
their data on maps quickly, easily and most importantly, accurately. 
Cartopy has been made freely available under the terms of the
`GNU Lesser General Public License <http://www.gnu.org/licenses/lgpl.html>`_. It is suitable to be used in a variety
of scientific fields and has an :doc:`active development community <contributors>`.

There are many ways to get involved in the development of cartopy:

 * If you write a paper which makes use of cartopy, please consider :doc:`citing <citation>` it.
 * Report bugs to https://github.com/SciTools/cartopy/issues (please look to see if there are any outstanding
   bugs which cover the issue before making a new one).
 * Help others with support questions on `stackoverflow <http://stackoverflow.com/questions/tagged/cartopy>`_. 
 * Contribute to the documentation fixing typos, adding examples, explaining things more clearly, or even
   re-designing its layout/logos etc..
 * Contribute bug fixes (:issues:`a list of outstanding bugs can be found on github <bug>`).
 * Contribute enhancements and new features 
   (:issues:`a wish list of features can also be found on github <wishlist>`).

Development discussion takes place on the `matplotlib-devel mailing list
<http://matplotlib.1069221.n5.nabble.com/matplotlib-devel-f28077.html>`_ with all discussion
subjects prefixed with "cartopy: ".



.. toctree::
    :glob:
    :hidden:

    gallery.rst
    examples/*.rst
    copyright.rst
