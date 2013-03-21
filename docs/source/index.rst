.. toctree::
    :hidden:

    whats_new.rst
    building_from_source/index.rst
    matplotlib/intro.rst
    tutorials/using_the_shapereader.rst
    matplotlib/feature_interface.rst
    projections/table.rst
    examples/*.rst


Introduction
============

Cartopy is a Python package designed to make drawing maps for data analysis easy and enjoyable.

Primarily the non-drawing functionality of cartopy makes use of proj.4, numpy, geos and shapely; 
with simple and intuitive interfacing to matplotlib for map drawing.

Some of the key features of cartopy are:

 * object oriented projection definitions
 * point, line, polygon and image transformations between projections
 * integration to expose advanced mapping in matplotlib with a simple and intuitive interface
 * work-in-progress mechanisms for accessing specialist data such as those from the "Shuttle Radar Topography 
   Mission" (SRTM) and the "Global Self-consistent, Hierarchical, High-resolution Shoreline" database (GSHHS).

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

To start using cartopy with matplotlib see the :doc:`matplotlib/intro` section.

View the list of cartopy projections on the :doc:`projections/table`.

For shapefile interaction with cartopy see :doc:`tutorials/using_the_shapereader`
and :doc:`matplotlib/feature_interface`.

Also, check out our :doc:`what's new <whats_new>` page.


Getting involved
================

Cartopy was originally developed at the UK Met Office to allow scientists to visualise 
their data on maps quickly, easily and most importantly, accurately. 
Cartopy has since been made freely available under the terms of the 
`GNU Lesser General Public License <http://www.gnu.org/licenses/lgpl.html>`_ and is suitable to be used in a variety
of scientific fields.

There are many ways to get involved in the development of cartopy:

 * Report bugs to https://github.com/SciTools/cartopy/issues (please look to see if there are any outstanding
   bugs which cover the issue before making a new one).
 * Helping others with support questions on `stackoverflow <http://stackoverflow.com/questions/tagged/cartopy>`_. 
 * Contribute to the documentation fixing typos, adding examples, explaining things more clearly, or even
   re-designing its layout/logos etc..
 * Contribute bug fixes (:issues:`a list of outstanding bugs can be found on github <bug>`).
 * Contributing enhancements and new features 
   (:issues:`a wish list of features can also be found on github <wishlist>`).
                          
If you or your organisation use cartopy, we would love to hear about it. 

Development discussion takes place on the `matplotlib-devel mailing list
<http://matplotlib.1069221.n5.nabble.com/matplotlib-devel-f28077.html>`_ with all discussion
subjects prefixed with "cartopy: ".


Citing  Cartopy
===============
To cite Cartopy see the :doc:`citation` section.


