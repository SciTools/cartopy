.. toctree::
    :hidden:

    whats_new.rst
    installing.rst
    cartopy_outline.rst


Introduction
============

Cartopy is a Python package designed for geospatial data processing in order to produce maps and other geospatial data analyses.

Cartopy makes use of the powerful PROJ.4, NumPy and Shapely libraries and includes a programmatic interface
built on top of Matplotlib for the creation of publication quality maps.

Key features of cartopy are its object oriented :ref:`projection definitions <cartopy_projections>`, and its
ability to transform points, lines, vectors, polygons and images between those projections.

You will find cartopy especially useful for large area / small scale data, where Cartesian
assumptions of spherical data traditionally break down. If you've ever experienced a singularity
at the pole or a cut-off at the dateline, it is likely you will appreciate cartopy's unique features!


.. _getting-started-with-cartopy:

Getting started
===============

The :ref:`installation guide <installing>` provides information on getting up and running. 
Cartopy's documentation is arranged in userguide form, with reference documentation available inline.

.. toctree::
    :maxdepth: 1

    crs/index.rst
    crs/projections.rst
    matplotlib/intro.rst
    matplotlib/feature_interface.rst
    tutorials/understanding_transform.rst
    tutorials/using_the_shapereader.rst
    developer_interfaces.rst


The :doc:`outline <cartopy_outline>` link found above the cartopy logo on all pages can be used to
quickly find the reference documentation for known classes or functions.

For those updating from an older version of cartopy, the :doc:`what's new <whats_new>` page
outlines recent changes, new features, and future development plans.


..
  All the stuff that isn't in chapter form.

.. toctree::
    :hidden:

    cartopy.rst
    cartopy/geodesic.rst
    cartopy/io/img_tiles.rst
    cartopy/io/ogc_clients.rst
    cartopy/trace.rst
    cartopy/util/util.rst


Getting involved
================

Cartopy was originally developed at the UK Met Office to allow scientists to visualise
their data on maps quickly, easily and most importantly, accurately.
Cartopy has been made freely available under the terms of the
`GNU Lesser General Public License <https://www.gnu.org/licenses/lgpl.html>`_.
It is suitable to be used in a variety
of scientific fields and has an :doc:`active development community <contributors>`.

There are many ways to get involved in the development of cartopy:

 * If you write a paper which makes use of cartopy, please consider :doc:`citing <citation>` it.
 * If you publish the maps and plots, please consider the required :ref:`attribution of copyright <referencing_copyright>` for the data.
 * Report bugs and problems with the code or documentation to https://github.com/SciTools/cartopy/issues (please look to see if there are any outstanding
   bugs which cover the issue before making a new one).
 * Help others with cartopy questions on `StackOverflow <https://stackoverflow.com/questions/tagged/cartopy>`_.
 * Contribute to the documentation fixing typos, adding examples, explaining things more clearly, or even
   re-designing its layout/logos etc.. The `documentation source <https://github.com/SciTools/cartopy>`_ is kept in the same repository as the source code.
 * Contribute bug fixes (:issues:`a list of outstanding bugs can be found on github <bug>`).
 * Contribute enhancements and new features on the issue tracker.
 * Chat with users and developers in the `Gitter chat room <https://gitter.im/SciTools/cartopy>`_.



.. toctree::
    :glob:
    :hidden:

    gallery/index
    citation
    copyright
    contributors
