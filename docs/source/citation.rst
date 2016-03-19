.. _Citing_Cartopy:

Citing cartopy
==============

If cartopy played an important part in your research then please add us to your reference list by using one of the recommendations below.

************
BibTeX entry 
************

For example::

 @manual{Cartopy,
 author = {{Met Office}},
 title = {Cartopy: a cartographic python library with a matplotlib interface},
 year = {2010 - 2015},
 address = {Exeter, Devon },
 url = {http://scitools.org.uk/cartopy}
 } 


*******************
Downloaded Software
*******************

ProductName. Version. ReleaseDate. Publisher. Location. DOIorURL. DownloadDate.

For example::

 Cartopy. v0.11.2. 22-Aug-2014. Met Office. UK. https://github.com/SciTools/cartopy/archive/v0.11.2.tar.gz


********************
Checked out Software
********************

ProductName. Publisher. URL. CheckoutDate. RepositorySpecificCheckoutInformation.

For example::

 Cartopy. Met Office. git@github.com:SciTools/cartopy.git. 2015-02-18. 7b2242e.

.. _How to cite and describe software: http://software.ac.uk/so-exactly-what-software-did-you-use


.. [Jackson] Jackson, M. 2012. `How to cite and describe software`_. Accessed 2013-03-06.


.. _referencing:

**********************
Referencing (the data)
**********************


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

An :ref:`example <data_copyright_plot>` showing such annotation for Natural Earth data is shown in the `feature_creation example </examples/feature_creation.html#feature-creation-example>`_