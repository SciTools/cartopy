.. _cartopy_developer_interfaces:

Cartopy developer interfaces
============================

Cartopy exposes several interfaces to help make it easier to add new
functionality quickly and easily.

..
    Feature API
    -----------
    TODO

    Image API
    ---------
    TODO


Data/Download API
-----------------

In order to keep the size of a cartopy release low, the majority of data is
not included as standard. This means that, when developing new features, it is
often necessary to provide interfaces which can acquire data from external
sources (normally via HTTP). The :class:`~cartopy.io.Downloader` class
has been designed to make this process as easy as possible for developers to
extend, whilst still making it possible for users to configure in their own
way.

.. autoclass:: cartopy.io.Downloader
    :members: FORMAT_KEYS, path, url, target_path,
              pre_downloaded_path, acquire_resource, from_config

An example of specialising this class can be found in
:mod:`cartopy.io.shapereader.NEShpDownloader` which enables the downloading of
zipped shapefiles from the `<http://NaturalEarthData.com>`_ website. All known
subclasses of :class:`~cartopy.io.Downloader` are listed below for
reference:

   * :class:`cartopy.io.shapereader.NEShpDownloader`
   * :class:`cartopy.io.srtm.SRTM3Downloader`