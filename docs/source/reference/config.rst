.. _config:

Cartopy configuration
---------------------

.. module:: cartopy

The top level cartopy module contains the :attr:`~cartopy.config` dictionary which controls various aspects of cartopy's behaviour.

..
    n.b. cartopy.config docstring should be mirrored in lib/cartopy/__init__.py.


.. py:data:: config

    The config dictionary stores global configuration values for cartopy.

    In the first instance, the config is defined in ``cartopy/__init__.py``. It
    is possible to provide site wide customisations by including a
    ``siteconfig.py`` file along with the cartopy source code. ``siteconfig.py``
    should contain a function called ``update_config`` which takes the config
    dictionary instance as its first and only argument (from where it is
    possible to update the dictionary howsoever desired).

    For users without write permission to the cartopy source directory, a package
    called ``cartopy_userconfig`` should be made importable (consider putting it
    in ``site.getusersitepackages()``) and should expose a
    function called ``update_config`` which takes the config dictionary as its
    first and only argument.


    Keys in the config dictionary:

    ``pre_existing_data_dir``
        The absolute path to a directory where standard data (such as that from
        NaturalEarth) can be found. If it is not found in this location the
        ``data_dir`` config item will be used.

    ``data_dir``
        The absolute path to a directory where standard data (such as that from
        NaturalEarth) can be found. If it is not found and the item is
        downloadable cartopy will download the appropriate file(s) to a
        subdirectory of this directory, therefore ``data_dir`` should be
        writable by the user.

    ``repo_data_dir``
        The absolute path to the directory where the data delivered with the
        cartopy repository is stored. Typically this will only be set by OS
        packagers and system administrators for site wide deployments.

    ``downloaders``
        A dictionary mapping standard "specifications" to the appropriate
        :class:`~cartopy.io.Downloader`. For further documentation and an
        example see :func:`cartopy.io.Downloader.from_config`.
