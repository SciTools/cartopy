.. _using_the_shapereader:

Using the cartopy shapereader
=============================

Cartopy provides an object oriented shapefile reader based on top of the 
`pyshp`_ module to provide easy, programmatic, access to standard vector datasets.

Cartopy's wrapping of pyshp has the benefit of being pure python, and is therefore
easy to install and extremely portable. However, for heavy duty shapefile I/O `Fiona`_ and
`GeoPandas`_ are highly recommended.

.. _pyshp: https://github.com/GeospatialPython/pyshp
.. _Fiona: https://fiona.readthedocs.io/
.. _GeoPandas: http://geopandas.org/

Detailed API for the shapereader functionality can be found in the :ref:`documentation <api.io.shapereader>`

.. currentmodule:: cartopy.io.shapereader

Helper functions for shapefile acquisition
-------------------------------------------

Cartopy provides an interface for access to frequently used data such as the
`GSHHS <https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html>`_ dataset and from
the `NaturalEarthData <https://www.naturalearthdata.com/>`_ website.
These interfaces allow the user to define the data programmatically, and if the data does not exist
on disk, it will be retrieved from the appropriate source (normally by
downloading the data from the internet). Currently the interfaces available are
:func:`natural_earth` and :func:`gshhs`.


Using the shapereader
---------------------

We can acquire the countries dataset from Natural Earth found at
https://www.naturalearthdata.com/downloads/110m-cultural-vectors/110m-admin-0-countries/
by using the :func:`natural_earth` function:


.. testcode:: countries

    import cartopy.io.shapereader as shpreader

    shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')


From here, we can make use of the :class:`Reader` to get the first country
in the shapefile:

.. testcode:: countries

    reader = shpreader.Reader(shpfilename)
    countries = reader.records()
    country = next(countries)

We can get the country's attributes dictionary with the
:data:`Record.attributes` attribute:

.. doctest:: countries
    :options: +ELLIPSIS

    >>> print type(country.attributes)
    <type 'dict'>
    >>> print sorted(country.attributes.keys())
    ['abbrev', ..., 'name_long', ... 'pop_est', ...]

We could now find the 5 least populated countries with:

.. testcode:: countries

    reader = shpreader.Reader(shpfilename)

    # define a function which returns the population given the country
    population = lambda country: country.attributes['pop_est']

    # sort the countries by population and get the first 5
    countries_by_pop = sorted(reader.records(), key=population)[:5]

Which we can print with

.. doctest:: countries

    >>> ', '.join([country.attributes['name_long']
    ...            for country in countries_by_pop])
    'Western Sahara, French Southern and Antarctic Lands, Falkland Islands, Antarctica, Greenland'


**Exercises**:

 * **SHP.1**: Repeat the last example to show the 4 most populated African countries in to the shapefile.

   Hint: Look at the possible attributes to find out which continent a country belongs.

   Answer:

    .. testcode:: countries
        :hide:

        reader = shpreader.Reader(shpfilename)

        # define a function which can return the population of a given country
        population = lambda country: country.attributes['pop_est']

        # sort the countries by population
        countries_by_pop = sorted(reader.records(), key=population)

        # define a function which can return whether a country belongs to Africa
        is_african = lambda country: country.attributes['continent'] == 'Africa'

        # remove non-African countries
        african_countries = filter(is_african, countries_by_pop)

        print ', '.join([country.attributes['name_long']
                         for country in african_countries[-4:]])

    .. testoutput:: countries

        Democratic Republic of the Congo, Egypt, Ethiopia, Nigeria

 * **SHP.2**: Using the countries shapefile, find the most populated country grouped 
   by the first letter of the "name_long".

   Hint: :func:`itertools.groupby` can help with the grouping.

   Answer:

   .. testcode:: countries
        :hide:

        import itertools

        reader = shpreader.Reader(shpfilename)

        # define a function which returns the first letter of a country's name
        first_letter = lambda country: country.attributes['name_long'][0]
        # define a function which returns the population of a country
        population = lambda country: country.attributes['pop_est']

        # sort the countries so that the groups come out alphabetically
        countries = sorted(reader.records(), key=first_letter)

        # group the countries by first letter
        for letter, countries in itertools.groupby(countries, key=first_letter):
            # print the letter and least populated country
            print letter, sorted(countries, key=population)[-1].attributes['name_long']

   .. testoutput:: countries

            A Argentina
            B Brazil
            C China
            D Democratic Republic of the Congo
            E Ethiopia
            F France
            G Germany
            H Hungary
            I India
            J Japan
            K Kenya
            L Lao PDR
            M Mexico
            N Nigeria
            O Oman
            P Pakistan
            Q Qatar
            R Russian Federation
            S South Africa
            T Turkey
            U United States
            V Vietnam
            W Western Sahara
            Y Yemen
            Z Zimbabwe

