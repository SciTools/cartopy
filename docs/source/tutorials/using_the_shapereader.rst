.. _using_the_shapereader:

Using the cartopy shapereader
=============================

Cartopy provides an object oriented shapefile reader based on top of the 
`pyshp <http://code.google.com/p/pyshp/>`_ module to provide easy, programmatic,
access to standard vector datasets.

.. currentmodule:: cartopy.io.shapereader

.. autoclass:: Reader
    :members:
    :undoc-members:

.. autoclass:: Record
    :members:
    :undoc-members:


Helper functions for shapefile acquisition
------------------------------------------

Cartopy provides an interface for access to frequently used data such as the
`GSHHS <http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html>`_ dataset and from
the `NaturalEarthData <http://www.naturalearthdata.com/>`_ website. 
These interfaces allow the user to define the data programmatically, and if the data does not exist
on disk, it will be retrieved from the appropriate source (normally by
downloading the data from the internet). Currently the interfaces available are:


.. autofunction:: natural_earth

.. autofunction:: gshhs


Using the shapereader
---------------------




We can acquire the countries dataset from Natural Earth found at
http://www.naturalearthdata.com/downloads/110m-cultural-vectors/110m-admin-0-countries/
by using the :func:`natural_earth` function:


.. testcode:: countries

    import cartopy.io.shapereader as shpreader

    shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')


By the nature of this shapfile from Natural Earth, it consists of a country
for each field.  Similarly, there are shapefiles available that consists of
boundaries, or prefectures etc.  Here, we can make use of the :class:`Reader`
to get the first country in the shapefile:

.. testcode:: countries

    reader = shpreader.Reader(shpfilename)
    countries = reader.records()
    country = countries.next()

The :func:`~Reader.records` method returns a generator,
which we can then succesively apply the next method to return successive
records.  The instance of the :class:`~Reader` class also contains a
:func:`~Reader.geometries` method which returns a generator for our
gemoetries.  Geometry is associated with the individual records,
however, if the metadata of the records is not neccessary to the user, then
they can be accessed directly from this generator.

.. note::
    :func:`~Reader.records` is usefull to extract/filter based on the properties
    of what the shape represents.

    :func:`~Reader.geometries` is usefull to extract/filter based on the
    properties of the shape.

We now get a country's attributes dictionary through
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

    # Define a function which returns the population given the country.
    population = lambda country: country.attributes['pop_est']

    # Sort the countries by population and get the first 5.
    countries_by_pop = sorted(reader.records(), key=population)[:5]

Which we can print with

.. doctest:: countries

    >>> ', '.join([country.attributes['name_long']
    ...            for country in countries_by_pop])
    'Western Sahara, French Southern and Antarctic Lands, Falkland Islands, Antarctica, Greenland'

Extracting the records for these entries then follows the same logic:

.. testcode:: countries

    >>> extracted_countries = [country for country in countries_by_pop]

See :doc:`plotting_shapefiles` for what we can do with these records/geometries.

**Excercises**:

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

 * **SHP.2**: Find the most populated country in each first letter group of
        "name_long" according to the population attribute "pop_est" in the
        shapefile.

        Hint: :func:`itertools.groupby` can help with the grouping and
        attributes["name_long"].

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
        sd_countries = sorted(reader.records(), key=first_letter)

        # group the countries by first letter
        for letter, countries in itertools.groupby(sd_countries, key=first_letter):
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

 * **SHP.3**: Extract the records of interest from the shapefile corresponding
        to the list of countries identified above.

    .. testcode:: countries
        :hide:

        # extract the records of interest
        extracted_countries = [
            sorted(countries, key=population)[-1] for
                letter, countries in itertools.groupby(sd_countries,
                                                       key=first_letter)]

        for record in extracted_countries:
            # print type and long name of each entry int the list
            print record.attributes['name_long'], type(record)

    .. testoutput:: countries

            Argentina <class 'cartopy.io.shapereader.Record'>
            Brazil <class 'cartopy.io.shapereader.Record'>
            China <class 'cartopy.io.shapereader.Record'>
            Democratic Republic of the Congo <class 'cartopy.io.shapereader.Record'>
            Ethiopia <class 'cartopy.io.shapereader.Record'>
            France <class 'cartopy.io.shapereader.Record'>
            Germany <class 'cartopy.io.shapereader.Record'>
            Hungary <class 'cartopy.io.shapereader.Record'>
            India <class 'cartopy.io.shapereader.Record'>
            Japan <class 'cartopy.io.shapereader.Record'>
            Kenya <class 'cartopy.io.shapereader.Record'>
            Lao PDR <class 'cartopy.io.shapereader.Record'>
            Mexico <class 'cartopy.io.shapereader.Record'>
            Nigeria <class 'cartopy.io.shapereader.Record'>
            Oman <class 'cartopy.io.shapereader.Record'>
            Pakistan <class 'cartopy.io.shapereader.Record'>
            Qatar <class 'cartopy.io.shapereader.Record'>
            Russian Federation <class 'cartopy.io.shapereader.Record'>
            South Africa <class 'cartopy.io.shapereader.Record'>
            Turkey <class 'cartopy.io.shapereader.Record'>
            United States <class 'cartopy.io.shapereader.Record'>
            Vietnam <class 'cartopy.io.shapereader.Record'>
            Western Sahara <class 'cartopy.io.shapereader.Record'>
            Yemen <class 'cartopy.io.shapereader.Record'>
            Zimbabwe <class 'cartopy.io.shapereader.Record'>

Each record has an associated geometry object.
For examples of plotting shapefiles, see :doc:`plotting_shapefiles`
