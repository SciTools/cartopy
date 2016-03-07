# (C) British Crown Copyright 2011 - 2016, Met Office
#
# This file is part of cartopy.
#
# cartopy is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cartopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cartopy.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division, print_function)

import contextlib
import os
import shutil
import tempfile
import warnings

from nose.tools import assert_equal, assert_raises

import cartopy
import cartopy.io as cio
from cartopy.io.shapereader import NEShpDownloader
from cartopy.tests.mpl.test_caching import CallCounter


def test_Downloader_data():
    di = cio.Downloader('https://testing.com/{category}/{name}.zip',
                        os.path.join('{data_dir}', '{category}',
                                     'shape.shp'),
                        '/project/foobar/{category}/sample.shp')

    replacement_dict = {'category': 'example',
                        'name': 'test',
                        'data_dir': os.path.join('/wibble', 'foo', 'bar')}

    assert_equal(di.url(replacement_dict),
                 'https://testing.com/example/test.zip')

    assert_equal(di.target_path(replacement_dict),
                 os.path.join('/wibble', 'foo', 'bar', 'example', 'shape.shp')
                 )

    assert_equal(di.pre_downloaded_path(replacement_dict),
                 '/project/foobar/example/sample.shp'
                 )


@contextlib.contextmanager
def config_replace(replacement_dict):
    """
    Provides a context manager to replace the ``cartopy.config['downloaders']``
    dict with the given dictionary. Great for testing purposes!

    """
    downloads_orig = cartopy.config['downloaders']
    cartopy.config['downloaders'] = replacement_dict
    yield
    cartopy.config['downloaders'] = downloads_orig


@contextlib.contextmanager
def download_to_temp():
    """
    Context manager which defaults the "data_dir" to a temporary directory
    which is automatically cleaned up on exit.

    """
    old_downloads_dict = cartopy.config['downloaders'].copy()
    old_dir = cartopy.config['data_dir']

    tmp_dir = tempfile.mkdtemp(suffix='_cartopy_data')
    cartopy.config['data_dir'] = tmp_dir
    try:
        yield tmp_dir
        cartopy.config['downloaders'] = old_downloads_dict
        cartopy.config['data_dir'] = old_dir
    finally:
        shutil.rmtree(tmp_dir)


def test_from_config():
    generic_url = 'https://example.com/generic_ne/{name}.zip'

    land_downloader = cio.Downloader(generic_url, '', '')
    generic_ne_downloader = cio.Downloader(generic_url, '', '')

    ocean_spec = ('shapefile', 'natural_earth', '110m', 'physical', 'ocean')
    land_spec = ('shapefile', 'natural_earth', '110m', 'physical', 'land')
    generic_spec = ('shapefile', 'natural_earth')

    target_config = {land_spec: land_downloader,
                     generic_spec: generic_ne_downloader,
                     }

    with config_replace(target_config):
        # ocean spec is not explicitly in the config, but a subset of it is,
        # so check that an appropriate downloader is returned
        r = cio.Downloader.from_config(ocean_spec)

        # check the resulting download item produces a sensible url.
        assert_equal(r.url({'name': 'ocean'}),
                     'https://example.com/generic_ne/ocean.zip')

        downloaders = cio.config['downloaders']

        r = cio.Downloader.from_config(land_spec)
        assert r is land_downloader


def test_downloading_simple_ascii():
    # downloads a file from the Google APIs. (very high uptime and file will
    # always be there - if this goes down, most of the internet would break!)
    # to test the downloading mechanisms.
    file_url = 'https://ajax.googleapis.com/ajax/libs/jquery/1.8.2/{name}.js'

    format_dict = {'name': 'jquery'}

    with download_to_temp() as tmp_dir:
        target_template = os.path.join(tmp_dir, '{name}.txt')
        tmp_fname = target_template.format(**format_dict)

        dnld_item = cio.Downloader(file_url, target_template)

        assert_equal(dnld_item.target_path(format_dict), tmp_fname)

        with warnings.catch_warnings(record=True) as w:
            assert_equal(dnld_item.path(format_dict), tmp_fname)

            assert len(w) == 1, ('Expected a single download warning to be '
                                 'raised. Got {}.'.format(len(w)))
            assert issubclass(w[0].category, cio.DownloadWarning)

        with open(tmp_fname, 'r') as fh:
            _ = fh.readline()
            assert_equal(" * jQuery JavaScript Library v1.8.2\n",
                         fh.readline())

        # check that calling path again doesn't try re-downloading
        with CallCounter(dnld_item, 'acquire_resource') as counter:
            assert_equal(dnld_item.path(format_dict), tmp_fname)
        assert counter.count == 0, 'Item was re-downloaded.'


def test_natural_earth_downloader():
    # downloads a file to a temporary location, and uses that temporary
    # location, then:
    #   * Checks that the file is only downloaded once even when calling
    #     multiple times
    #   * Checks that shapefiles have all the necessary files when downloaded
    #   * Checks that providing a path in a download item gets used rather
    #     than triggering another download

    tmp_dir = tempfile.mkdtemp()

    shp_path_template = os.path.join(tmp_dir,
                                     '{category}_{resolution}_{name}.shp')

    # picking a small-ish file to speed up download times, the file itself
    # isn't important - it is the download mechanism that is.
    format_dict = {'category': 'physical',
                   'name': 'rivers_lake_centerlines',
                   'resolution': '110m'}

    try:
        dnld_item = NEShpDownloader(target_path_template=shp_path_template)

        # check that the file gets downloaded the first time path is called
        with CallCounter(dnld_item, 'acquire_resource') as counter:
            shp_path = dnld_item.path(format_dict)
        assert counter.count == 1, 'Item not downloaded.'

        assert_equal(shp_path_template.format(**format_dict), shp_path)

        # check that calling path again doesn't try re-downloading
        with CallCounter(dnld_item, 'acquire_resource') as counter:
            assert_equal(dnld_item.path(format_dict), shp_path)
        assert counter.count == 0, 'Item was re-downloaded.'

        # check that we have the shp and the shx
        exts = ['.shp', '.shx']
        for ext in exts:
            stem = os.path.splitext(shp_path)[0]
            msg = "Shapefile's {0} file doesn't exist in {1}{0}".format(ext,
                                                                        stem)
            assert os.path.exists(stem + ext), msg

        # check that providing a pre downloaded path actually works
        pre_dnld = NEShpDownloader(target_path_template='/not/a/real/file.txt',
                                   pre_downloaded_path_template=shp_path
                                   )
        # check that the pre_dnld downloader doesn't re-download, but instead
        # uses the path of the previously downloaded item

        with CallCounter(pre_dnld, 'acquire_resource') as counter:
            assert_equal(pre_dnld.path(format_dict), shp_path)
        assert counter.count == 0, 'Aquire resource called more than once.'

    finally:
        shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
