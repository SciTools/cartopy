# (C) British Crown Copyright 2011 - 2012, Met Office
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
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.


import contextlib
import os
import shutil
import tempfile
import warnings

from nose.tools import assert_equal

import cartopy
import cartopy.io as cio
from cartopy.io.shapereader import NEShpDownloader


def test_DownloadableItem_data():
    di = cio.DownloadableItem('http://testing.com/{category}/{name}.zip',
                              os.path.join('{data_dir}', '{category}', 
                                           'shape.shp'),
                              '/project/foobar/{category}/sample.shp',
                              )

    replacement_dict = {'category': 'example', 
                        'name': 'test',
                        'data_dir': os.path.join('/wibble', 'foo', 'bar'),
                        } 

    assert_equal(di.url(replacement_dict), 
                 'http://testing.com/example/test.zip')
    
    # not having all the appropriate keys should fail gracefully
    assert_equal(di.url({}), 
                 'http://testing.com/{category}/{name}.zip')
    
    assert_equal(di.target_path(replacement_dict),
                 os.path.join('/wibble', 'foo', 'bar', 'example', 'shape.shp')
                 )
                 
    assert_equal(di.pre_downloaded_path(replacement_dict),
                 '/project/foobar/example/sample.shp'
                 )


@contextlib.contextmanager
def config_replace(replacement_dict):
    """
    Provides a context manager to replace the ``cartopy.config['downloads']``
    dict with the given dictionary. Great for testing purposes!
    
    """
    downloads_orig = cartopy.config['downloads']
    cartopy.config['downloads'] = replacement_dict
    yield
    cartopy.config['downloads'] = downloads_orig


@contextlib.contextmanager
def download_to_temp():
    """
    Context manager which defaults the "data_dir" to a temporary directory
    which is automatically cleaned up on exit.
     
    """
    old_downloads_dict = cartopy.config['downloads'].copy()
    old_dir = cartopy.config['data_dir']
    
    tmp_dir = tempfile.mkdtemp(suffix='cartopy_data')
    cartopy.config['data_dir'] = tmp_dir
    try:
        yield tmp_dir
        cartopy.config['downloads'] = old_downloads_dict
        cartopy.config['data_dir'] = old_dir
    finally:
        shutil.rmtree(tmp_dir)



def test_from_config():
    generic_url = 'http://example.com/generic_ne/{name}.zip'
    
    land_downloader = cio.DownloadableItem(generic_url, '', '')
    generic_ne_downloader = cio.DownloadableItem(generic_url, '', '')
    
    ocean_spec = ('shapefile', 'natural_earth', '110m', 'physical', 'ocean')
    land_spec = ('shapefile', 'natural_earth', '110m', 'physical', 'land')
    generic_spec = ('shapefile', 'natural_earth')
    
    target_config = {land_spec: land_downloader,
                     generic_spec: generic_ne_downloader,
                     }
     
    with config_replace(target_config):
        # ocean spec is not explicitly in the config, but a subset of it is,
        # so check that an appropriate downloader is returned
        r = cio.DownloadableItem.from_config(ocean_spec)
        
        # check the resulting download item produces a sensible url.
        assert_equal(r.url({'name': 'ocean'}), 
                     'http://example.com/generic_ne/ocean.zip')
        
        downloads = cio.config['downloads']
        
        # check there has been a copy operation
        assert downloads[generic_spec] is not downloads[ocean_spec]

        r = cio.DownloadableItem.from_config(land_spec)
        assert r is land_downloader


def test_downloading_simple_ascii():
    # downloads a file from the Google APIs. (very high uptime and file will 
    # always be there - if this goes down, most of the internet would break!)
    # to test the downloading mechanisms.
    file_url = 'http://ajax.googleapis.com/ajax/libs/jquery/1.8.2/{name}.js'
    
    format_dict = {'name': 'jquery'}
    
    suffix = format_dict['name'] + '.txt'
    
    # make a temp file, get its name, then immediately delete it. we use the
    # filename for saving a DownloadableItem
    with tempfile.NamedTemporaryFile('r', suffix='_' + suffix) as tmp_fh:
        tmp_fname = tmp_fh.name
        
    try:
        # turn the explicit filename from tempfile into a template,
        # removing the "test1.txt" from the filename.
        target_template = tmp_fname[:-len(suffix)] + '{name}.txt'
        
        dnld_item = cio.DownloadableItem(file_url, target_template)
        
        assert_equal(dnld_item.target_path(format_dict), tmp_fname)
        
        with warnings.catch_warnings(record=True) as w:
            result_path = dnld_item.path(format_dict)
            assert len(w) == 1
            assert issubclass(w[0].category, cio._DownloadWarning)
            assert_equal(str(w[0].message), 
                         'Downloading: ' + file_url.format(name='jquery')) 
            
        assert_equal(result_path, tmp_fname)
        
        with open(tmp_fname, 'r') as fh:
            _ = fh.readline()
            assert_equal(" * jQuery JavaScript Library v1.8.2\n",
                             fh.readline())
        
        # check that calling path again doesn't try re-downloading    
        def no_way(*args, **kwargs):
            raise ValueError('Repeated path calls result in repeated '
                             'downloading of the file.')
        dnld_item.__class__.acquire_resource = no_way
        
        assert_equal(dnld_item.path(format_dict), result_path)
    
    # remove the tmp_fname no matter what happens
    finally:
        os.remove(tmp_fname)
    
    
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
                   'resolution': '110m'
                  } 
    
    try:
        dnld_item = NEShpDownloader(target_path_template=shp_path_template)
        
        with warnings.catch_warnings(record=True) as w:
            shp_path = dnld_item.path(format_dict)
            assert len(w) == 1
            assert issubclass(w[0].category, cio._DownloadWarning)
            
        assert_equal(shp_path_template.format(**format_dict), shp_path)
        
        # check that calling path again doesn't try re-downloading    
        def no_way(*args, **kwargs):
            raise ValueError('Subsequent calls to path resulted in '
                             'downloading the resource multiple times')
        dnld_item.__class__.acquire_resource = no_way
        
        assert_equal(dnld_item.path(format_dict), shp_path)
        
        # check that we have the shp and the shx
        exts = ['.shp', '.shx']
        for ext in exts:
            stem = os.path.splitext(shp_path)[0]
            assert os.path.exists(stem + ext), ("Shapefile's {} file doesn't "
                                                "exist in {1}{0}".format(ext,
                                                                      stem))
        
        # check that providing a pre downloaded path actually works 
        pre_dnld = NEShpDownloader(target_path_template=shp_path_template,
                                   pre_downloaded_path_template=shp_path
                                   )
        # check that the pre_dnld downloader doesn't re-download, but instead
        # uses the path of the previously downloaded item
        pre_dnld.__class__.acquire_resource = no_way
        assert_equal(pre_dnld.path(format_dict), shp_path)
            
    finally:
        shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
