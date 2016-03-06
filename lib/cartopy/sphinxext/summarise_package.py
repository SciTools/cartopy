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

import inspect
import itertools
import os
import sys
import warnings
import six

import cartopy.tests


def walk_module(mod_name, exclude_folders=None):
    """
    Recursively walks the given module name.

    Returns:

        A generator of::

            (fully_qualified_import_name,
             root_directory_of_subpackage,
             fname_in_root_directory,
             sub_folders_in_root_directory)

    """
    __import__(mod_name)
    mod = sys.modules[mod_name]
    mod_dir = os.path.dirname(mod.__file__)
    exclude_folders = exclude_folders or []

    for root, folders, files in os.walk(mod_dir):
        for folder in exclude_folders:
            try:
                folders.remove(folder)
            except ValueError:
                pass

        # only allow python packages
        if '__init__.py' not in files:
            del folders[:]
            continue

        # Sort by filename and folder name.
        files.sort()
        folders.sort()

        def is_py_src(fname):
            root, ext = os.path.splitext(fname)
            return ext in ('.py', '.so')

        files = filter(is_py_src, files)

        for fname in files:
            sub_mod_name = mod_name
            relpath = os.path.relpath(root, mod_dir)
            if relpath == '.':
                relpath = ''

            for sub_mod in [f for f in relpath.split(os.path.sep) if f]:
                sub_mod_name += '.' + sub_mod

            if fname != '__init__.py':
                sub_mod_name += '.' + os.path.splitext(fname)[0]

            yield sub_mod_name, root, fname, folders


def objects_to_document(module_name):
    """
    Creates a generator of (obj_name, obj) that the given module of the
    given name should document.

    The module name may be any importable, including submodules
    (e.g. ``cartopy.io``)

    """
    try:
        __import__(module_name)
    except ImportError:
        warnings.warn('Failed to document {}'.format(module_name))
        return []
    module = sys.modules[module_name]
    elems = dir(module)

    if '__all__' in elems:
        document_these = [(obj, getattr(module, obj))
                          for obj in module.__all__]
    else:
        document_these = [(obj, getattr(module, obj)) for obj in elems
                          if not inspect.ismodule(getattr(module, obj)) and
                          not obj.startswith('_')]

        def is_from_this_module(x):
            return getattr(x[1], '__module__', '') == module_name

        document_these = filter(is_from_this_module, document_these)
        document_these = sorted(document_these,
                                key=lambda x: (str(type(x[1])),
                                               not x[0].isupper(),
                                               x[0]))

    # allow a developer to add other things to the documentation
    if hasattr(module, '__document_these__'):
        extra_objects_to_document = tuple((obj, getattr(module, obj))
                                          for obj in module.__document_these__)
        document_these = extra_objects_to_document + tuple(document_these)

    return document_these


def main(package_name, exclude_folders=None):
    """
    Return a string containing the rst that documents the given package name.

    """
    result = ''
    mod_walk = walk_module(package_name, exclude_folders=exclude_folders)
    for mod, _, fname, folders in mod_walk:
        for folder in folders:
            if folder.startswith('_'):
                folders.remove(folder)
        if fname.startswith('_') and fname != '__init__.py':
            continue

        result += '\n'
        result += mod + '\n'
        result += '*' * len(mod) + '\n'

        result += '\n'
        result += '.. currentmodule:: {}\n'.format(mod) + '\n'

        mod_objects = objects_to_document(mod)
        if mod_objects:
            result += '.. csv-table::\n' + '\n'

        table_elements = itertools.cycle(('\n\t', ) + (', ', ) * 3)
        for table_elem, (obj_name, _) in zip(table_elements, mod_objects):
            result += '{}:py:obj:`{}`'.format(table_elem, obj_name)

        result += '\n'

    return result


def gen_summary_rst(app):
    """
    Creates the rst file to summarise the desired packages.

    """
    package_names = app.config.summarise_package_names
    exclude_dirs = app.config.summarise_package_exclude_directories
    fnames = app.config.summarise_package_fnames

    if isinstance(package_names, six.string_types):
        package_names = [package_names]

    if package_names is None:
        raise ValueError('Please define a config value containing a list '
                         'of packages to summarise.')

    if exclude_dirs is None:
        exclude_dirs = [None] * len(package_names)
    else:
        exception = ValueError('Please provide a list of exclude directory '
                               'lists (one list for each package to '
                               'summarise).')

        if len(exclude_dirs) != len(package_names):
            raise exception

        for exclude_dirs_individual in exclude_dirs:
            if isinstance(exclude_dirs_individual, six.string_types):
                raise exception

    if fnames is None:
        fnames = ['outline_of_{}.rst'.format(package_name)
                  for package_name in package_names]
    else:
        if isinstance(fnames, six.string_types) or \
                len(fnames) != len(package_names):
            raise TypeError('Please provide a list of filenames for each of '
                            'the packages which are to be summarised.')

    outdir = app.builder.srcdir

    for package_name, out_fname, exclude_folders in zip(package_names,
                                                        fnames,
                                                        exclude_dirs):
        out_fpath = os.path.join(outdir, out_fname)
        content = main(package_name, exclude_folders=exclude_folders)
        with open(out_fpath, 'w') as fh:
            fh.write(content)


@cartopy.tests.not_a_nose_fixture
def setup(app):
    """
    Defined the Sphinx application interface for the summary generation.

    """
    app.connect('builder-inited', gen_summary_rst)

    # Allow users to define a config value to determine the names to summarise
    app.add_config_value('summarise_package_names', None, 'env')

    # Allow users to define a config value to determine the folders to exclude
    app.add_config_value('summarise_package_exclude_directories', None, 'env')

    # Allow users to define a config value to determine name of the output file
    app.add_config_value('summarise_package_fnames', None, 'env')
