# (C) British Crown Copyright 2013 - 2016, Met Office
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

import os.path
import sys

from cartopy.sphinxext.summarise_package import walk_module
import cartopy.tests


def out_of_date(original_fname, target_fname):
    """
    Checks to see if the ``target_fname`` exists, and if so, whether
    the modification timestamp suggests that ``original_fname`` is newer
    than ``target_fname``.

    """
    return (not os.path.exists(target_fname) or
            os.path.getmtime(original_fname) > os.path.getmtime(target_fname))


def same_contents(fname, contents_str):
    """
    Checks to see if the given fname contains the contents given by
    ``contents_str``. The result could be used to determine if the
    contents need to be written to the given fname.

    """
    if os.path.exists(fname):
        with open(fname, 'r') as fh:
            return fh.read() == contents_str
    return False


def parent_module(module):
    """
    Returns the direct module ascendent of the given module.

    For example, giving a module ``a.b.c`` would return the
    ``a.b`` module.

    If the module is a top level package, None will be returned.

    .. note::

        Requires the __name__ attribute on the given module
        to be correctly defined and the parent module to be importable.

    >>> import numpy.ma.core
    >>> from cartopy.sphinxext.gallery import parent_module
    >>> parent_module(numpy.ma.core) # doctest: +ELLIPSIS
    <module 'numpy.ma' from '...'>

    """
    result = None

    name = module.__name__
    by_dot = name.split('.')
    if len(by_dot) > 1:
        parent_name = '.'.join(by_dot[:-1])
        result = get_module(parent_name)

    return result


def get_module(mod_name):
    """
    Return the module instance of the given name, by importing
    it into the system.

    """
    __import__(mod_name)
    return sys.modules[mod_name]


def safe_mod_name_and_fname(module_name, ancestor_module_name):
    """
    Returns a safe module name (for linking too) and safe filename (suffixed
    with ".rst") for the given module name, relative to the given ancestor
    module.


    >>> from cartopy.sphinxext.gallery import safe_mod_name_and_fname
    >>> safe_mod_name_and_fname('numpy.ma.core', 'numpy')
    ('ma-core', 'ma/core.rst')

    """
    mod = get_module(module_name)
    ancestor_package = get_module(ancestor_module_name)

    safe_fname = os.path.relpath(os.path.splitext(mod.__file__)[0],
                                 os.path.dirname(ancestor_package.__file__))
    safe_name = safe_fname.replace(os.path.sep, '-')
    safe_fname = safe_fname + '.rst'
    return safe_name, safe_fname


def examples_code(examples_mod_name,
                  source_directory,
                  output_directory='examples'):
    """
    Generates the rst code for the given examples module.

    examples_mod_name - the name of the parent (importable) module which
                        should be recursively walked when looking for
                        examples

    source_directory - the path to the top level source directory containing
                       the rst content of the documentation

    output_directory - the directory for the output to be put in. Should be
                       relative to the source_directory

    """
    for mod_name, root_dir, fname, _ in walk_module(examples_mod_name):
        if fname.startswith('__init__'):
            continue

        rst, rst_fname = individual_example_rst(mod_name, examples_mod_name,
                                                output_directory)
        rst_fname = os.path.join(source_directory, output_directory, rst_fname)

        py_fname = os.path.join(source_directory, output_directory,
                                os.path.splitext(rst_fname)[0] + '.py')

        if not os.path.isdir(os.path.dirname(py_fname)):
            os.makedirs(os.path.dirname(py_fname))

        if out_of_date(os.path.join(root_dir, fname), py_fname):
            with open(os.path.join(root_dir, fname), 'r') as in_fh:
                with open(py_fname, 'w') as out_fh:
                    for line in in_fh:
                        # Crudely remove the __tags__ line.
                        if line.startswith('__tags__ = '):
                            continue
                        out_fh.write(line)

        if not same_contents(rst_fname, rst):
            with open(rst_fname, 'w') as fh:
                fh.write(rst)


def gallery_code(examples_mod_name):
    """
    Returns rst code suitable for generating a html gallery using sphinx.

    examples_mod_name - the name of the importable (sub)module which contains
                        the examples

    """
    # Store a dictionary mapping tag_name to (mod_name, mod_instance, tags)
    examples_by_tag = {}

    for mod_name, _, _, _ in walk_module(examples_mod_name):
        if mod_name != examples_mod_name:
            __import__(mod_name)
            mod = sys.modules[mod_name]
            tags = getattr(mod, '__tags__', ['Miscellanea'])

            for tag in tags:
                examples_by_tag.setdefault(tag, []).append((mod_name, mod))

    result = ['Gallery',
              '=======',
              'Tags:\n',
              '.. container:: inline-paragraphs\n'
              ]

    examples_by_tag = sorted(iter(examples_by_tag.items()),
                             key=lambda pair: (pair[0] == 'Miscellanea',
                                               pair[0]))

    for tag, _ in examples_by_tag:
        result.append('\t:ref:`gallery-tag-{}`\n'.format(tag))

    for tag, mod_list in examples_by_tag:
        result.extend(['.. _gallery-tag-{}:\n'.format(tag),
                       '{}'.format(tag),
                       '-' * len(tag) + '\n',
                       '.. container:: gallery_images\n'])

        for (mod_name, mod) in mod_list:
            safe_name, _ = safe_mod_name_and_fname(mod_name,
                                                   examples_mod_name)

            # XXX The path is currently determined out of process by
            # the plot directive. It would be nice to figure out the
            # naming scheme to handle multiple plots in a single example.
            img_path = 'examples/{}_00_00.png'.format(
                mod_name.split('.')[-1])
            thumb_path = 'examples/{}_00_00.thumb.png'.format(
                mod_name.split('.')[-1])

            entry = ["|image_{}|_\n".format(safe_name),
                     ".. |image_{}| image:: {}".format(safe_name, thumb_path),
                     # XXX Hard-codes the html - rst cannot do nested inline
                     # elements (very well).
                     ".. _image_{}: examples/{}.html".format(
                         safe_name, safe_name)]

            result.extend(['\n\n\t' + line for line in entry])

    return '\n'.join(result)


def individual_example_rst(example_mod_name, examples_mod_name,
                           output_directory):
    """
    Generates the rst code for the given example and suggests a sensible
    rst filename.

    example_mod_name - the name of the importable submodule which contains
                       the individual example which is to be documented

    examples_mod_name - the name of the importable (sub)module which contains
                        the examples

    output_directory - is a path to the examples output directory, relative
                       to the source directory.

    """
    mod = get_module(example_mod_name)
    safe_name, safe_fname = safe_mod_name_and_fname(example_mod_name,
                                                    examples_mod_name)
    fname_base = os.path.splitext(safe_fname)[0] + '.py'
    example_code_fname = os.path.join(output_directory, fname_base)

    result = ['.. _examples-{}:\n'.format(safe_name)]

    if mod.__doc__:
        result.append(mod.__doc__ + '\n')
    else:
        result.extend(['{} example'.format(safe_name),
                       '-' * (len(safe_name) + 8) + '\n'])

    result.extend(['.. plot:: {}\n'.format(example_code_fname),
                   '.. literalinclude:: {}\n'.format(fname_base)])

    return '\n'.join(result), safe_fname


def gen_gallery(app):
    """Produces the gallery rst file."""
    example_package_name = app.config.examples_package_name
    fname = app.config.gallery_name

    if example_package_name is None:
        raise ValueError('A config value for gallery_package_name should be '
                         'defined.')

    outdir = app.builder.srcdir
    fname = os.path.join(outdir, fname)

    gallery_rst = gallery_code(example_package_name)

    if not same_contents(fname, gallery_rst):
        with open(fname, 'w') as fh:
            fh.write(gallery_rst)


def gen_examples(app):
    """Produces the examples directory."""
    example_package_name = app.config.examples_package_name

    source_dir = app.builder.srcdir

    examples_code(example_package_name, source_dir, 'examples')


@cartopy.tests.not_a_nose_fixture
def setup(app):
    app.connect('builder-inited', gen_gallery)
    app.connect('builder-inited', gen_examples)

    # Allow users to define a config value to determine the name of the
    # gallery rst file (with file extension included)
    app.add_config_value('gallery_name', 'gallery.rst', 'env')

    # Allow users to define a config value to determine the name of the
    # importable examples (sub)package
    app.add_config_value('examples_package_name', None, 'env')
