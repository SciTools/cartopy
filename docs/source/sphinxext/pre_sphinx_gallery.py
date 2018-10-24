# -*- coding: utf-8 -*-
"""
Override sphinx_gallery's treatment of groups (folders) with cartopy's
``__tags__`` semantics. This is tightly bound to the sphinx_gallery
implementation, hence the explicit version checking.

This code was modified from sphinx-gallery:

Copyright (c) 2015, Óscar Nájera
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of sphinx-gallery nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
from collections import OrderedDict
import os.path
import shutil
import tempfile
import textwrap

import sphinx_gallery
if sphinx_gallery.__version__ not in ['0.1.12']:  # noqa: E402
    raise RuntimeError('not tested with this version of sphinx_gallery ({}). '
                       'Please modify this check, and validate sphinx_gallery'
                       ' behaves as expected.'
                       ''.format(sphinx_gallery.__version__))

import sphinx_gallery.gen_gallery
import sphinx_gallery.gen_rst
from sphinx_gallery.gen_rst import (
    write_backreferences, extract_intro, _thumbnail_div,
    generate_file_rst, sphinx_compatibility)


GALLERY_HEADER = textwrap.dedent("""

    Cartopy Gallery
    ---------------

    The following visual examples demonstrate some of the functionality of
    Cartopy, particularly its matplotlib interface.

    For a structured introduction to cartopy, including some of these
    examples, see :ref:`getting-started-with-cartopy`.

""")


def example_groups(src_dir):
    """Return a dictionary of {tag: [example filenames]} for the given dir."""

    sorted_listdir = [fname for fname in sorted(os.listdir(src_dir))
                      if fname.endswith('.py') and not fname.startswith('_')]
    tagged_examples = {}

    for fname in sorted_listdir:
        fpath = os.path.join(src_dir, fname)
        with open(fpath, 'r') as fh:
            for line in fh:
                # Crudely remove the __tags__ line.
                if line.startswith('__tags__ = '):
                    exec(line.strip(), locals(), globals())
                    for tag in __tags__:  # noqa:
                        tagged_examples.setdefault(tag, []).append(fname)
                    break
            else:
                tag = 'Miscellanea'
                tagged_examples.setdefault(tag, []).append(fname)
    return tagged_examples


def order_examples(tagged_examples):
    """Order the tags and their examples."""
    preferred_tag_order = ['Introductory',
                           'Lines and polygons',
                           'Scalar data',
                           'Vector data',
                           'Web services']

    def sort_key(item):
        tag = item[0]
        try:
            index = preferred_tag_order.index(tag)
        except ValueError:
            index = len(preferred_tag_order) + 1

        return (index, tag.lower())
    sorted_items = sorted(tagged_examples.items(), key=sort_key)
    return OrderedDict(sorted_items)


def write_example(src_fpath, target_dir):
    target_fpath = os.path.join(target_dir, os.path.basename(src_fpath))
    with open(src_fpath, 'r') as fh:
        with open(target_fpath, 'w') as fh_out:
            for line in fh:
                # Crudely remove the __tags__ line.
                if line.startswith('__tags__ = '):
                    continue
                fh_out.write(line)


def generate_dir_rst(src_dir, target_dir, gallery_conf, seen_backrefs):
    """Generate the gallery reStructuredText for an example directory"""

    fhindex = GALLERY_HEADER

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    tagged_examples = example_groups(src_dir)
    tagged_examples = order_examples(tagged_examples)

    computation_times = []
    build_target_dir = os.path.relpath(target_dir, gallery_conf['src_dir'])

    seen = set()
    tmp_dir = tempfile.mkdtemp()

    for tag, examples in tagged_examples.items():
        sorted_listdir = examples

        entries_text = []
        iterator = sphinx_compatibility.status_iterator(
            sorted_listdir,
            'Generating gallery for %s ' % tag,
            length=len(sorted_listdir))
        for fname in iterator:
            write_example(os.path.join(src_dir, fname), tmp_dir)
            amount_of_code, time_elapsed = generate_file_rst(
                fname, target_dir, tmp_dir, gallery_conf)

            if fname not in seen:
                seen.add(fname)
                computation_times.append((time_elapsed, fname))

            new_fname = os.path.join(src_dir, fname)
            intro = extract_intro(new_fname)
            this_entry = _thumbnail_div(build_target_dir, fname, intro) + textwrap.dedent("""

                .. toctree::
                   :hidden:

                   /%s

                   """) % os.path.join(build_target_dir, fname[:-3]).replace(os.sep, '/')  # noqa: E501

            entries_text.append((amount_of_code, this_entry))

            if gallery_conf['backreferences_dir']:
                write_backreferences(seen_backrefs, gallery_conf,
                                     target_dir, fname, intro)

        # sort to have the smallest entries in the beginning
        entries_text.sort()

        fhindex += textwrap.dedent("""

        {tag}
        {tag_underline}

        .. container:: gallery_images

        """.format(tag=tag, tag_underline='-' * len(tag)))

        for _, entry_text in entries_text:
            fhindex += '\n    '.join(entry_text.split('\n'))

        # clear at the end of the section
        fhindex += """.. raw:: html\n
        <div style='clear:both'></div>\n\n"""

    # Tidy up the temp directory
    shutil.rmtree(tmp_dir)

    return fhindex, computation_times


# Monkey-patch sphinx_gallery to handle cartopy's example format.
sphinx_gallery.gen_rst.generate_dir_rst = generate_dir_rst
sphinx_gallery.gen_gallery.generate_dir_rst = generate_dir_rst


def setup(app):
    pass
