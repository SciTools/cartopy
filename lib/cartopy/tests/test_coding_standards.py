# (C) British Crown Copyright 2012 - 2014, Met Office
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

from __future__ import (absolute_import, division, print_function)

from fnmatch import fnmatch
from itertools import chain
import os
import re
import unittest

import pep8

import cartopy


class TestCodeFormat(unittest.TestCase):
    def test_pep8_conformance(self):
        # Tests the cartopy codebase against the "pep8" tool.
        #
        # Users can add their own excluded files (should files exist in the
        # local directory which is not in the repository) by adding a
        # ".pep8_test_exclude.txt" file in the same directory as this test.
        # The file should be a line separated list of filenames/directories
        # as can be passed to the "pep8" tool's exclude list.
        pep8style = pep8.StyleGuide(quiet=False)
        pep8style.options.exclude.extend(['trace.py', '_crs.py'])

        # allow users to add their own exclude list
        extra_exclude_file = os.path.join(os.path.dirname(__file__),
                                          '.pep8_test_exclude.txt')
        if os.path.exists(extra_exclude_file):
            with open(extra_exclude_file, 'r') as fh:
                extra_exclude = [line.strip() for line in fh if line.strip()]
            pep8style.options.exclude.extend(extra_exclude)

        result = pep8style.check_files([os.path.dirname(cartopy.__file__)])
        self.assertEqual(result.total_errors, 0, "Found code syntax "
                                                 "errors (and warnings).")


class TestFutureImports(unittest.TestCase):
    excluded = (
        '*/cartopy/examples/*.py',
        '*/docs/source/examples/*.py',
        '*/cartopy/_crs.py',   # A file created by setuptools for so loading.
        '*/cartopy/trace.py',  # Ditto.
    )

    future_imports_pattern = re.compile(
        r"^from __future__ import \(absolute_import,\s*division,\s*"
        r"print_function(,\s*unicode_literals)?\)$",
        flags=re.MULTILINE)

    def test_future_imports(self):
        # Tests that every single Python file includes the appropriate
        # __future__ import to enforce consistent behaviour.
        check_paths = [os.path.dirname(cartopy.__file__)]

        failed = False
        for dirpath, _, files in chain.from_iterable(os.walk(path)
                                                     for path in check_paths):
            for fname in files:
                full_fname = os.path.join(dirpath, fname)
                if not full_fname.endswith('.py'):
                    continue
                if not os.path.isfile(full_fname):
                    continue
                if any(fnmatch(full_fname, pat) for pat in self.excluded):
                    continue

                with open(full_fname, "r") as fh:
                    content = fh.read()

                    if re.search(self.future_imports_pattern, content) is None:
                        print('The file {} has no valid __future__ imports '
                              'and has not been excluded from the imports '
                              'test.'.format(full_fname))
                        failed = True

        if failed:
            raise ValueError('There were __future__ import check failures. '
                             'See stdout.')


if __name__ == '__main__':
    unittest.main()
