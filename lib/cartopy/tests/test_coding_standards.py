# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from __future__ import (absolute_import, division, print_function)

from datetime import datetime
from fnmatch import fnmatch
import io
from itertools import chain
import os
import re
import subprocess

import pytest

import cartopy


LICENSE_TEMPLATE = """
# (C) British Crown Copyright {YEARS}, Met Office
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
# along with cartopy.  If not, see <https://www.gnu.org/licenses/>.""".strip()


LICENSE_RE_PATTERN = re.escape(LICENSE_TEMPLATE).replace(r'\{YEARS\}', '(.*?)')
# Add shebang possibility or C comment starter to the LICENSE_RE_PATTERN
SHEBANG_PATTERN = r'((\#\!.*|\/\*)\n)?'
LICENSE_RE = re.compile(SHEBANG_PATTERN + LICENSE_RE_PATTERN, re.MULTILINE)


LICENSE_TEMPLATE_v2 = """
# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
""".strip()
LICENSE_RE_PATTERN_v2 = re.escape(LICENSE_TEMPLATE_v2)
LICENSE_RE_v2 = re.compile(SHEBANG_PATTERN + LICENSE_RE_PATTERN_v2,
                           re.MULTILINE)


# Guess cartopy repo directory of cartopy - realpath is used to mitigate
# against Python finding the cartopy package via a symlink.
CARTOPY_DIR = os.path.realpath(os.path.dirname(cartopy.__file__))
REPO_DIR = os.getenv('CARTOPY_GIT_DIR',
                     os.path.dirname(os.path.dirname(CARTOPY_DIR)))


class TestLicenseHeaders(object):
    @staticmethod
    def years_of_license_in_file(content, fname):
        """
        Using :data:`LICENSE_RE` look for the years defined in the license
        header of the given file handle.

        If the license cannot be found in the given content, None will be
        returned, else a tuple of (start_year, end_year) will be returned.

        """
        license_matches = LICENSE_RE.match(content)
        if not license_matches:
            # no license found in file.
            return None

        years = license_matches.groups()[-1]
        if len(years) == 4:
            start_year = end_year = int(years)
        elif len(years) == 11:
            start_year, end_year = int(years[:4]), int(years[7:])
        else:
            raise ValueError("Unexpected year(s) string in {}'s copyright "
                             "notice: {!r}".format(fname, years))
        return (start_year, end_year)

    @staticmethod
    def last_change_by_fname():
        """
        Return a dictionary of all the files under git which maps to
        the datetime of their last modification in the git history.

        .. note::

            This function raises a ValueError if the repo root does
            not have a ".git" folder. If git is not installed on the system,
            or cannot be found by subprocess, an IOError may also be raised.

        """
        # Check the ".git" folder exists at the repo dir.
        if not os.path.isdir(os.path.join(REPO_DIR, '.git')):
            raise ValueError('{} is not a git repository.'.format(REPO_DIR))

        # Call "git whatchanged" to get the details of all the files and when
        # they were last changed.
        output = subprocess.check_output(['git', 'ls-tree', '-z', '-r',
                                          '--name-only', 'HEAD'],
                                         cwd=REPO_DIR)
        output = output.rstrip(b'\0').split(b'\0')
        res = {}
        for fname in output:
            fname = fname.decode()
            dt = subprocess.check_output(['git', 'log', '-1', '--pretty=%ct',
                                          '--', fname],
                                         cwd=REPO_DIR)
            dt = datetime.fromtimestamp(int(dt))
            res[fname] = dt

        return res

    def test_license_headers(self):
        exclude_patterns = ('build/*',
                            'dist/*',
                            'docs/build/*',
                            'docs/source/examples/*.py',
                            'docs/source/sphinxext/*.py',
                            'lib/cartopy/examples/*.py',
                            'lib/cartopy/_version.py',
                            'versioneer.py',
                            )
        try:
            last_change_by_fname = self.last_change_by_fname()
        except ValueError as e:
            # Caught the case where this is not a git repo.
            return pytest.skip('cartopy installation did not look like a git '
                               'repo: ' + str(e))

        failed = False
        for fname, last_change in sorted(last_change_by_fname.items()):
            full_fname = os.path.join(REPO_DIR, fname)
            root, ext = os.path.splitext(full_fname)
            if ext in ('.py', '.pyx', '.c', '.cpp', '.h') and \
                    os.path.isfile(full_fname) and \
                    not any(fnmatch(fname, pat) for pat in exclude_patterns):

                is_empty = os.path.getsize(full_fname) == 0

                with io.open(full_fname, encoding='utf-8') as fh:
                    content = fh.read()

                is_yearless_license = bool(LICENSE_RE_v2.match(content))
                years = TestLicenseHeaders.years_of_license_in_file(
                    content, full_fname)

                if is_empty:
                    # Allow completely empty files (e.g. ``__init__.py``)
                    pass
                elif is_yearless_license:
                    # Allow new style license (v2).
                    pass

                # What is left is the old-style (pre 2019) header.
                elif years is None:
                    print('The file {} has no valid header license and '
                          'has not been excluded from the license header '
                          'test.'.format(fname))
                    failed = True
                elif last_change.year > years[1]:
                    print('The file header at {} is out of date. The last'
                          ' commit was in {}, but the copyright states it'
                          ' was {}.'.format(fname, last_change.year,
                                            years[1]))
                    failed = True

        if failed:
            raise ValueError('There were license header failures. See stdout.')


class TestFutureImports(object):
    excluded = (
        '*/cartopy/examples/*.py',
        '*/docs/source/examples/*.py',
        '*/cartopy/_crs.py',   # A file created by setuptools for so loading.
        '*/cartopy/trace.py',  # Ditto.
        '*/cartopy/geodesic/_geodesic.py',
        '*/cartopy/_version.py',
    )

    future_imports_pattern = re.compile(
        r"^from __future__ import \(absolute_import,\s*division,\s*"
        r"print_function(,\s*unicode_literals)?\)$",
        flags=re.MULTILINE)

    def test_future_imports(self):
        # Tests that every single Python file includes the appropriate
        # __future__ import to enforce consistent behaviour.
        check_paths = [CARTOPY_DIR]

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

                is_empty = os.path.getsize(full_fname) == 0

                with io.open(full_fname, "r", encoding='utf-8') as fh:
                    content = fh.read()

                has_future_import = re.search(
                    self.future_imports_pattern, content) is not None

                if is_empty:
                    pass
                elif not has_future_import:
                    print('The file {} has no valid __future__ imports '
                          'and has not been excluded from the imports '
                          'test.'.format(full_fname))
                    failed = True

        if failed:
            raise ValueError('There were __future__ import check failures. '
                             'See stdout.')
