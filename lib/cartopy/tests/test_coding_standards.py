# (C) British Crown Copyright 2012 - 2016, Met Office
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

from datetime import datetime
from fnmatch import fnmatch
from itertools import chain
import os
import re
import subprocess
import unittest

import pep8

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


LICENSE_RE_PATTERN = re.escape(LICENSE_TEMPLATE).replace('\{YEARS\}', '(.*?)')
# Add shebang possibility or C comment starter to the LICENSE_RE_PATTERN
LICENSE_RE_PATTERN = r'((\#\!.*|\/\*)\n)?' + LICENSE_RE_PATTERN
LICENSE_RE = re.compile(LICENSE_RE_PATTERN, re.MULTILINE)


# Guess cartopy repo directory of cartopy - realpath is used to mitigate
# against Python finding the cartopy package via a symlink.
CARTOPY_DIR = os.path.realpath(os.path.dirname(cartopy.__file__))
REPO_DIR = os.getenv('CARTOPY_GIT_DIR',
                     os.path.dirname(os.path.dirname(CARTOPY_DIR)))


class TestLicenseHeaders(unittest.TestCase):
    @staticmethod
    def years_of_license_in_file(fh):
        """
        Using :data:`LICENSE_RE` look for the years defined in the license
        header of the given file handle.

        If the license cannot be found in the given fh, None will be returned,
        else a tuple of (start_year, end_year) will be returned.

        """
        license_matches = LICENSE_RE.match(fh.read())
        if not license_matches:
            # no license found in file.
            return None

        years = license_matches.groups()[-1]
        if len(years) == 4:
            start_year = end_year = int(years)
        elif len(years) == 11:
            start_year, end_year = int(years[:4]), int(years[7:])
        else:
            fname = getattr(fh, 'name', 'unknown filename')
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
        output = subprocess.check_output(['git', 'ls-tree', '-r',
                                          '--name-only', 'HEAD'],
                                         cwd=REPO_DIR)
        output = output.decode().split('\n')
        res = {}
        for fname in output:
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
                            'lib/cartopy/examples/*.py')

        try:
            last_change_by_fname = self.last_change_by_fname()
        except ValueError as e:
            # Caught the case where this is not a git repo.
            return self.skipTest('cartopy installation did not look like a '
                                 'git repo: ' + str(e))

        failed = False
        for fname, last_change in sorted(last_change_by_fname.items()):
            full_fname = os.path.join(REPO_DIR, fname)
            root, ext = os.path.splitext(full_fname)
            if ext in ('.py', '.pyx', '.c', '.cpp', '.h') and \
                    os.path.isfile(full_fname) and \
                    not any(fnmatch(fname, pat) for pat in exclude_patterns):
                with open(full_fname) as fh:
                    years = TestLicenseHeaders.years_of_license_in_file(fh)
                    if years is None:
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
        pep8style.options.exclude.extend(['trace.py', '_crs.py',
                                          '*/cartopy/geodesic/_geodesic.py'])

        # Ignore E402 module level import not at top of file
        pep8style.options.ignore += ('E402', )

        # allow users to add their own exclude list
        extra_exclude_file = os.path.join(os.path.dirname(__file__),
                                          '.pep8_test_exclude.txt')
        if os.path.exists(extra_exclude_file):
            with open(extra_exclude_file, 'r') as fh:
                extra_exclude = [line.strip() for line in fh if line.strip()]
            pep8style.options.exclude.extend(extra_exclude)

        result = pep8style.check_files([CARTOPY_DIR])
        self.assertEqual(result.total_errors, 0, "Found code syntax "
                                                 "errors (and warnings).")


class TestFutureImports(unittest.TestCase):
    excluded = (
        '*/cartopy/examples/*.py',
        '*/docs/source/examples/*.py',
        '*/cartopy/_crs.py',   # A file created by setuptools for so loading.
        '*/cartopy/trace.py',  # Ditto.
        '*/cartopy/geodesic/_geodesic.py',
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
