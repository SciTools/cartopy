import glob
import os
import unittest

import pep8

import cartopy


class TestCodeFormat(unittest.TestCase):
    def test_pep8_conformance(self):
        pep8style = pep8.StyleGuide(quiet=False)
        pep8style.options.exclude.extend(['gshhs.py', 'examples'])
        result = pep8style.check_files([os.path.dirname(cartopy.__file__)])
        self.assertEqual(result.total_errors, 0, "Found code syntax errors (and warnings).")


if __name__ == '__main__':
    unittest.main()