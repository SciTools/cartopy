# (C) British Crown Copyright 2011 - 2018, Met Office
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

import matplotlib.pyplot as plt
import pytest

from cartopy.tests.mpl import MPL_VERSION, ImageTesting


class ExampleImageTesting(ImageTesting):
    """Subclasses ImageTesting to nullify the plt.show commands."""
    def __call__(self, test_func):
        fn = ImageTesting.__call__(self, test_func)

        def new_fn(*args, **kwargs):
            try:
                show = plt.show
                plt.show = lambda *args, **kwargs: None
                r = fn(*args, **kwargs)
            finally:
                plt.show = show
            return r

        new_fn.__name__ = fn.__name__
        return new_fn


@pytest.mark.natural_earth
@ExampleImageTesting(['global_map'],
                     tolerance=4 if MPL_VERSION < '2' else 0)
def test_global_map():
    import cartopy.examples.global_map as c
    c.main()
