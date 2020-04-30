# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

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
                     tolerance=4.5 if MPL_VERSION < '2' else 0.5)
def test_global_map():
    import cartopy.examples.global_map as example
    example.main()


if MPL_VERSION < '2':
    contour_labels_tolerance = 7.5
elif MPL_VERSION <= '2.0.2':
    contour_labels_tolerance = 1.24
elif MPL_VERSION <= '2.1.2':
    contour_labels_tolerance = 0.63
else:
    contour_labels_tolerance = 0


@pytest.mark.natural_earth
@ExampleImageTesting(['contour_label'], tolerance=contour_labels_tolerance)
def test_contour_label():
    import cartopy.examples.contour_labels as example
    example.main()
