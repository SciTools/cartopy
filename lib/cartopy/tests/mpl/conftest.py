# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.
from contextlib import ExitStack

import matplotlib.pyplot as plt
import pytest

from cartopy.mpl import _MPL_311


@pytest.fixture(autouse=True)
def mpl_test_cleanup(request):
    """Run tests in a context manager and close figures after each test."""
    with ExitStack() as stack:
        # At exit, close all open figures and switch backend back to original.
        stack.callback(plt.switch_backend, plt.get_backend())
        stack.callback(plt.close, 'all')

        # Run each test in a context manager so that state does not leak out
        plt.switch_backend('Agg')
        stack.enter_context(plt.rc_context())
        yield


def pytest_itemcollected(item):
    mpl_marker = item.get_closest_marker('mpl_image_compare')
    if mpl_marker is None:
        return

    # Matches old ImageTesting class default tolerance.
    mpl_marker.kwargs.setdefault('tolerance', 0.5)

    for path in item.fspath.parts(reverse=True):
        if path.basename == 'cartopy':
            return
        elif path.basename == 'tests':
            subdir = item.fspath.relto(path)[:-len(item.fspath.ext)]
            mpl_marker.kwargs.setdefault('baseline_dir',
                                         f'baseline_images/{subdir}')
            break


@pytest.fixture
def text_placeholders(monkeypatch):
    """
    Replace texts with placeholder rectangles.

    The rectangle size only depends on the font size and the number of characters. It is
    thus insensitive to font properties and rendering details. This should be used for
    tests that depend on text geometries but not the actual text rendering, e.g. layout
    tests.
    """
    from matplotlib.patches import Rectangle

    def patched_get_sfnt_table(font, name):
        """
        Replace ``FT2Font.get_sfnt_table`` with empty results.

        This forces ``Text._get_layout`` to fall back to
        ``get_text_width_height_descent``, which produces results from the patch below.
        """
        return None

    def patched_get_text_metrics_with_cache(renderer, text, fontprop, ismath, dpi):
        """
        Replace ``_get_text_metrics_with_cache`` with fixed results.

        The usual ``renderer.get_text_width_height_descent`` would depend on font
        metrics; instead the fixed results are based on font size and the length of the
        string only.
        """
        # While get_window_extent returns pixels and font size is in points, font size
        # includes ascenders and descenders. Leaving out this factor and setting
        # descent=0 ends up with a box that is relatively close to DejaVu Sans.
        height = fontprop.get_size()
        width = len(text) * height / 1.618  # Golden ratio for character size.
        descent = 0
        return width, height, descent

    def patched_text_draw(self, renderer):
        """
        Replace ``Text.draw`` with a fixed bounding box Rectangle.

        The bounding box corresponds to ``Text.get_window_extent``, which ultimately
        depends on the above patched ``_get_text_metrics_with_cache``.
        """
        if renderer is not None:
            self._renderer = renderer
        if not self.get_visible():
            return
        if self.get_text() == '':
            return
        bbox = self.get_window_extent()
        # Keep the placeholder aligned with rotated labels as well.
        rect = Rectangle(
            bbox.p0, bbox.width, bbox.height,
            angle=self.get_rotation(), rotation_point='center',
            facecolor=self.get_color(), edgecolor='none',
        )
        rect.draw(renderer)

    if _MPL_311:
        monkeypatch.setattr('matplotlib.ft2font.FT2Font.get_sfnt_table',
                            patched_get_sfnt_table)
    monkeypatch.setattr('matplotlib.text._get_text_metrics_with_cache',
                        patched_get_text_metrics_with_cache)
    monkeypatch.setattr('matplotlib.text.Text.draw', patched_text_draw)
