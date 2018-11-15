# (C) British Crown Copyright 2018, Met Office
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

"""
Handles matplotlib styling in a single consistent place.

"""
import warnings


# Define the matplotlib style aliases that cartopy can expand.
# Note: This should not contain the plural aliases
# (e.g. linewidths -> linewidth).
# This is an intended duplication of
# https://github.com/matplotlib/matplotlib/blob/\
#   2d2dab511d22b6cc9c812cfbcca6df3f9bf3094a/lib/matplotlib/patches.py#L20-L26
# Duplication intended to simplify readability, given the small number of
# aliases.
_ALIASES = {
    'lw': 'linewidth',
    'ls': 'linestyle',
    'fc': 'facecolor',
    'ec': 'edgecolor',
}


def merge(*style_dicts):
    """
    Merge together multiple matplotlib style dictionaries in a predictable way

    The approach taken is:

        For each style:
            * Expand aliases, such as "lw" -> "linewidth", but always prefer
              the full form if over-specified (i.e. lw AND linewidth
              are both set)
            * "color" overwrites "facecolor" and "edgecolor" (as per
              matplotlib), UNLESS facecolor == "never", which will be expanded
              at finalization to 'none'

    >>> merge({'lw': 1, "edgecolor": "black", "facecolor": "never"},
    ...       {"linewidth": 2, "color": "gray"})
    {'edgecolor': 'gray', 'facecolor': 'none', 'linewidth': 2}

    """
    style = {}
    facecolor = None

    for orig_style in style_dicts:
        this_style = orig_style.copy()

        for alias_from, alias_to in _ALIASES.items():
            alias = this_style.pop(alias_from, None)
            if alias:
                # n.b. alias_from doesn't trump alias_to
                # (e.g. 'lw' doesn't trump 'linewidth').
                this_style.setdefault(alias_to, alias)

        color = this_style.pop('color', None)
        if color:
            this_style['edgecolor'] = color
            this_style['facecolor'] = color

        if facecolor == 'never':
            this_style.pop('facecolor', None)
            if 'fc' in orig_style or 'facecolor' in orig_style:
                warnings.warn('facecolor may not be set, as it has been '
                              'defined as "never".')
        else:
            facecolor = this_style.get('facecolor', facecolor)

        # Push the remainder of the style into the merged style.
        style.update(this_style)

    return style


def finalize(style):
    """
    Update the given matplotlib style according to cartopy's style rules.

    Rules:

        1. A facecolor of 'never' is replaced with 'none'.

    """
    # Expand 'never' to 'none' if we have it.
    facecolor = style.get('facecolor', None)
    if facecolor == 'never':
        style['facecolor'] = 'none'
    return style
