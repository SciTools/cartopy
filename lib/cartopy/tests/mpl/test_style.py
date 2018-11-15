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

import pytest

from ...mpl import style


d = dict


@pytest.mark.parametrize(
    ('styles', 'expected'),
    [([], {}),
     ([{}, {}, {}], {}),
     ([{}, d(a=2), d(a=1)], d(a=1)),
     ([d(fc='red')], d(facecolor='red')),
     ([d(fc='red', color='blue')], d(facecolor='blue', edgecolor='blue')),
     ([d(fc='red', facecolor='blue')], d(facecolor='blue')),
     ([d(color='red')],
      d(edgecolor='red', facecolor='red')),
     ([d(edgecolor='blue'), d(color='red')],
      d(edgecolor='red', facecolor='red')),
     ([d(edgecolor='blue'), d(color='red')],
      d(edgecolor='red', facecolor='red')),
     ([d(color='blue'), d(edgecolor='red')],
      d(edgecolor='red', facecolor='blue')),
     # Even if you set an edgecolor, color should trump it.
     ([d(color='blue'), d(edgecolor='red', color='yellow')],
      d(edgecolor='yellow', facecolor='yellow')),
     # Support for 'never' being honoured.
     ([d(facecolor='never'), d(color='yellow')],
      d(edgecolor='yellow', facecolor='never')),
     ([d(lw=1, linewidth=2)],
      d(linewidth=2)),
     ([d(lw=1, linewidth=2), d(lw=3)],
      d(linewidth=3)),
     ]
)
def test_merge(styles, expected):
    merged_style = style.merge(*styles)
    assert merged_style == expected


def test_merge_warning():
    with pytest.warns(UserWarning, match=r'defined as \"never\"') as record:
        style.merge({'facecolor': 'never'}, {'fc': 'red'})
    assert len(record) == 1


@pytest.mark.parametrize(
    ('style_d', 'expected'),
    [
     # Support for 'never' being honoured.
     (d(facecolor='never', edgecolor='yellow'),
      d(edgecolor='yellow', facecolor='none')),
    ])
def test_finalize(style_d, expected):
    assert style.finalize(style_d) == expected
    # Double check we are updating in-place
    assert style_d == expected
