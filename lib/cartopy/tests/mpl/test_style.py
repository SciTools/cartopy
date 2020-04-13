# (C) British Crown Copyright 2018 - 2019, Met Office
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

from cartopy.mpl import style


@pytest.mark.parametrize(
    ('styles', 'expected'),
    [([], {}),
     ([{}, {}, {}], {}),
     ([{}, dict(a=2), dict(a=1)], dict(a=1)),
     ([dict(fc='red')], dict(facecolor='red')),
     ([dict(fc='red', color='blue')],
      dict(facecolor='blue', edgecolor='blue')),
     ([dict(fc='red', facecolor='blue')], dict(facecolor='blue')),
     ([dict(color='red')],
      dict(edgecolor='red', facecolor='red')),
     ([dict(edgecolor='blue'), dict(color='red')],
      dict(edgecolor='red', facecolor='red')),
     ([dict(edgecolor='blue'), dict(color='red')],
      dict(edgecolor='red', facecolor='red')),
     ([dict(color='blue'), dict(edgecolor='red')],
      dict(edgecolor='red', facecolor='blue')),
     # Even if you set an edgecolor, color should trump it.
     ([dict(color='blue'), dict(edgecolor='red', color='yellow')],
      dict(edgecolor='yellow', facecolor='yellow')),
     # Support for 'never' being honoured.
     ([dict(facecolor='never'), dict(color='yellow')],
      dict(edgecolor='yellow', facecolor='never')),
     ([dict(lw=1, linewidth=2)], dict(linewidth=2)),
     ([dict(lw=1, linewidth=2), dict(lw=3)], dict(linewidth=3)),
     ([dict(color=None), dict(facecolor='red')],
      dict(facecolor='red', edgecolor=None)),
     ([dict(linewidth=1), dict(lw=None)], dict(linewidth=None)),
     ]
)
def test_merge(styles, expected):
    merged_style = style.merge(*styles)
    assert merged_style == expected


@pytest.mark.parametrize(
    ('case', 'should_warn'),
    [[{'fc': 'red'}, True], [{'fc': 'NoNe'}, False], [{'fc': 1}, True]])
def test_merge_warning(case, should_warn):
    warn_type = UserWarning if should_warn else None
    with pytest.warns(warn_type, match=r'defined as \"never\"') as record:
        style.merge({'facecolor': 'never'}, case)
    assert len(record) == (1 if should_warn else 0)


@pytest.mark.parametrize(
    ('style_d', 'expected'),
    [
        # Support for 'never' being honoured.
        (dict(facecolor='never', edgecolor='yellow'),
         dict(edgecolor='yellow', facecolor='none')),
    ])
def test_finalize(style_d, expected):
    assert style.finalize(style_d) == expected
    # Double check we are updating in-place
    assert style_d == expected
