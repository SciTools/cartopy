
import pytest
from ...mpl import style

d = dict

@pytest.mark.parametrize(('styles', 'expected'),
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
    ])
def test_merge(styles, expected):
    merged_style = style.merge(*styles)
    assert merged_style == expected


@pytest.mark.parametrize(('style_d', 'expected'),
   [
    # Support for 'never' being honoured.
    (d(facecolor='never', edgecolor='yellow'),
     d(edgecolor='yellow', facecolor='none')),
   ])
def test_finalize(style_d, expected):
    assert style.finalize(style_d) == expected
    # Double check we are updating in-place
    assert style_d == expected
