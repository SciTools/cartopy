# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from unittest import mock

from matplotlib.testing.decorators import cleanup
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import pytest

import cartopy.crs as ccrs

from cartopy.mpl.scalebar import fancy_scalebar


def test_scalebar():
    """Test that we do not trigger matplotlib's line subslice optimization."""
    # This behavior caused lines with > 1000 points and
    # sorted data to disappear

    fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.Mercator()})
    
    projections = [ccrs.Mercator(), ccrs.PlateCarree()]
    
    axes = axes.ravel()
    
    for proj,  ax in zip(projections, axes):
        ax.projection = projection
        fancy_scalebar(ax, 
                       location=(0.5, 0.2), 
                       length=10000,
                       metres_per_unit=1000, 
                       unit_name='km',
                       tol=0.01, 
                           angle=0,
                           max_stripes=5,
                           fontsize= 8,
                           dy = 0.05)

