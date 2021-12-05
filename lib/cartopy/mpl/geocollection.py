# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
from matplotlib.collections import QuadMesh
import numpy as np


class GeoQuadMesh(QuadMesh):
    """
    A QuadMesh designed to help handle the case when the mesh is wrapped.

    """
    # No __init__ method here - most of the time a GeoQuadMesh will
    # come from GeoAxes.pcolormesh. These methods morph a QuadMesh by
    # fiddling with instance.__class__.

    def get_array(self):
        # Retrieve the array - use copy to avoid any chance of overwrite
        A = super(QuadMesh, self).get_array().copy()
        # If the input array has a mask, retrieve the associated data
        if hasattr(self, '_wrapped_mask'):
            A[self._wrapped_mask] = self._wrapped_collection_fix.get_array()
        return A

    def set_array(self, A):
        # raise right away if A is 2-dimensional.
        if A.ndim > 1:
            raise ValueError('Collections can only map rank 1 arrays. '
                             'You likely want to call with a flattened array '
                             'using collection.set_array(A.ravel()) instead.')

        # Only use the mask attribute if it is there.
        if hasattr(self, '_wrapped_mask'):
            # Update the pcolor data with the wrapped masked data
            self._wrapped_collection_fix.set_array(A[self._wrapped_mask])
            # If the input array was a masked array, keep that data masked
            if hasattr(A, 'mask'):
                A = np.ma.array(A, mask=self._wrapped_mask | A.mask)
            else:
                A = np.ma.array(A, mask=self._wrapped_mask)

        # Now that we have prepared the collection data, call on
        # through to the underlying implementation.
        super(QuadMesh, self).set_array(A)

    def set_clim(self, vmin=None, vmax=None):
        # Update _wrapped_collection_fix color limits if it is there.
        if hasattr(self, '_wrapped_collection_fix'):
            self._wrapped_collection_fix.set_clim(vmin, vmax)

        # Update color limits for the rest of the cells.
        super().set_clim(vmin, vmax)

    def get_datalim(self, transData):
        # Return the corners that were calculated in
        # the pcolormesh routine.
        return self._corners
