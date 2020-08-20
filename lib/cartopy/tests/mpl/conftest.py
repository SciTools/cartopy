# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.
import matplotlib.pyplot as plt
import pytest


@pytest.fixture(autouse=True)
def mpl_test_cleanup(request):
    """Runs tests in a context manager and closes figures after each test"""
    try:
        # Run each test in a context manager so that state does not leak out
        orig_backend = plt.get_backend()
        plt.switch_backend('Agg')
        with plt.rc_context():
            yield
    finally:
        # Closes all open figures and switches backend back to original
        plt.switch_backend(orig_backend)
