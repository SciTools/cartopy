# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.


def pytest_configure(config):
    # Register additional markers.
    config.addinivalue_line('markers',
                            'natural_earth: mark tests that use Natural Earth '
                            'data, and the network, if not cached.')
    config.addinivalue_line('markers',
                            'network: mark tests that use the network.')
