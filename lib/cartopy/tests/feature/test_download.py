# Copyright Crown and Cartopy Contributors
#
# This file is part of Cartopy and is released under the BSD 3-clause license.
# See LICENSE in the root of the repository for full licensing details.

from cartopy.feature.download.__main__ import download_features


def test_cultural_extra_includes_50m_states_provinces_lines(capsys):
    download_features(['cultural-extra'])

    urls = capsys.readouterr().out.splitlines()
    assert ('URL: https://naturalearth.s3.amazonaws.com/50m_cultural/'
            'ne_50m_admin_1_states_provinces_lines.zip') in urls
