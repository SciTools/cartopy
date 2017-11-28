# (C) British Crown Copyright 2011 - 2017, Met Office
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

from bs4 import BeautifulSoup as bs

path_html = 'build/html/gallery/index.html'
with open(path_html, 'r') as ff:
    html = bs(ff.read(), "lxml")

    # Search for thumbnails that link to __init__ files
    thumbs = html.find_all('div', attrs={'class': 'sphx-glr-thumbcontainer'})
    for div in thumbs:
        links = div.find_all('a')

        # If you find any, add a class that will make it invisible
        any_init = any('__init__' in ilink.attrs['href'] for ilink in links)
        if any_init:
            div.attrs['class'] += ['dontshow']

# Now update the HTML
with open(path_html, 'w') as ff:
    out_html = html.prettify(html.original_encoding)
    ff.write(out_html)
