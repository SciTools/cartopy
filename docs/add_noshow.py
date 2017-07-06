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
