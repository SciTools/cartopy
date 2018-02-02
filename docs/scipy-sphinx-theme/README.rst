scipy-sphinx-theme
==================

`Sphinx <http://sphinx-doc.org>`__ theme for `Scipy <http://scipy.org>`__.


Theme options
-------------

The theme takes the followin options in the `html_options`
configuration variable:

- ``edit_link``

  ``True`` or ``False``. Determines if an "edit this page" link is displayed
  in the left sidebar.

- ``rootlinks``

  List of tuples ``(url, link_name)`` to show in the beginning of the
  breadcrumb list on the top left. You can override it by defining an
  `edit_link` block in ``searchbox.html``.

- ``sidebar``

  One of ``"left"``, ``"right"``, ``"none"``.  Defines where the sidebar
  should appear.

- ``scipy_org_logo``

  ``True`` or ``False``. Whether to plaster the scipy.org logo on top.

  You can use your own logo by overriding the :attr:`layout.html:header`
  block.

- ``navigation_links``

  ``True`` or ``False``. Whether to display "next", "prev", "index", etc.
  links.

The following blocks are defined:

- ``layout.html:header``
   
  Block at the top of the page, for logo etc.

- ``searchbox.html:edit_link``

  Edit link HTML code to paste in the left sidebar, if `edit_link` is enabled.
