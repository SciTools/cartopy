.. _contribute:

Contributing to Cartopy
=======================

Cartopy is an open-source project that welcomes contributions from anyone. More information
about the copyright and license agreement can be found on the `SciTools Governance page <https://scitools.org.uk/organisation.html>`_.
A good place to start contributing is by fixing typos, adding examples, or improving the explanations
of the functions.

To contribute new features or bug fixes, you can create a pull request with your code changes
on GitHub. The easiest way to get started contributing code is to `Fork Cartopy <https://github.com/scitools/cartopy/fork>`_
to your GitHub account and then clone the repository to your
local system. After that you can install the conda development environment included with the source code
to install all of the required dependencies

.. code::

  # Clone your forked version of Cartopy
  git clone git@github.com:YOUR-USERNAME/cartopy.git

  # Set up the development environment and install the dependencies
  cd cartopy
  conda env create -f environment.yml
  conda activate cartopy-dev

  # Install cartopy
  pip install -e .

This will install all of the required dependencies for compiling the code, running the tests, and
building the documentation. Remember to activate the `cartopy-dev` environment.

A development and test workflow may look something like the following code block.

.. code::

  # Activate your development environment
  conda activate cartopy-dev

  # Checkout a new feature branch
  git checkout -b my-cool-feature

  # Make your code changes here, including tests for the changes!
  # Add and commit the changes with a useful commit message
  # describing what changes you made
  git add lib/cartopy/my-changed-file.py
  git commit

  # Test your code on 4 processors
  pytest -n 4 lib/cartopy

  # push your changes to your GitHub account
  git push --set-upstream origin my-cool-feature

  # Go to GitHub and make a Pull Request to Cartopy
