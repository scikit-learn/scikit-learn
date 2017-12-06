# Documentation for scikit-learn

This directory contains the full manual and web site as displayed at
http://scikit-learn.org.

## Building manually

Building the website requires the sphinx, sphinx-gallery and matplotlib
packages:

    pip install sphinx numpydoc sphinx-gallery matplotlib

It also requires having the version of scikit-learn installed that corresponds
to the documentation, e.g.:

    pip install --editable ..

To generate the full web site, including the example gallery (this might take a
while):

    make html

Or, if you'd rather not build the example gallery:

    make html-noplot

That should create all the doc in directory `_build/html`.  Set the environment
variable `USE_MATHJAX=1` if you intend to view the documentation in an online
setting.

To build the PDF manual, run

    make latexpdf

Make sure you first have the correct version of scikit-learn

## Hosting and automatic builds

The website is hosted at github, but should rarely be updated manually
by pushing to the https://github.com/scikit-learn/scikit-learn.github.io repository.

Most updates can be made by pushing to master (for /dev) or a release branch
like 0.99.X, from which Circle CI builds and uploads documentation. (See the
Developer Documentation for further details.)
