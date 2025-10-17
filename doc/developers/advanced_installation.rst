
.. _advanced-installation:

==================================================
Installing the development version of scikit-learn
==================================================

.. _install_nightly_builds:

Installing nightly builds
=========================

The continuous integration servers of the scikit-learn project build, test
and upload wheel packages for the most recent Python version on a nightly
basis.

Installing a nightly build is the quickest way to:

- try a new feature that will be shipped in the next release (that is, a
  feature from a pull-request that was recently merged to the main branch);

- check whether a bug you encountered has been fixed since the last release.

You can install the nightly build of scikit-learn using the `scientific-python-nightly-wheels`
index from the PyPI registry of `anaconda.org`:

.. prompt:: bash $

  pip install --pre --extra-index https://pypi.anaconda.org/scientific-python-nightly-wheels/simple scikit-learn

Note that first uninstalling scikit-learn might be required to be able to
install nightly builds of scikit-learn.
