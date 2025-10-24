
.. _advanced-installation:

==================================================
More advanced installation / Troubleshooting
==================================================

Here, you find some more advanced notes and troubleshooting tips related to
:ref:`install_bleeding_edge`.

.. _openMP_notes:

Notes on OpenMP
===============

Even though the default C compiler on macOS (Apple clang) is confusingly aliased
as `/usr/bin/gcc`, it does not directly support OpenMP.

.. note::

  If OpenMP is not supported by the compiler, the build will be done with
  OpenMP functionalities disabled. This is not recommended since it will force
  some estimators to run in sequential mode instead of leveraging thread-based
  parallelism. Setting the ``SKLEARN_FAIL_NO_OPENMP`` environment variable
  (before cythonization) will force the build to fail if OpenMP is not
  supported.

To check if `scikit-learn` has been built correctly with OpenMP, run

.. prompt:: bash $

  python -c "import sklearn; sklearn.show_versions()"

and check if it contains `Built with OpenMP: True`.

When using conda on Mac, you can also check that the custom compilers
are properly installed from conda-forge using the following command:

.. prompt:: bash $

    conda list

which should include ``compilers`` and ``llvm-openmp``.

The compilers meta-package will automatically set custom environment
variables:

.. prompt:: bash $

    echo $CC
    echo $CXX
    echo $CFLAGS
    echo $CXXFLAGS
    echo $LDFLAGS

They point to files and folders from your ``sklearn-dev`` conda environment
(in particular in the bin/, include/ and lib/ subfolders). For instance
``-L/path/to/conda/envs/sklearn-dev/lib`` should appear in ``LDFLAGS``.

Notes on Conda
==============

Sometimes it can be necessary to open a new prompt before activating a newly
created conda environment.

If you get any conflicting dependency error messages on Mac or Linux, try commenting out
any custom conda configuration in the ``$HOME/.condarc`` file. In
particular the ``channel_priority: strict`` directive is known to cause
problems for this setup.

Notes on Meson
==============

When :ref:`building scikit-learn from source <install_from_source>`, existing
scikit-learn installations and meson builds can lead to conflicts.
You can use the `Makefile` provided in the `scikit-learn repository <https://github.com/scikit-learn/scikit-learn/>`__
to remove conflicting builds by calling:

.. prompt:: bash $

    make clean

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
