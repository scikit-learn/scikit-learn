
.. _misc-info:

==================================================
Miscellaneous information / Troubleshooting
==================================================

Here, you find some more advanced notes and troubleshooting tips related to
:ref:`setup_development_environment`.

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
(in particular in the `bin/`, `include/` and `lib/` subfolders). For instance
``-L/path/to/conda/envs/sklearn-dev/lib`` should appear in ``LDFLAGS``.

Notes on Conda
==============

Sometimes it can be necessary to open a new prompt before activating a newly
created conda environment.

If you get any conflicting dependency error messages on Mac or Linux, try commenting out
any custom conda configuration in the ``$HOME/.condarc`` file. In
particular the ``channel_priority: strict`` directive is known to cause
problems for this setup.

Note on dependencies for other Linux distributions
==================================================

When precompiled wheels of the runtime dependencies are not available for your
architecture (e.g. **ARM**), you can install the system versions:

.. prompt::

  sudo apt-get install cython3 python3-numpy python3-scipy


Notes on Meson
==============

When :ref:`building scikit-learn from source <install_from_source>`, existing
scikit-learn installations and meson builds can lead to conflicts.
You can use the `Makefile` provided in the `scikit-learn repository <https://github.com/scikit-learn/scikit-learn/>`__
to remove conflicting builds by calling:

.. prompt:: bash $

    make clean
