
.. _advanced-installation:

===================================
Advanced installation instructions
===================================

There are different ways to get scikit-learn installed:

  * :ref:`Install an official release <install_official_release>`. This
    is the best approach for most users. It will provide a stable version
    and pre-build packages are available for most platforms.

  * Install the version of scikit-learn provided by your
    :ref:`operating system or Python distribution <install_by_distribution>`.
    This is a quick option for those who have operating systems
    that distribute scikit-learn. It might not provide the latest release
    version.

  * :ref:`Building the package from source
    <install_bleeding_edge>`. This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code. This document describes how to build from source.

.. note::

    If you wish to contribute to the project, you need to
    :ref:`install the latest development version<install_bleeding_edge>`.

.. _install_nightly_builds:

Installing nightly builds
=========================

The continuous integration servers of the scikit-learn project build, test
and upload wheel packages for the most recent Python version on a nightly
basis to help users test bleeding edge features or bug fixes::

  pip install --pre -f https://sklearn-nightly.scdn8.secure.raxcdn.com scikit-learn


.. _install_bleeding_edge:

Building from source
=====================

Scikit-learn requires:

- Python (>= 3.5),
- NumPy (>= 1.11),
- SciPy (>= 0.17).

.. note::

   For installing on PyPy, PyPy3-v5.10+, Numpy 1.14.0+, and scipy 1.1.0+
   are required. For PyPy, only installation instructions with pip apply.


Building Scikit-learn also requires

- Cython >=0.28.5
- OpenMP

Running tests requires

.. |PytestMinVersion| replace:: 3.3.0

- pytest >=\ |PytestMinVersion|

Some tests also require `pandas <https://pandas.pydata.org>`_.

.. _git_repo:

Retrieving the latest code
--------------------------

We use `Git <https://git-scm.com/>`_ for version control and
`GitHub <https://github.com/>`_ for hosting our main repository.

You can check out the latest sources with the command::

    git clone git://github.com/scikit-learn/scikit-learn.git

If you want to build a stable version, you can ``git checkout <VERSION>``
to get the code for that particular version, or download an zip archive of
the version from github.

If you have all the build requirements installed (see below for details), you
can build and install the package in the following way.

If you run the development version, it is cumbersome to reinstall the
package each time you update the sources. Therefore it's recommended that you
install in editable, which allows you to edit the code in-place. This
builds the extension in place and creates a link to the development directory
(see `the pip docs <https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs>`_)::

    pip install --editable .

.. note::

    This is fundamentally similar to using the command ``python setup.py develop``
    (see `the setuptool docs <https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode>`_).
    It is however preferred to use pip.

.. note::

    If you decide to do an editable install you have to rerun::

        pip install --editable .

    every time the source code of a compiled extension is
    changed (for instance when switching branches or pulling changes from upstream).

On Unix-like systems, you can simply type ``make`` in the top-level folder to
build in-place and launch all the tests. Have a look at the ``Makefile`` for
additional utilities.

Mac OSX
-------

The default C compiler, Apple-clang, on Mac OSX does not directly support
OpenMP. The first solution to build scikit-learn is to install another C
compiler such as gcc or llvm-clang. Another solution is to enable OpenMP
support on the default Apple-clang. In the following we present how to
configure this second option.

You first need to install the OpenMP library::

    brew install libomp

Then you need to set the following environment variables::

    export CC=/usr/bin/clang
    export CXX=/usr/bin/clang++
    export CPPFLAGS="$CPPFLAGS -Xpreprocessor -fopenmp"
    export CFLAGS="$CFLAGS -I/usr/local/opt/libomp/include"
    export CXXFLAGS="$CXXFLAGS -I/usr/local/opt/libomp/include"
    export LDFLAGS="$LDFLAGS -L/usr/local/opt/libomp/lib -lomp"
    export DYLD_LIBRARY_PATH=/usr/local/opt/libomp/lib

Finally you can build the package using the standard command.

Installing build dependencies
=============================

Linux
-----

Installing from source requires you to have installed the scikit-learn runtime
dependencies, Python development headers and a working C/C++ compiler.
Under Debian-based operating systems, which include Ubuntu::

    sudo apt-get install build-essential python3-dev python3-setuptools \
                     python3-numpy python3-scipy \
                     libatlas-dev libatlas3-base

On recent Debian and Ubuntu (e.g. Ubuntu 14.04 or later) make sure that ATLAS
is used to provide the implementation of the BLAS and LAPACK linear algebra
routines::

    sudo update-alternatives --set libblas.so.3 \
        /usr/lib/atlas-base/atlas/libblas.so.3
    sudo update-alternatives --set liblapack.so.3 \
        /usr/lib/atlas-base/atlas/liblapack.so.3

.. note::

    In order to build the documentation and run the example code contains in
    this documentation you will need matplotlib::

        sudo apt-get install python-matplotlib

.. note::

    The above installs the ATLAS implementation of BLAS
    (the Basic Linear Algebra Subprograms library).
    Ubuntu 11.10 and later, and recent (testing) versions of Debian,
    offer an alternative implementation called OpenBLAS.

    Using OpenBLAS can give speedups in some scikit-learn modules,
    but can freeze joblib/multiprocessing prior to OpenBLAS version 0.2.8-4,
    so using it is not recommended unless you know what you're doing.

    If you do want to use OpenBLAS, then replacing ATLAS only requires a couple
    of commands. ATLAS has to be removed, otherwise NumPy may not work::

        sudo apt-get remove libatlas3gf-base libatlas-dev
        sudo apt-get install libopenblas-dev

        sudo update-alternatives  --set libblas.so.3 \
            /usr/lib/openblas-base/libopenblas.so.0
        sudo update-alternatives --set liblapack.so.3 \
            /usr/lib/lapack/liblapack.so.3

On Red Hat and clones (e.g. CentOS), install the dependencies using::

    sudo yum -y install gcc gcc-c++ numpy python-devel scipy


Windows
-------

To build scikit-learn on Windows you need a working C/C++ compiler in
addition to numpy, scipy and setuptools.

The building command depends on the architecture of the Python interpreter,
32-bit or 64-bit. You can check the architecture by running the following in
``cmd`` or ``powershell`` console::

    python -c "import struct; print(struct.calcsize('P') * 8)"

The above commands assume that you have the Python installation folder in your
PATH environment variable.

You will need `Build Tools for Visual Studio 2017
<https://visualstudio.microsoft.com/de/downloads/>`_.

For 64-bit Python, configure the build environment with::

    SET DISTUTILS_USE_SDK=1
    "C:\Program Files (x86)\Microsoft Visual Studio\2017\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64

And build scikit-learn from this environment::

    python setup.py install

Replace ``x64`` by ``x86`` to build for 32-bit Python.


Building binary packages and installers
---------------------------------------

The ``.whl`` package and ``.exe`` installers can be built with::

    pip install wheel
    python setup.py bdist_wheel bdist_wininst -b doc/logos/scikit-learn-logo.bmp

The resulting packages are generated in the ``dist/`` folder.


Using an alternative compiler
-----------------------------

It is possible to use `MinGW <http://www.mingw.org>`_ (a port of GCC to Windows
OS) as an alternative to MSVC for 32-bit Python. Not that extensions built with
mingw32 can be redistributed as reusable packages as they depend on GCC runtime
libraries typically not installed on end-users environment.

To force the use of a particular compiler, pass the ``--compiler`` flag to the
build step::

    python setup.py build --compiler=my_compiler install

where ``my_compiler`` should be one of ``mingw32`` or ``msvc``.


.. _testing:

Testing
=======

Testing scikit-learn once installed
-----------------------------------

Testing requires having `pytest <https://docs.pytest.org>`_ >=\ |PytestMinVersion|\ .
Some tests also require having `pandas <https://pandas.pydata.org/>` installed.
After installation, the package can be tested by executing *from outside* the
source directory::

    $ pytest sklearn

This should give you a lot of output (and some warnings) but
eventually should finish with a message similar to::

    =========== 8304 passed, 26 skipped, 4659 warnings in 557.76 seconds ===========

Otherwise, please consider posting an issue into the `GitHub issue tracker
<https://github.com/scikit-learn/scikit-learn/issues>`_ or to the
:ref:`mailing_lists` including the traceback of the individual failures
and errors. Please include your operating system, your version of NumPy, SciPy
and scikit-learn, and how you installed scikit-learn.


Testing scikit-learn from within the source folder
--------------------------------------------------

Scikit-learn can also be tested without having the package
installed. For this you must compile the sources inplace from the
source directory::

    python setup.py build_ext --inplace

Test can now be run using pytest::

    pytest sklearn

This is automated by the commands::

    make in

and::

    make test


You can also install a symlink named ``site-packages/scikit-learn.egg-link``
to the development folder of scikit-learn with::

    pip install --editable .
