
.. _advanced-installation:

==================================================
Installing the development version of scikit-learn
==================================================

This section introduces how to install the **master branch** of scikit-learn.
This is required:

- to test a feature that will be shipped in the next release (that is, a
  feature from a pull-request that was recently merged to the master branch);

- to check whether a bug you encountered has been fixed since the last
  release;

- or to work on a contribution of your own.

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

.. _git_repo:

Use `Git <https://git-scm.com/>`_ to check out the latest source from the
`scikit-learn repository <https://github.com/scikit-learn/scikit-learn>`_ on
`Github <https://github.com>`_ and pip to build the project::

    git clone git://github.com/scikit-learn/scikit-learn.git
    cd scikit-learn

    pip install cython pytest flake8
    pip install --editable .

Check that the installed scikit-learn with a version number ending with
`.dev0`::

    python -c "import sklearn; sklearn.show_versions()"

Run the tests with::

    pytest sklearn

Please read below for more details and :ref:`platform system specific
instructions <platform_specific_instructions>` (e.g. in case compilation error
message).

.. note::

    You will have to re-run the ``pip install --editable .`` command every
    time the source code of a compiled extension is changed (for instance when
    switching branches or pulling changes from upstream). Compiled extensions
    are Cython files (ending in `.pyx` or `.pxd`). See :ref:`editable_mode`
    for more details.

Dependencies
------------

Runtime dependencies
~~~~~~~~~~~~~~~~~~~~

Scikit-learn requires the following dependencies both at build time and at
runtime:

- Python (>= 3.5),
- NumPy (>= 1.11),
- SciPy (>= 0.17),
- Joblib (>= 0.11).

Those dependencies are **automatically installed by pip and conda** when
installing scikit-learn or building scikit-learn from source.

.. note::

   For running on PyPy, PyPy3-v5.10+, Numpy 1.14.0+, and scipy 1.1.0+
   are required. For PyPy, only installation instructions with pip apply.

Build dependencies
~~~~~~~~~~~~~~~~~~

Building Scikit-learn also requires:

- Cython >= 0.28.5
- A C/C++ compiler and a matching `OpenMP
  <https://en.wikipedia.org/wiki/OpenMP>`_ runtime library. See the
  :ref:`platform system specific instructions
  <platform_specific_instructions>` for more details.

.. note::

   It is possible to build scikit-learn without OpenMP support by setting the
   ``SKLEARN_NO_OPENMP`` environment variable (before cythonization). This is
   not recommended since it will force some estimators to run in sequential
   mode.

Test dependencies
~~~~~~~~~~~~~~~~~

Running tests requires:

.. |PytestMinVersion| replace:: 3.3.0

- pytest >=\ |PytestMinVersion|

Some tests also require `pandas <https://pandas.pydata.org>`_.


Building a specific version from a tag
--------------------------------------

If you want to build a stable version, you can ``git checkout <VERSION>``
to get the code for that particular version, or download an zip archive of
the version from github.

.. _editable_mode:

Editable mode
-------------

If you run the development version, it is cumbersome to reinstall the package
each time you update the sources. Therefore it is recommended that you install
in with the ``pip install --editable .`` command, which allows you to edit the
code in-place. This builds the extension in place and creates a link to the
development directory (see `the pip docs
<https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs>`_).

This is fundamentally similar to using the command ``python setup.py develop``
(see `the setuptool docs
<https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode>`_).
It is however preferred to use pip.

On Unix-like systems, you can equivalently type ``make in`` from the top-level
folder. Have a look at the ``Makefile`` for additional utilities.

.. _platform_specific_instructions:

Platform-specific instructions
==============================

Windows
-------

Using the official Microsoft compilers (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To build scikit-learn on Windows you need a working C/C++ compiler to build
the Cython extensions.

The building command depends on the architecture of the Python interpreter,
32-bit or 64-bit. You can check the architecture by running the following in
``cmd`` or ``powershell`` console::

    python -c "import struct; print(struct.calcsize('P') * 8)"

The above commands assume that you have the Python installation folder in your
PATH environment variable.

You will need `Build Tools for Visual Studio 2017
<https://visualstudio.microsoft.com/downloads/>`_.

.. warning::

    You DO NOT need to install Visual Studio 2019. You only need the "Build
    Tools for Visual Studio 2019", under "All downloads" -> "Tools for Visual
    Studio 2019".

For 64-bit Python, configure the build environment with::

    SET DISTUTILS_USE_SDK=1
    "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64

Please be aware that the path above might be different from user to user. The
aim is to point to the "vcvarsall.bat" file.

And build scikit-learn from this environment::

    pip install -e .

Replace ``x64`` by ``x86`` to build for 32-bit Python.

Using the MinGW compiler
~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to use `MinGW <http://www.mingw.org>`_ (a port of GCC to Windows
OS) as an alternative to MSVC for 32-bit Python. Not that extensions built with
mingw32 can be redistributed as reusable packages as they depend on GCC runtime
libraries typically not installed on end-users environment.

To force the use of a particular compiler, pass the ``--compiler`` flag to the
build step::

    python setup.py build --compiler=my_compiler develop

where ``my_compiler`` should be one of ``mingw32`` or ``msvc``.

macOS
-----

The default C compiler, Apple-clang, on Mac OSX does not directly support
OpenMP. We present two solutions to enable OpenMP support (you need to do only
one).

.. note::

    First, clean any previously built files in the source folder of
    scikit-learn::

        make clean

Using conda
~~~~~~~~~~~

One solution is to install another compiler which supports OpenMP. If you use
the conda package manager, you can install the ``compilers`` meta-package from
the conda-forge channel, which provides OpenMP-enabled C compilers.

It is recommended to use a dedicated conda environment to build scikit-learn
from source::

    conda create -n sklearn-dev python numpy scipy cython joblib pytest \
        conda-forge::compilers conda-forge::llvm-openmp
    conda activate sklearn-dev
    pip install --verbose --editable .

.. note::

    If you get any conflicting dependency error message, try commenting out
    any custom conda configuration in the ``$HOME/.condarc`` file. In
    particular the ``channel_priority: strict`` directive is known to cause
    problems for this setup.

You can check that the custom compilers are properly installed from conda
forge using the following command::

    conda list compilers llvm-openmp

The compilers meta-package will automatically set custom environment
variables::

    echo $CC
    echo $CXX
    echo $CFLAGS
    echo $CXXFLAGS
    echo $LDFLAGS

They point to files and folders from your sklearn-dev conda environment
(in particular in the bin/, include/ and lib/ subfolders).

The compiled extensions should be built with the clang and clang++ compilers
with the ``-fopenmp`` command line flag.

Using homebrew
~~~~~~~~~~~~~~

Another solution is to enable OpenMP support for the clang compiler shipped
by default on macOS.

You first need to install the OpenMP library::

    brew install libomp

Then you need to set the following environment variables::

    export CC=/usr/bin/clang
    export CXX=/usr/bin/clang++
    export CPPFLAGS="$CPPFLAGS -Xpreprocessor -fopenmp"
    export CFLAGS="$CFLAGS -I/usr/local/opt/libomp/include"
    export CXXFLAGS="$CXXFLAGS -I/usr/local/opt/libomp/include"
    export LDFLAGS="$LDFLAGS -Wl,-rpath,/usr/local/opt/libomp/lib -L/usr/local/opt/libomp/lib -lomp"

Finally, build scikit-learn in verbose mode::

    pip install --verbose --editable .

Linux
-----

Installing from source without conda requires you to have installed the
scikit-learn runtime dependencies, Python development headers and a working
C/C++ compiler. Under Debian-based operating systems, which include Ubuntu::

    sudo apt-get install build-essential python3-dev python3-setuptools \
                     python3-pip

and then::

    pip3 install numpy scipy cython

.. note::

    In order to build the documentation and run the example code contains in
    this documentation you will need matplotlib::

        pip3 install matplotlib

When precompiled wheels are not avalaible for your architecture, you can
install the system versions::

    sudo apt-get install cython3 python3-numpy python3-scipy python3-matplotlib

On Red Hat and clones (e.g. CentOS), install the dependencies using::

    sudo yum -y install gcc gcc-c++ python-devel numpy scipy

.. note::

    Scikit-learn will automatically detect the BLAS library (e.g. OpenBlas or
    MKL) used by SciPy at runtime. See `scipy installation instructions
    <https://docs.scipy.org/doc/scipy/reference/building/linux.html>`_ for
    more details.


FreeBSD
-------

The clang compiler included in FreeBSD 12.0 and 11.2 base systems does not
include OpenMP support. You need to install the `openmp` library from packages
(or ports)::

    sudo pkg install openmp

This will install header files in ``/usr/local/include`` and libs in
``/usr/local/lib``. Since these directories are not searched by default, you
can set the environment variables to these locations::

    export CFLAGS="$CFLAGS -I/usr/local/include"
    export CXXFLAGS="$CXXFLAGS -I/usr/local/include"
    export LDFLAGS="$LDFLAGS -Wl,-rpath,/usr/local/lib -L/usr/local/lib -lomp"

Finally you can build the package using the standard command.

For the upcomming FreeBSD 12.1 and 11.3 versions, OpenMP will be included in
the base system and these steps will not be necessary.
