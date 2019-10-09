
.. _advanced-installation:

==================================================================
Installing the development version of scikit-learn (master branch)
==================================================================

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

In the vast majority of cases, building scikit-learn for development purposes
can be done with::

    pip install cython pytest flake8

Then, in the main repository::

    pip install --editable .

Please read below for details and more advanced instructions.

Dependencies
------------

Scikit-learn requires:

- Python (>= 3.5),
- NumPy (>= 1.11),
- SciPy (>= 0.17),
- Joblib (>= 0.11).

.. note::

   For installing on PyPy, PyPy3-v5.10+, Numpy 1.14.0+, and scipy 1.1.0+
   are required. For PyPy, only installation instructions with pip apply.


Building Scikit-learn also requires

- Cython >=0.28.5
- OpenMP

.. note::

   It is possible to build scikit-learn without OpenMP support by setting the
   ``SKLEARN_NO_OPENMP`` environment variable (before cythonization). This is
   not recommended since it will force some estimators to run in sequential
   mode and their ``n_jobs`` parameter will be ignored.


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

Once you have all the build requirements installed (see below for details),
you can build and install the package in the following way.

If you run the development version, it is cumbersome to reinstall the
package each time you update the sources. Therefore it's recommended that you
install in editable mode, which allows you to edit the code in-place. This
builds the extension in place and creates a link to the development directory
(see `the pip docs <https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs>`_)::

    pip install --editable .

.. note::

    This is fundamentally similar to using the command ``python setup.py develop``
    (see `the setuptool docs <https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode>`_).
    It is however preferred to use pip.

.. note::

    You will have to re-run::

        pip install --editable .

    every time the source code of a compiled extension is changed (for
    instance when switching branches or pulling changes from upstream).
    Compiled extensions are Cython files (ending in `.pyx` or `.pxd`).

On Unix-like systems, you can equivalently type ``make in`` from the
top-level folder. Have a look at the ``Makefile`` for additional utilities.

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
    export LDFLAGS="$LDFLAGS -L/usr/local/lib -lomp"
    export DYLD_LIBRARY_PATH=/usr/local/lib

Finally you can build the package using the standard command.

For the upcomming FreeBSD 12.1 and 11.3 versions, OpenMP will be included in 
the base system and these steps will not be necessary.


Installing build dependencies
=============================

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

    To use a high performance BLAS library (e.g. OpenBlas) see 
    `scipy installation instructions
    <https://docs.scipy.org/doc/scipy/reference/building/linux.html>`_.

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
<https://visualstudio.microsoft.com/downloads/>`_.

.. warning::
	You DO NOT need to install Visual Studio 2019. 
	You only need the "Build Tools for Visual Studio 2019", 
	under "All downloads" -> "Tools for Visual Studio 2019". 

For 64-bit Python, configure the build environment with::

    SET DISTUTILS_USE_SDK=1
    "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" x64

Please be aware that the path above might be different from user to user. 
The aim is to point to the "vcvarsall.bat" file.

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
