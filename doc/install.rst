.. _installation-instructions:

=======================
Installing scikit-learn
=======================

There are different ways to get scikit-learn installed:

  * Install the version of scikit-learn provided by your
    :ref:`operating system or Python distribution <install_by_distribution>`.
    This is the quickest option for those who have operating systems that
    distribute scikit-learn.

  * :ref:`Install an official release <install_official_release>`. This
    is the best approach for users who want a stable version number
    and aren't concerned about running a slightly older version of
    scikit-learn.

  * :ref:`Install the latest development version
    <install_bleeding_edge>`. This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code.

.. note::

    If you wish to contribute to the project, it's recommended you
    :ref:`install the latest development version<install_bleeding_edge>`.


.. _install_official_release:

Installing an official release
==============================

Scikit-learn requires:

- Python (>= 2.6 or >= 3.3),
- NumPy (>= 1.6.1),
- SciPy (>= 0.9).


Windows
-------

First you need to install `numpy <http://numpy.scipy.org/>`_ and `scipy
<http://www.scipy.org/>`_ from their own official installers.

Wheel packages (.whl files) for scikit-learn from `PyPI
<https://pypi.python.org/pypi/scikit-learn/>`_ can be installed with the `pip
<http://pip.readthedocs.org/en/latest/installing.html>`_ utility.
Open a console and type the following to install or upgrade scikit-learn to the
latest stable release::

    pip install -U scikit-learn

If there are no binary packages matching your Python version you might
to try to install scikit-learn and its dependencies from `Christoph Gohlke
Unofficial Windows installers
<http://www.lfd.uci.edu/~gohlke/pythonlibs/#scikit-learn>`_
or from a :ref:`Python distribution <install_by_distribution>` instead.


Mac OSX
-------

Scikit-learn and its dependencies are all available as wheel packages for OSX::

    pip install -U numpy scipy scikit-learn


Linux
-----

At this time scikit-learn does not provide official binary packages for Linux
so you have to build from source.


Installing build dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installing from source requires you to have installed the scikit-learn runtime
dependencies, Python development headers and a working C/C++ compiler.
Under Debian-based operating systems, which include Ubuntu, if you have
Python 2 you can install all these requirements by issuing::

    sudo apt-get install build-essential python-dev python-setuptools \
                         python-numpy python-scipy \
                         libatlas-dev libatlas3gf-base

If you have Python 3::

    sudo apt-get install build-essential python3-dev python3-setuptools \
                         python3-numpy python3-scipy \
                         libatlas-dev libatlas3gf-base

On recent Debian and Ubuntu (e.g. Ubuntu 13.04 or later) make sure that ATLAS
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


Building scikit-learn with pip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is usually the fastest way to install or upgrade to the latest stable
release::

    pip install --user --install-option="--prefix=" -U scikit-learn

The ``--user`` flag ask pip to install scikit-learn in the ``$HOME/.local``
folder therefore not requiring root permission. This flag should make pip
ignore any old version of scikit-learn previously installed on the system while
benefitting from system packages for numpy and scipy. Those dependencies can
be long and complex to build correctly from source.

The ``--install-option="--prefix="`` flag is only required if Python has a
``distutils.cfg`` configuration with a predefined ``prefix=`` entry.


From source package
~~~~~~~~~~~~~~~~~~~

Download the source package from http://pypi.python.org/pypi/scikit-learn/
, unpack the sources and cd into the source directory.

This packages uses distutils, which is the default way of installing
python modules. The install command is::

    python setup.py install


.. _install_by_distribution:

Third party distributions of scikit-learn
=========================================

Some third-party distributions are now providing versions of
scikit-learn integrated with their package-management systems.

These can make installation and upgrading much easier for users since
the integration includes the ability to automatically install
dependencies (numpy, scipy) that scikit-learn requires.

The following is an incomplete list of Python and OS distributions
that provide their own version of scikit-learn.


Debian and derivatives (Ubuntu)
-------------------------------

The Debian package is named ``python-sklearn``
(formerly ``python-scikits-learn``)
and can be installed using the following command::

      sudo apt-get install python-sklearn

Additionally, backport builds of the most recent release of
scikit-learn for existing releases of Debian and Ubuntu are available
from the `NeuroDebian repository
<http://neuro.debian.net/pkgs/python-sklearn.html>`_ .

A quick-'n'-dirty way of rolling your own ``.deb`` package
is to `use stdeb <https://github.com/scikit-learn/scikit-learn/wiki/Quick-packaging-for-Debian-Ubuntu>`_.


Python(x,y) for Windows
-----------------------

The `Python(x,y) <https://code.google.com/p/pythonxy/>`_ project distributes
scikit-learn as an additional plugin, which can be found in the `Additional
plugins <http://code.google.com/p/pythonxy/wiki/AdditionalPlugins>`_ page.


Canopy and Anaconda for all supported platforms
-----------------------------------------------

`Canopy
<http://www.enthought.com/products/canopy>`_ and `Anaconda
<https://store.continuum.io/cshop/anaconda/>`_ ships a recent
version, in addition to a large set of scientific python library.


MacPorts for Mac OSX
--------------------

The MacPorts package is named ``py<XY>-scikits-learn``,
where ``XY`` denotes the Python version.
It can be installed by typing the following
command::

    sudo port install py26-scikit-learn

or::

    sudo port install py27-scikit-learn


Arch Linux
----------

Arch Linux's package is provided through the `official repositories
<https://www.archlinux.org/packages/?q=scikit-learn>`_ as
``python-scikit-learn`` for Python 3 and ``python2-scikit-learn`` for Python 2.
It can be installed by typing the following command:

.. code-block:: none

     # pacman -S python-scikit-learn

or:

.. code-block:: none

     # pacman -S python2-scikit-learn

depending on the version of Python you use.


NetBSD
------

scikit-learn is available via `pkgsrc-wip <http://pkgsrc-wip.sourceforge.net/>`_:

    http://pkgsrc.se/wip/py-scikit_learn

Fedora
------

The Fedora package is called ``python-scikit-learn`` for the Python 2 version
and ``python3-scikit-learn`` for the Python 3 version. Both versions can
be installed using ``yum``::

    $ sudo yum install python-scikit-learn

or::

    $ sudo yum install python3-scikit-learn


Building on windows
===================

To build scikit-learn on Windows you need a working C/C++ compiler in
addition to numpy, scipy and setuptools.

Picking the right compiler depends on the version of Python (2 or 3)
and the architecture of the Python interpreter, 32-bit or 64-bit.
You can check the Python version by running the following in ``cmd`` or
``powershell`` console::

    python --version

and the architecture with::

    python -c "import struct; print(struct.calcsize('P') * 8)"

The above commands assume that you have the Python installation folder in your
PATH environment variable.


32-bit Python
-------------

For 32-bit Python it is possible use the standalone installers for
`Microsoft Visual C++ Express 2008 <http://go.microsoft.com/?linkid=7729279>`_
for Python 2 or
`Microsoft Visual C++ Express 2010 <http://go.microsoft.com/?linkid=9709949>`_
or Python 3.

Once installed you should be able to build scikit-learn without any
particular configuration by running the following command in the scikit-learn
folder::

   python setup.py install


64-bit Python
-------------

For the 64-bit architecture, you either need the full Visual Studio or
the free Windows SDKs that can be downloaded from the links below.

The Windows SDKs include the MSVC compilers both for 32 and 64-bit
architectures. They come as a ``GRMSDKX_EN_DVD.iso`` file that can be mounted
as a new drive with a ``setup.exe`` installer in it.

- For Python 2 you need SDK **v7.0**: `MS Windows SDK for Windows 7 and .NET
  Framework 3.5 SP1
  <http://www.microsoft.com/en-us/download/details.aspx?id=18950>`_

- For Python 3 you need SDK **v7.1**: `MS Windows SDK for Windows 7 and .NET
  Framework 4
  <https://www.microsoft.com/en-us/download/details.aspx?id=8442>`_

Both SDKs can be installed in parallel on the same host. To use the Windows
SDKs, you need to setup the environment of a ``cmd`` console launched with the
following flags (at least for SDK v7.0)::

    cmd /E:ON /V:ON /K

Then configure the build environment with::

    SET DISTUTILS_USE_SDK=1
    SET MSSdk=1
    "C:\Program Files\Microsoft SDKs\Windows\v7.0\Setup\WindowsSdkVer.exe" -q -version:v7.0
    "C:\Program Files\Microsoft SDKs\Windows\v7.0\Bin\SetEnv.cmd" /x64 /release

Finally you can build scikit-learn in the same ``cmd`` console::

    python setup.py install

Replace ``v7.0`` by the ``v7.1`` in the above commands to do the same for
Python 3 instead of Python 2.

Replace ``/x64`` by ``/x86``  to build for 32-bit Python instead of 64-bit
Python.


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


.. _install_bleeding_edge:

Bleeding Edge
=============

See section :ref:`git_repo` on how to get the development version. Then follow
the previous instructions to build from source depending on your platform.


.. _testing:

Testing
=======

Testing scikit-learn once installed
-----------------------------------

Testing requires having the `nose
<http://somethingaboutorange.com/mrl/projects/nose/>`_ library. After
installation, the package can be tested by executing *from outside* the
source directory::

    $ nosetests -v sklearn

Under Windows, it is recommended to use the following command (adjust the path
to the ``python.exe`` program) as using the ``nosetests.exe`` program can badly
interact with tests that use ``multiprocessing``::

    C:\Python34\python.exe -c "import nose; nose.main()" -v sklearn

This should give you a lot of output (and some warnings) but
eventually should finish with a message similar to::

    Ran 3246 tests in 260.618s
    OK (SKIP=20)

Otherwise, please consider posting an issue into the `bug tracker
<https://github.com/scikit-learn/scikit-learn/issues>`_ or to the
:ref:`mailing_lists` including the traceback of the individual failures
and errors.


Testing scikit-learn from within the source folder
--------------------------------------------------

Scikit-learn can also be tested without having the package
installed. For this you must compile the sources inplace from the
source directory::

    python setup.py build_ext --inplace

Test can now be run using nosetests::

    nosetests -v sklearn/

This is automated by the commands::

    make in

and::

    make test


You can also install a symlink named ``site-packages/scikit-learn.egg-link``
to the development folder of scikit-learn with::

    pip install --editable .
