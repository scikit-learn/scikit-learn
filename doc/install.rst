=========================
Installing `scikit-learn`
=========================

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
    <install_bleeding_edge>`.  This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code.

.. note::

    If you wish to contribute to the project, it's recommended you
    :ref:`install the latest development version<install_bleeding_edge>`.


.. _install_official_release:

Installing an official release
==============================


Getting the dependencies
------------------------

Installing from source requires you to have installed Python (>= 2.6),
NumPy (>= 1.3), SciPy (>= 0.7), setuptools, Python development headers
and a working C++ compiler.
Under Debian-based operating systems, which include Ubuntu,
you can install all these requirements by issuing::

    sudo apt-get install build-essential python-dev python-numpy python-setuptools python-scipy libatlas-dev libatlas3-base

.. note::

    In order to build the documentation and run the example code contains in
    this documentation you will need matplotlib::

        sudo apt-get install python-matplotlib

.. note::

    On older versions of Ubuntu,
    you might need to ``apt-get install python-numpy-dev``
    to get the header files for NumPy.

    On Ubuntu 10.04 LTS, the package `libatlas-dev` is called `libatlas-headers`.

.. note::

    The above installs the ATLAS implementation of BLAS
    (the Basic Linear Algebra Subprograms library).
    Ubuntu 11.10 and later, and recent (testing) versions of Debian,
    offer an alternative implementation called OpenBLAS.

    Using OpenBLAS can give speedups in some scikit-learn modules,
    but it doesn't play nicely with joblib/multiprocessing,
    so using it is not recommended unless you know what you're doing.

    If you do want to use OpenBLAS, then replacing ATLAS
    only requires two commands.
    ATLAS has to be removed, otherwise NumPy may not work::

        sudo apt-get remove libatlas3gf-base libatlas-dev
        sudo apt-get install libopenblas-dev

On Red Hat and clones (e.g. CentOS), install the dependencies using::

    sudo yum -y install gcc gcc-c++ numpy python-devel scipy

Easy install
~~~~~~~~~~~~

This is usually the fastest way to install the latest stable
release. If you have pip or easy_install, you can install or update
with the command::

    pip install -U scikit-learn

or::

    easy_install -U scikit-learn

for easy_install. Note that you might need root privileges to run
these commands.


From source package
~~~~~~~~~~~~~~~~~~~

Download the package from http://pypi.python.org/pypi/scikit-learn/
, unpack the sources and cd into the source directory.

This packages uses distutils, which is the default way of installing
python modules. The install command is::

  python setup.py install


Windows installer
-----------------

You can download a Windows installer from `downloads
<https://sourceforge.net/projects/scikit-learn/files/>`_ in the
project's web page. Note that must also have installed the packages
numpy and setuptools.

This package is also expected to work with python(x,y) as of 2.6.5.5.

.. topic:: **Installing on Windows 64-bit**

   To install a 64-bit version of scikit-learn, you can download the
   binaries from http://www.lfd.uci.edu/~gohlke/pythonlibs/#scikit-learn
   Note that this will require a compatible version of numpy, scipy and
   matplotlib. The easiest option is to also download them from the same
   URL.

Building on windows
-------------------

To build scikit-learn on windows you will need a C/C++ compiler in
addition to numpy, scipy and setuptools. At least
`MinGW <http://www.mingw.org>`_ (a port of GCC to Windows OS) and
Microsoft Visual C++ 2008 should work out of the box. To force the use
of a particular compiler, write a file named ``setup.cfg`` in the
source directory with the content::

    [build_ext]
    compiler=my_compiler

    [build]
    compiler=my_compiler

where ``my_compiler`` should be one of ``mingw32`` or ``msvc``.

When the appropriate compiler has been set, and assuming Python is
in your PATH (see
`Python FAQ for windows <http://docs.python.org/faq/windows.html>`_
for more details), installation is done by
executing the command::

    python setup.py install


To build a precompiled package like the ones distributed at
`the downloads section <https://sourceforge.net/projects/scikit-learn/files/>`_,
the command to execute is::

    python setup.py bdist_wininst -b doc/logos/scikit-learn-logo.bmp

This will create an installable binary under directory ``dist/``.


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
<http://neuro.debian.net/pkgs/python-sklearn.html>`__ .

A quick-'n'-dirty way of rolling your own ``.deb`` package
is to `use stdeb <https://github.com/scikit-learn/scikit-learn/wiki/Quick-packaging-for-Debian-Ubuntu>`_.

Python(x,y)
-----------

The `Python(x,y) <http://pythonxy.com>`_ distributes scikit-learn as an additional plugin, which can
be found in the `Additional plugins <http://code.google.com/p/pythonxy/wiki/AdditionalPlugins>`_
page.


Enthought Python distribution
-----------------------------

The `Enthought Python Distribution
<http://www.enthought.com/products/epd.php>`_ already ships a recent
version.


MacPorts
--------

The MacPorts package is named ``py<XY>-scikits-learn``,
where ``XY`` denotes the Python version.
It can be installed by typing the following
command::

    sudo port install py26-scikit-learn

or::

    sudo port install py27-scikit-learn


Archlinux
---------

Archlinux's package is provided at 
`Arch User Repository (AUR) <https://aur.archlinux.org/>`_ with name
`python2-scikit-learn` for latest stable version and `python2-scikit-learn-git`
for building from git version. If `yaourt` is available, it can be installed
by typing the following command::

     sudo yaourt -S python2-scikit-learn

or::

     sudo yaourt -S python2-scikit-learn-git

depending on the version of scikit-learn you want to use.


NetBSD
------

scikit-learn is available via `pkgsrc-wip <http://pkgsrc-wip.sourceforge.net/>`_:

    http://pkgsrc.se/wip/py-scikit_learn

.. _install_bleeding_edge:

Bleeding Edge
=============

See section :ref:`git_repo` on how to get the development version.


.. _testing:

Testing
=======

Testing requires having the `nose
<http://somethingaboutorange.com/mrl/projects/nose/>`_ library. After
installation, the package can be tested by executing *from outside* the
source directory::

    nosetests sklearn --exe

This should give you a lot of output (and some warnings) but
eventually should finish with a message similar to::

           Ran 601 tests in 27.920s
           OK (SKIP=2)

Otherwise, please consider posting an issue into the `bug tracker
<https://github.com/scikit-learn/scikit-learn/issues>`_ or to the
:ref:`mailing_lists`.

.. note:: **Alternative testing method**

   If for some reason the recommended method is failing for you, please try
   the alternate method::

    python -c "import sklearn; sklearn.test()"

   This method might display doctest failures because of nosetests issues.

scikit-learn can also be tested without having the package
installed. For this you must compile the sources inplace from the
source directory::

    python setup.py build_ext --inplace

Test can now be run using nosetests::

    nosetests sklearn/

This is automated by the commands::

    make in

and::

    make test
