===============================
Installing `scikits.learn`
===============================

There are different ways to get scikits.learn installed:

  * Install the version of scikits.learn provided by your
    :ref:`operating system distribution <install_by_distribution>` . This
    is the quickest option for those who have operating systems that
    distribute scikits.learn.

  * :ref:`Install an official release <install_official_release>`. This
    is the best approach for users who want a stable version number
    and aren't concerned about running a slightly older version of
    scikits.learn.

  * :ref:`Install the latest development version
    <install_bleeding_edge>`.  This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code.



.. _install_official_release:

Installing an official release
==============================


Installing from source
----------------------

Installing from source requires you to have installed numpy, 
scipy, setuptools, python development headers and a working C++
compiler. Under debian-like systems you can get all this by executing
with root privileges::

    sudo apt-get install python-dev python-numpy python-numpy-dev python-setuptools python-numpy-dev python-scipy libatlas-dev g++

.. note::

    In Order to build the documentation and run the example code contains in
    this documentation you will need matplotlib::

        sudo apt-get install python-matplotlib

.. note::

    On Ubuntu LTS (10.04) the package `libatlas-dev` is called `libatlas-headers`

Easy install
~~~~~~~~~~~~

This is usually the fastest way to install the latest stable
release. If you have pip or easy_install, you can install or update
with the command::

    pip install -U scikits.learn

or::

    easy_install -U scikits.learn

for easy_install. Note that you might need root privileges to run
these commands.


From source package
~~~~~~~~~~~~~~~~~~~

Download the package from http://sourceforge.net/projects/scikit-learn/files
, unpack the sources and cd into archive.

This packages uses distutils, which is the default way of installing
python modules. The install command is::

  python setup.py install


Windows installer
-----------------

You can download a windows installer from `downloads
<https://sourceforge.net/projects/scikit-learn/files/>`_ in the
project's web page. Note that must also have installed the packages
numpy and setuptools.

This package is also expected to work with python(x,y) as of 2.6.5.5.


.. _build_on_windows

Building on windows
-------------------

To build scikits.learn on windows you will need a C/C++ compiler in
addition to numpy, scipy and setuptools. At least
`MinGW <http://www.mingw.org>`_ (a port of GCC to Windows OS) and the
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

Third party distributions of scikits.learn
==========================================

Some third-party distributions are now providing versions of
scikits.learn integrated with their package-management systems. 

These can make installation and upgrading much easier for users since
the integration includes the ability to automatically install
dependencies (numpy, scipy) that scikits.learn requires.

The following is a list of linux distributions that provide their own
version of scikits.learn:


Debian and derivatives (Ubuntu)
-------------------------------

The Debian package is named python-scikits-learn and can be install
using the following commands with root privileges::

      apt-get install python-scikits-learn


Enthought python distribution
-----------------------------

The `Enthought Python Distribution
<http://www.enthought.com/products/epd.php>`_ already ships the latest
version.


Macports
--------

The macport's package is named py26-scikits-learn and can be installed
by typing the following command::

    sudo port install py26-scikits-learn

NetBSD
------

scikits.learn is available via `pkgsrc-wip <http://pkgsrc-wip.sourceforge.net/>`_:

    http://pkgsrc.se/wip/py-scikits_learn

.. _install_bleeding_edge:

Bleeding Edge
=============

See section :ref:`git_repo` on how to get the development version.


.. _testing:

Testing
=======

Testing requires having the `nose
<http://somethingaboutorange.com/mrl/projects/nose/>`_ library. After
installation, the package can be tested by executing from outside the
source directory::

    python -c "import scikits.learn as skl; skl.test()"

This should give you a lot of output (and some warnings) but
eventually should finish with the a text similar to::

           Ran 601 tests in 27.920s
           OK (SKIP=2)

otherwise please consider submitting a bug in the :ref:`bug_tracker`
or to the :ref:`mailing_lists`.

scikits.learn can also be tested without having the package
installed. For this you must compile the sources inplace from the
source directory::

    python setup.py build_ext --inplace

Test can now be run using nosetest::

     nosetests scikits/learn/

If you are running the deveopment version, this is automated in the
commands `make in` and `make test`.

.. warning::

   Because nosetest does not play well with multiprocessing on
   windows, this last approach is not recommended on such system.
