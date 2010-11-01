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


Easy install
------------

This is usually the fastest way to install the latest stable
release. If you have pip or easy_install, you can install or update
with the command::

    pip install -U scikits.learn

or::

    easy_install -U scikits.learn

for easy_install.


Windows installer
-----------------

You can download a windows installer from `downloads
<https://sourceforge.net/projects/scikit-learn/files/>`_
in the project's web page.


From Source
-----------
Download the package from http://sourceforge.net/projects/scikit-learn/files
, unpack the sources and cd into archive.

This packages uses distutils, which is the default way of installing
python modules. The install command is::

  python setup.py install


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



.. _install_bleeding_edge:

Bleeding Edge
=============

See section :ref:`git_repo` on how to get the development version.
