.. _installation-instructions:

=======================
Installing scikit-learn
=======================

There are different ways to get scikit-learn installed:

  * :ref:`Install the latest official release <install_official_release>`. This
    is the best approach for most users. It will provide a stable version
    and pre-build packages are available for most platforms.
    Note that :ref:`nightly builds <install_nightly_builds>` are also
    distributed.

  * Install the version of scikit-learn provided by your
    :ref:`operating system or Python distribution <install_by_distribution>`.
    This is a quick option for those who have operating systems or Python
    distributions that distribute scikit-learn.
    It might not provide the latest release version.

  * :ref:`Building the package from source
    <install_bleeding_edge>`. This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code.

.. note ::

   If you wish to contribute to the project, it's recommended you
   :ref:`install the latest development version<install_bleeding_edge>`.

.. _install_official_release:

Installing the latest release
=============================

Scikit-learn requires:

- Python (>= 3.5)
- NumPy (>= 1.11.0)
- SciPy (>= 0.17.0)
- joblib (>= 0.11)

Scikit-learn plotting capabilities (i.e., functions start with "plot\_"
and classes end with "Display") require Matplotlib (>= 1.5.1). For running the
examples Matplotlib >= 1.5.1 is required. A few examples require
scikit-image >= 0.12.3, a few examples require pandas >= 0.18.0.

.. warning::

    Scikit-learn 0.20 was the last version to support Python 2.7 and Python 3.4.
    Scikit-learn now requires Python 3.5 or newer.

If you already have a working installation of numpy and scipy,
the easiest way to install scikit-learn is using ``conda`` (see the
`instructions for downloading conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html>`_)::

    conda install scikit-learn

or ``pip``.
In that case, in order to avoid any OS dependent issue it is strongly
recommended to use python3 ``virtualenv``
(see `python3 virtualenv documentation
<https://docs.python.org/3/tutorial/venv.html>`_).

That gives, on Windows::

    python -m venv .myenv
    .myenv\Scripts\activate
    pip install -U scikit-learn

or on Linux::

    python3 -m venv .myenv
    source .myenv/bin/activate
    pip install -U scikit-learn


If you have not installed NumPy or SciPy yet, you can also install these using
conda or pip. When using pip, please ensure that *binary wheels* are used,
and NumPy and SciPy are not recompiled from source, which can happen when using
particular configurations of operating system and hardware (such as Linux on
a Raspberry Pi). 

If you must install scikit-learn and its dependencies with pip, you can install
it as ``scikit-learn[alldeps]``.

.. note::

   For installing on PyPy, PyPy3-v5.10+, Numpy 1.14.0+, and scipy 1.1.0+
   are required.

.. _install_nightly_builds:

Installing nightly builds
=========================

The continuous integration servers of the scikit-learn project build, test
and upload wheel packages for the most recent Python version on a nightly
basis to help users test bleeding edge features or bug fixes::

  pip install --pre -f https://sklearn-nightly.scdn8.secure.raxcdn.com scikit-learn

Again, in order to avoid any OS dependent issue it is strongly
recommended to use python3 ``virtualenv``
(see `python3 virtualenv documentation
<https://docs.python.org/3/tutorial/venv.html>`_).

.. _install_by_distribution:

Third party distributions of scikit-learn
=========================================

Some third-party distributions provide versions of
scikit-learn integrated with their package-management systems.

These can make installation and upgrading much easier for users since
the integration includes the ability to automatically install
dependencies (numpy, scipy) that scikit-learn requires.

The following is an incomplete list of OS and python distributions
that provide their own version of scikit-learn.

Arch Linux
----------

Arch Linux's package is provided through the `official repositories
<https://www.archlinux.org/packages/?q=scikit-learn>`_ as
``python-scikit-learn`` for Python.
It can be installed by typing the following command:

.. code-block:: none

     # pacman -S python-scikit-learn


Debian/Ubuntu
-------------

The Debian/Ubuntu package is splitted in three different packages called
``python3-sklearn`` (python modules), ``python3-sklearn-lib`` (low-level
implementations and bindings), ``python3-sklearn-doc`` (documentation).
Only the Python 3 version is available in the Debian Buster (the more recent
Debian distribution).
Packages can be installed using ``apt-get``::

    $ sudo apt-get install python3-sklearn python3-sklearn-lib
          python3-sklearn-doc


Fedora
------

The Fedora package is called ``python3-scikit-learn`` for the python 3 version,
the only one available in Fedora30.
It can be installed using ``dnf``::

    $ sudo dnf install python3-scikit-learn


NetBSD
------

scikit-learn is available via `pkgsrc-wip
<http://pkgsrc-wip.sourceforge.net/>`_:

    http://pkgsrc.se/wip/py-scikit_learn


MacPorts for Mac OSX
--------------------

The MacPorts package is named ``py<XY>-scikits-learn``,
where ``XY`` denotes the Python version.
It can be installed by typing the following
command::

    sudo port install py27-scikit-learn

or::

    sudo port install py36-scikit-learn


Canopy and Anaconda for all supported platforms
-----------------------------------------------

`Canopy
<https://www.enthought.com/products/canopy>`_ and `Anaconda
<https://www.anaconda.com/download>`_ both ship a recent
version of scikit-learn, in addition to a large set of scientific python
library for Windows, Mac OSX and Linux.

Anaconda offers scikit-learn as part of its free distribution.


WinPython for Windows
-----------------------

The `WinPython <https://winpython.github.io/>`_ project distributes
scikit-learn as an additional plugin.
