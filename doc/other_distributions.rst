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

The Debian/Ubuntu package is splitted in three difefrent packages called
``python3-sklearn`` (python modules), ``python3-sklearn-lib`` (low-level
implementations and bindings), ``python3-sklearn-doc`` (documentation).
Only the Python 3 version is available in the Debian Buster (the more recent one).
Packages can be installed using ``apt-get``::

    $ sudo apt-get install python3-sklearn python3-sklearn-lib python3-sklearn-doc


Fedora
------

The Fedora package is called ``python3-scikit-learn`` for the python 3 version,
the only one available in Fedora30.
It can be installed using ``dnf``::

    $ sudo dnf install python3-scikit-learn


NetBSD
------

scikit-learn is available via `pkgsrc-wip <http://pkgsrc-wip.sourceforge.net/>`_:

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

