.. _installation-instructions:

=======================
Installing scikit-learn
=======================

.. note::

    If you wish to contribute to the project, it's recommended you
    :ref:`install the latest development version<install_bleeding_edge>`.


Installing the latest release
=============================

Scikit-learn requires:

- Python (>= 2.6 or >= 3.3),
- NumPy (>= 1.6.1),
- SciPy (>= 0.9).

If you already have a working installation of numpy and scipy,
the easiest way to install scikit-learn is using ``pip`` ::

    pip install -U scikit-learn

or ``conda``::

    conda install scikit-learn

If you have not installed NumPy or SciPy yet, you can also install these using
conda or pip. When using pip, please ensure that *binary wheels* are used,
and NumPy and SciPy are not recompiled from source, which can happen when using
particular configurations of operating system and hardware (such as Linux on
a Raspberry Pi). 
Building numpy and scipy from source can be complex (especially on Windows) and
requires careful configuration to ensure that they link against an optimized
implementation of linear algebra routines.
Instead, use a third-party distribution as described below.

If you must install scikit-learn and its dependencies with pip, you can install
it as ``scikit-learn[alldeps]``. The most common use case for this is in a
``requirements.txt`` file used as part of an automated build process for a PaaS
application or a Docker image. This option is not intended for manual
installation from the command line.

Third-party Distributions
==========================
If you don't already have a python installation with numpy and scipy, we
recommend to install either via your package manager or via a python bundle.
These come with numpy, scipy, scikit-learn, matplotlib and many other helpful
scientific and data processing libraries.

Available options are:

Canopy and Anaconda for all supported platforms
-----------------------------------------------

`Canopy
<https://www.enthought.com/products/canopy>`_ and `Anaconda
<https://www.continuum.io/downloads>`_ both ship a recent
version of scikit-learn, in addition to a large set of scientific python
library for Windows, Mac OSX and Linux.

Anaconda offers scikit-learn as part of its free distribution.


.. warning::

    To upgrade or uninstall scikit-learn installed with Anaconda
    or ``conda`` you **should not use the pip command**. Instead:

    To upgrade ``scikit-learn``::

        conda update scikit-learn

    To uninstall ``scikit-learn``::

        conda remove scikit-learn

    Upgrading with ``pip install -U scikit-learn`` or uninstalling
    ``pip uninstall scikit-learn`` is likely fail to properly remove files
    installed by the ``conda`` command.

    pip upgrade and uninstall operations only work on packages installed
    via ``pip install``.


WinPython for Windows
-----------------------

The `WinPython <https://winpython.github.io/>`_ project distributes
scikit-learn as an additional plugin.


For installation instructions for particular operating systems or for compiling
the bleeding edge version, see the :ref:`advanced-installation`.
