.. _installation-instructions:

=======================
Installing scikit-learn
=======================

.. note::

    If you wish to contribute to the project, it's recommended you
    :ref:`install the latest development version<install_bleeding_edge>`.

There are different ways to get scikit-learn installed:

  * Install the version of scikit-learn provided by your
    :ref:`operating system or Python distribution <install_by_distribution>`.
    This is a quick option for those who have operating systems
    that distribute scikit-learn. It might not provide the latest release
    version.

  * :ref:`Install an official release <install_official_release>`. This
    is the best approach for most users. It will provide a stable version
    and pre-build packages are available for most platforms.

  * :ref:`Building the package from source
    <install_bleeding_edge>`. This is best for users who want the
    latest-and-greatest features and aren't afraid of running
    brand-new code.


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
the easiest way to install scikit-learn is using ``conda``::

    conda install scikit-learn

or ``pip``.
In that case, in order to avoid any OS dependent issue it is strongly
recommended to use python3 ``virtualenv``
(see `python3 virtualenv documentation <https://docs.python.org/3/tutorial/venv.html>`_).

That gives, on Linux::

    python -m venv .myenv
    source .myenv/bin/activate
    pip install -U scikit-learn

or on Windows::

    python -m venv .myenv
    .myenv\Scripts\activate
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

