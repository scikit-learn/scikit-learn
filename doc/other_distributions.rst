
.. _install_by_distribution:

Third party distributions of scikit-learn
=========================================

Some third-party distributions are now providing versions of
scikit-learn integrated with their package-management systems.
The most popular ones are listed on the :ref:`install` page.

These can make installation and upgrading much easier for users since
the integration includes the ability to automatically install
dependencies (numpy, scipy) that scikit-learn requires.

The following is an incomplete list of python and os distributions
that provide their own version of scikit-learn.


MacPorts for Mac OSX
--------------------

The MacPorts package is named ``py<XY>-scikits-learn``,
where ``XY`` denotes the Python version.
It can be installed by typing the following
command::

    sudo port install py27-scikit-learn

or::

    sudo port install py36-scikit-learn


Arch Linux
----------

Arch Linux's package is provided through the `official repositories
<https://www.archlinux.org/packages/?q=scikit-learn>`_ as
``python-scikit-learn`` for Python.
It can be installed by typing the following command:

.. code-block:: none

     # pacman -S python-scikit-learn



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

