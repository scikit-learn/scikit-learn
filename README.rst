.. -*- mode: rst -*-

|Travis|_ |AppVeyor|_ |Coveralls|_ |CircleCI|_ |Python27|_ |Python35|_ |PyPi|_ |DOI|_

.. |Travis| image:: https://api.travis-ci.org/scikit-learn/scikit-learn.svg?branch=master
.. _Travis: https://travis-ci.org/scikit-learn/scikit-learn

.. |AppVeyor| image:: https://ci.appveyor.com/api/projects/status/github/scikit-learn/scikit-learn?branch=master&svg=true
.. _AppVeyor: https://ci.appveyor.com/project/sklearn-ci/scikit-learn/history

.. |Coveralls| image:: https://coveralls.io/repos/scikit-learn/scikit-learn/badge.svg?branch=master&service=github
.. _Coveralls: https://coveralls.io/r/scikit-learn/scikit-learn

.. |CircleCI| image:: https://circleci.com/gh/scikit-learn/scikit-learn/tree/master.svg?style=shield&circle-token=:circle-token
.. _CircleCI: https://circleci.com/gh/scikit-learn/scikit-learn

.. |Python27| image:: https://img.shields.io/badge/python-2.7-blue.svg
.. _Python27: https://badge.fury.io/py/scikit-learn

.. |Python35| image:: https://img.shields.io/badge/python-3.5-blue.svg
.. _Python35: https://badge.fury.io/py/scikit-learn

.. |PyPi| image:: https://badge.fury.io/py/scikit-learn.svg
.. _PyPi: https://badge.fury.io/py/scikit-learn

.. |DOI| image:: https://zenodo.org/badge/21369/scikit-learn/scikit-learn.svg
.. _DOI: https://zenodo.org/badge/latestdoi/21369/scikit-learn/scikit-learn

scikit-learn
============

scikit-learn is a Python module for machine learning built on top of
SciPy and distributed under the 3-Clause BSD license.

The project was started in 2007 by David Cournapeau as a Google Summer
of Code project, and since then many volunteers have contributed. See
the `AUTHORS.rst <AUTHORS.rst>`_ file for a complete list of contributors.

It is currently maintained by a team of volunteers.

Website: http://scikit-learn.org

Installation
------------

Dependencies
~~~~~~~~~~~~

scikit-learn requires:

- Python (>= 2.6 or >= 3.3)
- NumPy (>= 1.6.1)
- SciPy (>= 0.9)

scikit-learn also uses CBLAS, the C interface to the Basic Linear Algebra
Subprograms library. scikit-learn comes with a reference implementation, but
the system CBLAS will be detected by the build system and used if present.
CBLAS exists in many implementations; see `Linear algebra libraries
<http://scikit-learn.org/stable/modules/computational_performance.html#linear-algebra-libraries>`_
for known issues.

User installation
~~~~~~~~~~~~~~~~~

If you already have a working installation of numpy and scipy,
the easiest way to install scikit-learn is using ``pip`` ::

    pip install -U scikit-learn

or ``conda``::

    conda install scikit-learn

The documentation includes more detailed `installation instructions <http://scikit-learn.org/stable/install.html>`_.


Development
-----------

We welcome new contributors of all experience levels. The scikit-learn
community goals are to be helpful, welcoming, and effective. The
`Development Guide <http://scikit-learn.org/stable/developers/index.html>`_
has detailed information about contributing code, documentation, tests, and
more. We've included some basic information in this README.

Important links
~~~~~~~~~~~~~~~

- Official source code repo: https://github.com/scikit-learn/scikit-learn
- Download releases: https://pypi.python.org/pypi/scikit-learn
- Issue tracker: https://github.com/scikit-learn/scikit-learn/issues

Source code
~~~~~~~~~~~

You can check the latest sources with the command::

    git clone https://github.com/scikit-learn/scikit-learn.git

Setting up a development environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quick tutorial on how to go about setting up your environment to
contribute to scikit-learn: https://github.com/scikit-learn/scikit-learn/blob/master/CONTRIBUTING.md

Testing
~~~~~~~

After installation, you can launch the test suite from outside the
source directory (you will need to have the ``nose`` package installed)::

    nosetests -v sklearn

Under Windows, it is recommended to use the following command (adjust the path
to the ``python.exe`` program) as using the ``nosetests.exe`` program can badly
interact with tests that use ``multiprocessing``::

    C:\Python34\python.exe -c "import nose; nose.main()" -v sklearn

See the web page http://scikit-learn.org/stable/install.html#testing
for more information.

    Random number generation can be controlled during testing by setting
    the ``SKLEARN_SEED`` environment variable.

Submitting a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~

Before opening a Pull Request, have a look at the
full Contributing page to make sure your code complies
with our guidelines: http://scikit-learn.org/stable/developers/index.html


Project History
---------------

The project was started in 2007 by David Cournapeau as a Google Summer
of Code project, and since then many volunteers have contributed. See
the AUTHORS.rst file for a complete list of contributors.

The project is currently maintained by a team of volunteers.

**Note**: `scikit-learn` was previously referred to as `scikits.learn`.


Help and Support
----------------

Documentation
~~~~~~~~~~~~~

- HTML documentation (stable release): http://scikit-learn.org
- HTML documentation (development version): http://scikit-learn.org/dev/
- FAQ: http://scikit-learn.org/stable/faq.html

Communication
~~~~~~~~~~~~~

- Mailing list: https://mail.python.org/mailman/listinfo/scikit-learn
- IRC channel: ``#scikit-learn`` at ``irc.freenode.net``
- Stack Overflow: http://stackoverflow.com/questions/tagged/scikit-learn
- Website: http://scikit-learn.org

Citation
~~~~~~~~

If you use scikit-learn in a scientific publication, we would appreciate citations: http://scikit-learn.org/stable/about.html#citing-scikit-learn
