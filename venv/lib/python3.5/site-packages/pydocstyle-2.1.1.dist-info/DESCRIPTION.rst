pydocstyle - docstring style checker
====================================

(formerly pep257)

**pydocstyle** is a static analysis tool for checking compliance with Python
docstring conventions.

**pydocstyle** supports most of
`PEP 257 <http://www.python.org/dev/peps/pep-0257/>`_ out of the box, but it
should not be considered a reference implementation.

**pydocstyle** supports Python 2.7, 3.3, 3.4, 3.5, 3.6 and pypy.

Quick Start
-----------

Install
^^^^^^^

.. code::

    pip install pydocstyle


Run
^^^^

.. code::

    $ pydocstyle test.py
    test.py:18 in private nested class `meta`:
            D101: Docstring missing
    test.py:27 in public function `get_user`:
        D300: Use """triple double quotes""" (found '''-quotes)
    test:75 in public function `init_database`:
        D201: No blank lines allowed before function docstring (found 1)
    ...


Links
-----

.. image:: https://travis-ci.org/PyCQA/pydocstyle.svg?branch=master
    :target: https://travis-ci.org/PyCQA/pydocstyle

.. image:: https://ci.appveyor.com/api/projects/status/40kkc366bmrrttca/branch/master?svg=true
    :target: https://ci.appveyor.com/project/Nurdok/pydocstyle/branch/master

.. image:: https://readthedocs.org/projects/pydocstyle/badge/?version=latest
    :target: https://readthedocs.org/projects/pydocstyle/?badge=latest
    :alt: Documentation Status

* `Read the full documentation here <http://pydocstyle.org>`_.

* `Fork pydocstyle on GitHub <http://github.com/PyCQA/pydocstyle>`_.

* `PyPI project page <https://pypi.python.org/pypi/pydocstyle>`_.


