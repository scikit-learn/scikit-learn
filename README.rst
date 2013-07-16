.. -*- mode: rst -*-

|Travis|_

.. |Travis| image:: https://api.travis-ci.org/scikit-learn/scikit-learn.png?branch=master
.. _Travis: https://travis-ci.org/scikit-learn/scikit-learn

scikit-learn
============

scikit-learn is a Python module for machine learning built on top of
SciPy and distributed under the 3-Clause BSD license.

The project was started in 2007 by David Cournapeau as a Google Summer
of Code project, and since then many volunteers have contributed. See
the AUTHORS.rst file for a complete list of contributors.

It is currently maintained by a team of volunteers.

**Note** `scikit-learn` was previously referred to as `scikits.learn`.


Important links
===============

- Official source code repo: https://github.com/scikit-learn/scikit-learn
- HTML documentation (stable release): http://scikit-learn.org
- HTML documentation (development version): http://scikit-learn.org/dev/
- Download releases: http://sourceforge.net/projects/scikit-learn/files/
- Issue tracker: https://github.com/scikit-learn/scikit-learn/issues
- Mailing list: https://lists.sourceforge.net/lists/listinfo/scikit-learn-general
- IRC channel: ``#scikit-learn`` at ``irc.freenode.net``

Dependencies
============

scikit-learn is tested to work under Python 2.6+ and Python 3.3+
(using the same codebase thanks to an embedded copy of [six](
http://pythonhosted.org/six/)).

The required dependencies to build the software Numpy >= 1.3, SciPy >= 0.7
and a working C/C++ compiler.

For running the examples Matplotlib >= 0.99.1 is required and for running the
tests you need nose >= 0.10.

This configuration matches the Ubuntu 10.04 LTS release from April 2010.


Install
=======

This package uses distutils, which is the default way of installing
python modules. To install in your home directory, use::

  python setup.py install --user

To install for all users on Unix/Linux::

  python setup.py build
  sudo python setup.py install


Development
===========

Code
----

GIT
~~~

You can check the latest sources with the command::

    git clone git://github.com/scikit-learn/scikit-learn.git

or if you have write privileges::

    git clone git@github.com:scikit-learn/scikit-learn.git


Testing
-------

After installation, you can launch the test suite from outside the
source directory (you will need to have nosetests installed)::

   $ nosetests --exe sklearn

See the web page http://scikit-learn.org/stable/install.html#testing
for more information.

    Random number generation can be controlled during testing by setting
    the ``SKLEARN_SEED`` environment variable.
