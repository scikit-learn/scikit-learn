.. -*- mode: rst -*-

About
=====

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

The required dependencies to build the software are Python >= 2.6,
setuptools, Numpy >= 1.3, SciPy >= 0.7 and a working C/C++ compiler.
This configuration matches the Ubuntu 10.04 LTS release from April 2010.

To run the tests you will also need nose >= 0.10.


Install
=======

This package uses distutils, which is the default way of installing
python modules. To install in your home directory, use::

  python setup.py install --home

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
source directory (you will need to have nosetest installed)::

    python -c "import sklearn; sklearn.test()"

See web page http://scikit-learn.sourceforge.net/install.html#testing
for more information.
