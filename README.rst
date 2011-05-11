.. -*- mode: rst -*-

About
=====

scikits.learn is a python module for machine learning built on top of
scipy.

The project was started in 2007 by David Cournapeau as a Google Summer
of Code project, and since then many volunteers have contributed. See
the AUTHORS.rst file for a complete list of contributors.

It is currently maintained by a team of volunteers.


Important links
===============

- Official source code repo: https://github.com/scikit-learn/scikit-learn
- HTML documentation (stable release): http://scikit-learn.sourceforge.net/
- HTML documentation (development version): http://scikit-learn.sourceforge.net/dev/
- Download releases: http://sourceforge.net/projects/scikit-learn/files/
- Issue tracker: https://github.com/scikit-learn/scikit-learn/issues
- Mailing list: https://lists.sourceforge.net/lists/listinfo/scikit-learn-general
- IRC channel: ``#scikit-learn`` at ``irc.freenode.net``

Dependencies
============

The required dependencies to build the software are python >= 2.5,
setuptools, NumPy >= 1.2, SciPy >= 0.7 and a working C++ compiler.

To run the tests you will also need nose >= 0.10.


Install
=======

This packages uses distutils, which is the default way of installing
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

    python -c "import scikits.learn as skl; skl.test()"

See web page http://scikit-learn.sourceforge.net/install.html#testing
for more information.