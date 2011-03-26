.. -*- mode: rst -*-

About
=====

scikits.learn is a python module for machine learning built on top of
scipy.

The project was started in 2007 by David Cournapeau as a Google Summer
of Code project, and since then many volunteers have contributed. See
the AUTHORS.rst file for a complete list of contributors.

It is currently maintained by a team of volunteers.


Download
========

You can download source code and Windows binaries from SourceForge:

http://sourceforge.net/projects/scikit-learn/files/


Dependencies
============

The required dependencies to build the software are python >= 2.5,
setuptools, NumPy >= 1.2, SciPy >= 0.7 and a working C++ compiler.

To run the tests you will also need nose >= 0.10.


Install
=======

This packages uses distutils, which is the default way of installing
python modules. The install command is::

  python setup.py install


Mailing list
============

There's a general and development mailing list, visit
https://lists.sourceforge.net/lists/listinfo/scikit-learn-general to
subscribe to the mailing list.


IRC channel
===========

Some developers tend to hang around the channel ``#scikit-learn``
at ``irc.freenode.net``, especially during the week preparing a new
release. If nobody is available to answer your questions there don't
hesitate to ask it on the mailing list to reach a wider audience.


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

Bugs
----

Please submit bugs you might encounter, as well as patches and feature
requests to the tracker located at github
https://github.com/scikit-learn/scikit-learn/issues


Testing
-------

After installation, you can launch the test suite from outside the
source directory (you will need to have nosetest installed)::

    python -c "import scikits.learn as skl; skl.test()"

See web page http://scikit-learn.sourceforge.net/install.html#testing
for more information.



