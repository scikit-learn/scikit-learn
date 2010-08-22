Installing `scikits.learn`
===============================

There are different ways to get scikits.learn installed. 

Linux distributions
-------------------

Some distributions (notably Debian) provide this package in its
distribution. However, depending on your distribution, the package can
be several months old.


Easy install
------------

This is usually the fastest way to install the latest stable
release. If you have pip or easy_install, you can install or update
with the command::

    pip install -U scikits.learn

or::

    easy_install -U scikits.learn

for easy_install.


Windows installer
-----------------

You can download a windows installer from `downloads
<https://sourceforge.net/projects/scikit-learn/files/>`_
in the project's web page.


From Source
-----------
Download the package from http://sourceforge.net/projects/scikit-learn/files
, unpack the sources and cd into archive.

This packages uses distutils, which is the default way of installing
python modules. The install command is::

  python setup.py install


Bleeding Edge
-------------

See section :ref:`git_repo` on how to get the development version.
