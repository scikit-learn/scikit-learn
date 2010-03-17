Install
=======

Binary Packages
---------------

There is a prebuild package for windows. See section downloads in the
project's web page.


From Source
-----------
Download the package from http://sourceforge.net/projects/scikit-learn/files
, unpack the sources and cd into archive.

This packages uses distutils, which is the default way of installing
python modules. The install command is::

  python setup.py install

If you have installed the boost libraries in a non-standard location
you might need to pass the appropriate --include argument so that it
find the correct headers. For example, if your headers reside in
/opt/local/include, (which is the case if you have installed them
through Mac Ports), you must issue the commands::

  python setup.py build_ext --include=/opt/local/include
  python setup.py install


Bleeding Edge
-------------

Latests source code can be found on subversion::

  svn co http://scikit-learn.svn.sourceforge.net/svnroot/scikit-learn/trunk scikit-learn
