Install
=======


Binary Packages
---------------

There are currently no binary packages available. If you would like to
make one, please see section :ref:`packaging <packaging>`.


From Source
-----------

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
