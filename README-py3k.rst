Support for Python 3
====================

There is experimental support for python3 in some modules, the status
of these is:

  It mostly works (it builds) but there are still failing tests and doctests,
  some of which are caused by changes in the way deprecation warnings are
  handled.

This codebase should work both with Python 3.2+ and 2.6+ without using 2to3.
There is a copy of the `six` helper library in `sklearn.externals.six` to be
used when needed.

If you would like to help with porting to python3, please propose
yourself in the scikit-learn mailing list:

    https://lists.sourceforge.net/lists/listinfo/scikit-learn-general
