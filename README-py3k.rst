Support for Python 3
====================

There is experimental support for python 3. To generate python3
compatible sources, run the 2to3 tool on the sources top directory::

    2to3 -wn --no-diffs scikits examples

To execute some examples and run the test suite You will also need
some other tools that as of July 2010 only have experimental python3
support, like scipy and nosetest3.
