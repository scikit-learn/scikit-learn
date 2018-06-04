The Cython language makes writing C extensions for the Python language as
easy as Python itself.  Cython is a source code translator based on Pyrex_,
but supports more cutting edge functionality and optimizations.

The Cython language is a superset of the Python language (almost all Python
code is also valid Cython code), but Cython additionally supports optional
static typing to natively call C functions, operate with C++ classes and
declare fast C types on variables and class attributes.  This allows the
compiler to generate very efficient C code from Cython code.

This makes Cython the ideal language for writing glue code for external
C/C++ libraries, and for fast C modules that speed up the execution of
Python code.

Note that for one-time builds, e.g. for CI/testing, on platforms that are not
covered by one of the wheel packages provided on PyPI, it is substantially faster
than a full source build to install an uncompiled (slower) version of Cython with::

    pip install Cython --install-option="--no-cython-compile"

.. _Pyrex: http://www.cosc.canterbury.ac.nz/greg.ewing/python/Pyrex/


