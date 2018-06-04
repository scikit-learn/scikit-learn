.. image:: https://img.shields.io/travis/jquast/wcwidth.svg
    :target: https://travis-ci.org/jquast/wcwidth
    :alt: Travis Continous Integration

.. image:: https://img.shields.io/coveralls/jquast/wcwidth.svg
    :target: https://coveralls.io/r/jquast/wcwidth
    :alt: Coveralls Code Coverage

.. image:: https://img.shields.io/pypi/v/wcwidth.svg
    :target: https://pypi.python.org/pypi/wcwidth/
    :alt: Latest Version

.. image:: https://img.shields.io/github/license/jquast/wcwidth.svg
    :target: https://pypi.python.org/pypi/wcwidth/
    :alt: License

.. image:: https://img.shields.io/pypi/wheel/wcwidth.svg
    :alt: Wheel Status

.. image:: https://img.shields.io/pypi/dm/wcwidth.svg
    :target: https://pypi.python.org/pypi/wcwidth/
    :alt: Downloads

============
Introduction
============

This Library is mainly for those implementing a Terminal Emulator, or programs
that carefully produce output to be interpreted by one.

**Problem Statement**: When printed to the screen, the length of the string is
usually equal to the number of cells it occupies.  However, there are
categories of characters that occupy 2 cells (full-wide), and others that
occupy 0.


**Solution**: POSIX.1-2001 and POSIX.1-2008 conforming systems provide
`wcwidth(3)`_ and `wcswidth(3)`_ C functions of which this python module's
functions precisely copy.  *These functions return the number of cells a
unicode string is expected to occupy.*

This library aims to be forward-looking, portable, and most correct.  The most
current release of this API is based on the Unicode Standard release files:

``DerivedGeneralCategory-9.0.0.txt``
  *Date: 2016-06-01, 10:34:26 GMT*
  © 2016 Unicode®, Inc.

``EastAsianWidth-9.0.0.txt``
  *Date: 2016-05-27, 17:00:00 GMT [KW, LI]*
  © 2016 Unicode®, Inc.


Installation
------------

The stable version of this package is maintained on pypi, install using pip::

    pip install wcwidth

Example
-------

To Display ``u'コンニチハ'`` right-adjusted on screen of 80 columns::

    >>> from wcwidth import wcswidth
    >>> text = u'コンニチハ'
    >>> text_len = wcswidth(text)
    >>> print(u' ' * (80 - text_len) + text)

wcwidth, wcswidth
-----------------
Use function ``wcwidth()`` to determine the length of a *single unicode
character*, and ``wcswidth()`` to determine the length of a several, or a
*string of unicode characters*.

Briefly, return values of function ``wcwidth()`` are:

``-1``
  Indeterminate (not printable).

``0``
  Does not advance the cursor, such as NULL or Combining.

``2``
  Characters of category East Asian Wide (W) or East Asian
  Full-width (F) which are displayed using two terminal cells.

``1``
  All others.

Function ``wcswidth()`` simply returns the sum of all values for each character
along a string, or ``-1`` when it occurs anywhere along a string.

More documentation is available using pydoc::

    $ pydoc wcwidth

=======
Caveats
=======

This library attempts to determine the printable width by an unknown targeted
terminal emulator.  It does not provide any ability to discern what the target
emulator software, version, of level of support is.  Results may vary!

A `crude method
<http://blessed.readthedocs.org/en/latest/examples.html#detect-multibyte-py>`_
of determining the level of unicode support by the target emulator may be
performed using the VT100 Query Cursor Position sequence.

The libc version of `wcwidth(3)`_ is often several unicode releases behind,
and therefor several levels of support lower than this python library.  You
may determine an exacting list of these discrepancies using the project
file `wcwidth-libc-comparator.py
<https://github.com/jquast/wcwidth/tree/master/bin/wcwidth-libc-comparator.py>`_.


==========
Developing
==========

Install wcwidth in editable mode::

   pip install -e.

Install developer requirements::

   pip install -r requirements-develop.txt

Execute unit tests using tox::

   tox

Updating Tables
---------------

The command ``python setup.py update`` will fetch the following resources:

- http://www.unicode.org/Public/UNIDATA/EastAsianWidth.txt
- http://www.unicode.org/Public/UNIDATA/extracted/DerivedGeneralCategory.txt

And generates the table files:

- `wcwidth/table_wide.py <https://github.com/jquast/wcwidth/tree/master/wcwidth/table_wide.py>`_
- `wcwidth/table_zero.py <https://github.com/jquast/wcwidth/tree/master/wcwidth/table_zero.py>`_

Uses
----

This library is used in:

- `jquast/blessed`_, a simplified wrapper around curses.

- `jonathanslenders/python-prompt-toolkit`_, a Library for building powerful
  interactive command lines in Python.

Additional tools for displaying and testing wcwidth are found in the `bin/
<https://github.com/jquast/wcwidth/tree/master/bin>`_ folder of this project's
source code.  They are not distributed.

=======
History
=======

0.1.7 *2016-07-01*
  * **Updated** tables to Unicode Specification 9.0.0. (`PR #18`_).

0.1.6 *2016-01-08 Production/Stable*
  * ``LICENSE`` file now included with distribution.

0.1.5 *2015-09-13 Alpha*
  * **Bugfix**:
    Resolution of "combining_ character width" issue, most especially
    those that previously returned -1 now often (correctly) return 0.
    resolved by `Philip Craig`_ via `PR #11`_.
  * **Deprecated**:
    The module path ``wcwidth.table_comb`` is no longer available,
    it has been superseded by module path ``wcwidth.table_zero``.

0.1.4 *2014-11-20 Pre-Alpha*
  * **Feature**: ``wcswidth()`` now determines printable length
    for (most) combining_ characters.  The developer's tool
    `bin/wcwidth-browser.py`_ is improved to display combining_
    characters when provided the ``--combining`` option
    (`Thomas Ballinger`_ and `Leta Montopoli`_ `PR #5`_).
  * **Feature**: added static analysis (prospector_) to testing
    framework.

0.1.3 *2014-10-29 Pre-Alpha*
  * **Bugfix**: 2nd parameter of wcswidth was not honored.
    (`Thomas Ballinger`_, `PR #4`_).

0.1.2 *2014-10-28 Pre-Alpha*
  * **Updated** tables to Unicode Specification 7.0.0.
    (`Thomas Ballinger`_, `PR #3`_).

0.1.1 *2014-05-14 Pre-Alpha*
  * Initial release to pypi, Based on Unicode Specification 6.3.0

This code was originally derived directly from C code of the same name,
whose latest version is available at
http://www.cl.cam.ac.uk/~mgk25/ucs/wcwidth.c::

 * Markus Kuhn -- 2007-05-26 (Unicode 5.0)
 *
 * Permission to use, copy, modify, and distribute this software
 * for any purpose and without fee is hereby granted. The author
 * disclaims all warranties with regard to this software.

.. _`prospector`: https://github.com/landscapeio/prospector
.. _`combining`: https://en.wikipedia.org/wiki/Combining_character
.. _`bin/wcwidth-browser.py`: https://github.com/jquast/wcwidth/tree/master/bin/wcwidth-browser.py
.. _`Thomas Ballinger`: https://github.com/thomasballinger
.. _`Leta Montopoli`: https://github.com/lmontopo
.. _`Philip Craig`: https://github.com/philipc
.. _`PR #3`: https://github.com/jquast/wcwidth/pull/3
.. _`PR #4`: https://github.com/jquast/wcwidth/pull/4
.. _`PR #5`: https://github.com/jquast/wcwidth/pull/5
.. _`PR #11`: https://github.com/jquast/wcwidth/pull/11
.. _`PR #18`: https://github.com/jquast/wcwidth/pull/18
.. _`jquast/blessed`: https://github.com/jquast/blessed
.. _`jonathanslenders/python-prompt-toolkit`: https://github.com/jonathanslenders/python-prompt-toolkit
.. _`wcwidth(3)`:  http://man7.org/linux/man-pages/man3/wcwidth.3.html
.. _`wcswidth(3)`: http://man7.org/linux/man-pages/man3/wcswidth.3.html


