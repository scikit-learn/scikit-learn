.. _int-ops:

Native integer operations
=========================

Operations on ``int`` values that are listed here have fast, optimized
implementations. Other integer operations use generic implementations
that are often slower. Some operations involving integers and other
types are documented elsewhere, such as list indexing.

Construction
------------

* Integer literal
* ``int(x: float)``
* ``int(x: str)``
* ``int(x: str, base: int)``

Operators
---------

* Arithmetic (``+``, ``-``, ``*``, ``//``, ``%``)
* Bitwise operations (``&``, ``|``, ``^``, ``<<``, ``>>``, ``~``)
* Comparisons (``==``, ``!=``, ``<``, etc.)
* Augmented assignment (``x += y``, etc.)

Statements
----------

For loop over range:

* ``for x in range(end)``
* ``for x in range(start, end)``
* ``for x in range(start, end, step)``
