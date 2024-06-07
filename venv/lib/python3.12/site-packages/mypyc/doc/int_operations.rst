.. _int-ops:

Native integer operations
=========================

Mypyc supports these integer types:

* ``int`` (arbitrary-precision integer)
* ``i64`` (64-bit signed integer)
* ``i32`` (32-bit signed integer)
* ``i16`` (16-bit signed integer)
* ``u8`` (8-bit unsigned integer)

``i64``, ``i32``, ``i16`` and ``u8`` are *native integer types* and
are available in the ``mypy_extensions`` module. ``int`` corresponds
to the Python ``int`` type, but uses a more efficient runtime
representation (tagged pointer). Native integer types are value types.

All integer types have optimized primitive operations, but the native
integer types are more efficient than ``int``, since they don't
require range or bounds checks.

Operations on integers that are listed here have fast, optimized
implementations. Other integer operations use generic implementations
that are generally slower. Some operations involving integers and other
types, such as list indexing, are documented elsewhere.

Construction
------------

``int`` type:

* Integer literal
* ``int(x: float)``
* ``int(x: i64)``
* ``int(x: i32)``
* ``int(x: i16)``
* ``int(x: u8)``
* ``int(x: str)``
* ``int(x: str, base: int)``
* ``int(x: int)`` (no-op)

``i64`` type:

* ``i64(x: int)``
* ``i64(x: float)``
* ``i64(x: i64)`` (no-op)
* ``i64(x: i32)``
* ``i64(x: i16)``
* ``i64(x: u8)``
* ``i64(x: str)``
* ``i64(x: str, base: int)``

``i32`` type:

* ``i32(x: int)``
* ``i32(x: float)``
* ``i32(x: i64)`` (truncate)
* ``i32(x: i32)`` (no-op)
* ``i32(x: i16)``
* ``i32(x: u8)``
* ``i32(x: str)``
* ``i32(x: str, base: int)``

``i16`` type:

* ``i16(x: int)``
* ``i16(x: float)``
* ``i16(x: i64)`` (truncate)
* ``i16(x: i32)`` (truncate)
* ``i16(x: i16)`` (no-op)
* ``i16(x: u8)``
* ``i16(x: str)``
* ``i16(x: str, base: int)``

Conversions from ``int`` to a native integer type raise
``OverflowError`` if the value is too large or small. Conversions from
a wider native integer type to a narrower one truncate the value and never
fail. More generally, operations between native integer types don't
check for overflow.

Implicit conversions
--------------------

``int`` values can be implicitly converted to a native integer type,
for convenience. This means that these are equivalent::

   from mypy_extensions import i64

   def implicit() -> None:
       # Implicit conversion of 0 (int) to i64
       x: i64 = 0

   def explicit() -> None:
       # Explicit conversion of 0 (int) to i64
       x = i64(0)

Similarly, a native integer value can be implicitly converted to an
arbitrary-precision integer. These two functions are equivalent::

   def implicit(x: i64) -> int:
       # Implicit conversion from i64 to int
       return x

   def explicit(x: i64) -> int:
       # Explicit conversion from i64 to int
       return int(x)

Operators
---------

* Arithmetic (``+``, ``-``, ``*``, ``//``, ``/``, ``%``)
* Bitwise operations (``&``, ``|``, ``^``, ``<<``, ``>>``, ``~``)
* Comparisons (``==``, ``!=``, ``<``, etc.)
* Augmented assignment (``x += y``, etc.)

If one of the above native integer operations overflows or underflows
with signed operands, the behavior is undefined. Signed native integer
types should only be used if all possible values are small enough for
the type. For this reason, the arbitrary-precision ``int`` type is
recommended for signed values unless the performance of integer
operations is critical.

Operations on unsigned integers (``u8``) wrap around on overflow.

It's a compile-time error to mix different native integer types in a
binary operation such as addition. An explicit conversion is required::

    from mypy_extensions import i64, i32

    def add(x: i64, y: i32) -> None:
        a = x + y  # Error (i64 + i32)
        b = x + i64(y)  # OK

You can freely mix a native integer value and an arbitrary-precision
``int`` value in an operation. The native integer type is "sticky"
and the ``int`` operand is coerced to the native integer type::

  def example(x: i64, y: int) -> None:
      a = x * y
      # Type of "a" is "i64"
      ...
      b = 1 - x
      # Similarly, type of "b" is "i64"

Statements
----------

For loop over a range is compiled efficiently, if the ``range(...)`` object
is constructed in the for statement (after ``in``):

* ``for x in range(end)``
* ``for x in range(start, end)``
* ``for x in range(start, end, step)``

If one of the arguments to ``range`` in a for loop is a native integer
type, the type of the loop variable is inferred to have this native
integer type, instead of ``int``::

  for x in range(i64(n)):
      # Type of "x" is "i64"
      ...
