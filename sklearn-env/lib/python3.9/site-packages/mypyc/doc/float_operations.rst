.. _float-ops:

Native float operations
========================

These ``float`` operations have fast, optimized implementations. Other
floating point operations use generic implementations that are often
slower.

.. note::

   At the moment, only a few float operations are optimized. This will
   improve in future mypyc releases.

Construction
------------

* Float literal
* ``float(string)``

Functions
---------

* ``abs(f)``
