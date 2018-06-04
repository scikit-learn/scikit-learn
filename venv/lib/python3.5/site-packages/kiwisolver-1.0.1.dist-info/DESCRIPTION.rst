Welcome to Kiwi
===============

.. image:: https://travis-ci.org/nucleic/kiwi.svg?branch=master
    :target: https://travis-ci.org/nucleic/kiwi

Kiwi is an efficient C++ implementation of the Cassowary constraint solving
algorithm. Kiwi is an implementation of the algorithm based on the seminal
Cassowary paper. It is *not* a refactoring of the original C++ solver. Kiwi
has been designed from the ground up to be lightweight and fast. Kiwi ranges
from 10x to 500x faster than the original Cassowary solver with typical use
cases gaining a 40x improvement. Memory savings are consistently > 5x.

In addition to the C++ solver, Kiwi ships with hand-rolled Python bindings.


