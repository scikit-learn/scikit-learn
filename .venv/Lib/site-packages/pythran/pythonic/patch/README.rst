=================
Pythran std patch
=================

The implementation of std::complex is very slow, due to the complex type not
implementing the limited range (see `-fcx-limited-range`) optimization. Numpy
does implement it, so we have to conform to numpy's version. the only way I
(SG) found to fix this is to monkey-patch `std::complex`. Inheritance or
defining a new class does not work because nt2 assumes we use std::complex.

The original source is libcxx, the diff is rather small (mostly removed libcxx
internal stuff and use numpy-compliant version of the multiply operator). The
speedup is very interesting!

GCC does provide the flag `-fcx-limited-range` to fix the issue in a more
elegant way, but it is not supported by clang.

The CPython impelmentation for complex division can be found in
``Objects/complexobject.c``, and the numpy implementation lies in
``numpy/core/src/umath/loops.c.src``, for those interested.
