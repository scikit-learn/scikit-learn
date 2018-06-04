NumPy is a general-purpose array-processing package designed to
efficiently manipulate large multi-dimensional arrays of arbitrary
records without sacrificing too much speed for small multi-dimensional
arrays.  NumPy is built on the Numeric code base and adds features
introduced by numarray as well as an extended C-API and the ability to
create arrays of arbitrary type which also makes NumPy suitable for
interfacing with general-purpose data-base applications.

There are also basic facilities for discrete fourier transform,
basic linear algebra and random number generation.

All numpy wheels distributed from pypi are BSD licensed.

Windows wheels are linked against the ATLAS BLAS / LAPACK library, restricted
to SSE2 instructions, so may not give optimal linear algebra performance for
your machine. See http://docs.scipy.org/doc/numpy/user/install.html for
alternatives.



