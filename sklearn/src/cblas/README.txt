This is a stripped-down version of CBLAS (C-interface to the Basic Linear
Algebra Subroutines), containing only those parts used by scikit-learn's
C/C++/Cython extensions. It is used when no CBLAS implementation is available
at build time.

Sources here are taken from the reference implementation in ATLAS. To add new
algorithms, the only thing that should be done is to copy the reference
implementation from ${ATLAS}/src/blas/reference/level* into this directory.

Header files are taken from ${ATLAS}/include, the only change being the
inclusion of "atlas_refalias*.h" into its respective "atlas_level*.h" file.
