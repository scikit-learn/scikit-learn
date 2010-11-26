
This is a stripped-down version of CBLAS.

Sources are taken from ATLAS 3.9.25. To avoid ATLAS' build system, we
only ship the reference implementation for the algorithms we need,
which can be found at
$(prefix)/ATLAS/src/blas/reference/level*. Header files are taken from
$(prefix)/ATLAS/include, the only change being the inclusion of
"atlas_refalias*.h" into its respective "atlas_level*.h" file.

To add new algorithms, the only thing that should be done is to copy
the reference implementation from
$(prefix)/ATLAS/src/blas/reference/level* into this directory.
