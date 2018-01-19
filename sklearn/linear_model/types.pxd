# Synopsis: Template declarations of some fundamental types
#
# Author: Elvis Dohmatob <gmdopp@gmail.com>
# License: BSD

cimport numpy as np

cdef fused floating:
    float
    double

cdef fused complexing:
    # real
    floating

    # complex
    float complex
    double complex
