cdef extern from "xsimd.hpp" namespace "xsimd::detail":
    ctypedef struct supported_arch:
        unsigned int sse2
        unsigned int sse3
        unsigned int ssse3
        unsigned int sse4_1
        unsigned int sse4_2
        unsigned int sse4a
        unsigned int fma3_sse
        unsigned int fma4
        unsigned int xop
        unsigned int avx
        unsigned int fma3_avx
        unsigned int avx2
        unsigned int fma3_avx2
        unsigned int avx512f
        unsigned int avx512cd
        unsigned int avx512dq
        unsigned int avx512bw
        unsigned int neon
        unsigned int neon64
        unsigned int sve

cdef extern from "xsimd.hpp" namespace "xsimd":
    cdef supported_arch available_architectures() noexcept nogil

cpdef get_available_architectures():
    return available_architectures()
