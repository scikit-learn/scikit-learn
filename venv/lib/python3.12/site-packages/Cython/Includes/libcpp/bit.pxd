cdef extern from "<bit>" namespace "std" nogil:
    # bit_cast (gcc >= 11.0, clang >= 14.0)
    cdef To bit_cast[To, From](From&)

    # byteswap (C++23)
    #cdef T byteswap[T](T)

    # integral powers of 2 (gcc >= 10.0, clang >= 12.0)
    cdef bint has_single_bit[T](T)
    cdef T bit_ceil[T](T)
    cdef T bit_floor[T](T)
    cdef int bit_width[T](T)

    # rotating (gcc >= 9.0, clang >= 9.0)
    cdef T rotl[T](T, int shift)
    cdef T rotr[T](T, int shift)

    # counting (gcc >= 9.0, clang >= 9.0)
    cdef int countl_zero[T](T)
    cdef int countl_one[T](T)
    cdef int countr_zero[T](T)
    cdef int countr_one[T](T)
    cdef int popcount[T](T)

    # endian
    cpdef enum class endian(int):
        little,
        big,
        native
