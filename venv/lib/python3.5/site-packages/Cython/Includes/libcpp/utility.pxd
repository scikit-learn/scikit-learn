cdef extern from "<utility>" namespace "std" nogil:
    cdef cppclass pair[T, U]:
        ctypedef T first_type
        ctypedef U second_type
        T first
        U second
        pair() except +
        pair(pair&) except +
        pair(T&, U&) except +
        bint operator==(pair&, pair&)
        bint operator!=(pair&, pair&)
        bint operator<(pair&, pair&)
        bint operator>(pair&, pair&)
        bint operator<=(pair&, pair&)
        bint operator>=(pair&, pair&)
