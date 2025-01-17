from libcpp cimport bool
from libcpp.typeinfo cimport type_info

cdef extern from "<any>" namespace "std" nogil:
    cdef cppclass any:
        any()
        any(any&) except +
        void reset()
        bool has_value()
        type_info& type()
        T& emplace[T](...) except +
        void swap(any&)
        any& operator=(any&) except +
        any& operator=[U](U&) except +

    cdef T any_cast[T](any&) except +
