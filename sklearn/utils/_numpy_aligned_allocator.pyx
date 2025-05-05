cdef extern from "src/numpy_allocator.h":
    void c_set_aligned_allocation()
    void c_set_numpy_default_allocation()


def set_aligned_allocation():
    c_set_aligned_allocation()


def set_numpy_default_allocation():
    c_set_numpy_default_allocation()
