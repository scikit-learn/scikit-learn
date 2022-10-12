from ..utils._typedefs cimport ITYPE_t

cdef class UnionFind:
    cdef ITYPE_t next_label
    cdef ITYPE_t[:] parent
    cdef ITYPE_t[:] size

    cdef void union(self, ITYPE_t m, ITYPE_t n)
    cdef ITYPE_t fast_find(self, ITYPE_t n)
