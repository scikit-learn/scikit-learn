from ._typedefs cimport float32_t, float64_t, intp_t, uint8_t, uint32_t, int32_t

ctypedef float64_t X_DTYPE_C
ctypedef uint8_t X_BINNED_DTYPE_C
ctypedef float64_t Y_DTYPE_C
ctypedef float32_t G_H_DTYPE_C
ctypedef uint32_t BITSET_INNER_DTYPE_C
ctypedef BITSET_INNER_DTYPE_C[8] BITSET_DTYPE_C

# TODO: add templating to allow each function to operate on uint8_t (X_BINNED_DTYPE_C), or some other types
# which would allow larger category sets in decision-tree splitting.

cdef void init_bitset(BITSET_DTYPE_C bitset) noexcept nogil

cdef void set_bitset(BITSET_DTYPE_C bitset, X_BINNED_DTYPE_C val) noexcept nogil

                                
cdef uint8_t in_bitset(BITSET_DTYPE_C bitset, X_BINNED_DTYPE_C val) noexcept nogil

cdef void set_bitset_memoryview(BITSET_INNER_DTYPE_C[:] bitset,  # OUT
                                X_BINNED_DTYPE_C val) noexcept nogil

cdef uint8_t in_bitset_memoryview(const BITSET_INNER_DTYPE_C[:] bitset,
                                  X_BINNED_DTYPE_C val) noexcept nogil

cpdef uint8_t py_in_bitset_memoryview(const BITSET_INNER_DTYPE_C[:] bitset,
                                      X_BINNED_DTYPE_C val) noexcept nogil

cdef uint8_t in_bitset_2d_memoryview(
    const BITSET_INNER_DTYPE_C[:, :] bitset,
    X_BINNED_DTYPE_C val,
    unsigned int row) noexcept nogil


cdef uint8_t set_bitset_2d_memoryview(
    const BITSET_INNER_DTYPE_C[:, :] bitset,
    X_BINNED_DTYPE_C val,
    unsigned int row) noexcept nogil




cdef void flip_bitset_up_to_nbits(
    BITSET_INNER_DTYPE_C[:] bitset,
    int32_t n_categories
) noexcept nogil

cdef void bs_from_template_memoryview(
    BITSET_INNER_DTYPE_C[:] bitset,   # OUT: a length‐8 uint32_t memoryview
    BITSET_INNER_DTYPE_C   template,  # template is now a 32-bit mask over “present categories”
    intp_t[:]             cat_offs,   # e.g. [2, 5, 7] if you have 3 present categories
    intp_t                ncats_present
) noexcept nogil