from sklearn.utils._typedefs cimport (
    float32_t, float64_t, int32_t, intp_t, uint8_t, uint16_t, uint32_t
)


ctypedef float64_t X_DTYPE_C
ctypedef uint16_t X_BINNED_DTYPE_C
ctypedef fused X_BINNED_DTYPE_FUSED_C:
    uint8_t
    uint16_t
ctypedef float64_t Y_DTYPE_C
ctypedef float32_t G_H_DTYPE_C
ctypedef uint32_t BITSET_INNER_DTYPE_C


cdef packed struct hist_struct:
    # Same as histogram dtype but we need a struct to declare views. It needs
    # to be packed since by default numpy dtypes aren't aligned
    Y_DTYPE_C sum_gradients
    Y_DTYPE_C sum_hessians
    unsigned int count


cdef packed struct node_struct:
    # Equivalent struct to PREDICTOR_RECORD_DTYPE to use in memory views. It
    # needs to be packed since by default numpy dtypes aren't aligned
    Y_DTYPE_C value
    unsigned int count
    intp_t feature_idx
    X_DTYPE_C num_threshold
    unsigned char missing_go_to_left
    unsigned int left
    unsigned int right
    Y_DTYPE_C gain
    unsigned int depth
    unsigned char is_leaf
    X_BINNED_DTYPE_C bin_threshold
    unsigned char is_categorical
    # The index of the corresponding bitsets in the Predictor's bitset arrays.
    # Only used if is_categorical is True
    unsigned int bitset_idx

cpdef enum MonotonicConstraint:
    NO_CST = 0
    POS = 1
    NEG = -1


cdef class Histograms:
    cdef:
        int n_features
        uint32_t [::1] bin_offsets_view
        hist_struct [::1] histograms_view
    cdef public:
        object bin_offsets
        # Only the attribute histograms will be used in the splitter.
        object histograms

    cdef inline uint32_t n_bins(self, int feature_idx) noexcept nogil

    cdef inline hist_struct* at(self, int feature_idx, uint32_t bin_idx) noexcept nogil


cdef class Bitsets:
    cdef:
        const uint32_t [::1] offsets_view
        BITSET_INNER_DTYPE_C [::1] bitsets_view
    cdef public:
        object offsets
        # Only the attribute histograms will be used in the splitter.
        object bitsets

    cdef inline uint32_t n_inner_bitsets(self, int feature_idx) noexcept nogil

    cdef inline BITSET_INNER_DTYPE_C* at(self, int feature_idx) noexcept nogil


cdef class BinnedData:
    # A data container for mixed uint8 and uint16 columns.
    cdef:
        uint8_t[::1, :] X8_view
        uint16_t[::1, :] X16_view
        uint8_t[::1] feature_is_8bit_view  # Python bool ~ unsigned char = uint8
        int32_t[::1] feature_index_view
    cdef public:
        object X8
        object X16
        object feature_is_8bit
        object feature_index

    cdef inline uint8_t[::1] get_feature_view8(self, int feature_idx) noexcept nogil

    cdef inline uint16_t[::1] get_feature_view16(self, int feature_idx) noexcept nogil

    cdef inline uint8_t get_item8(self, int i, int feature_idx) noexcept nogil

    cdef inline uint16_t get_item16(self, int i, int feature_idx) noexcept nogil
