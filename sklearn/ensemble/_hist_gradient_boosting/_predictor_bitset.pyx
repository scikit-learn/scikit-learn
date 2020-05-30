# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# cython: nonecheck=False
# distutils: language=c++
from ._bitset cimport in_bitset


cdef class PredictorBitSet:
    def insert_categories_bitset(self, unsigned int node_idx,
                                 floating[:] category_bins,
                                 BITSET_INNER_DTYPE_C[:] cat_bitset):
        # get cateogries from cat_bitset
        cdef:
            BITSET_INNER_DTYPE_C val
            unsigned int k, bit, offset
            int cat_bin

        for k in range(cat_bitset.shape[0]):
            # BITSET_INNER_DTYPE_C is bit 32
            offset = 32 * k
            val = cat_bitset[k]
            while val:
                bit = val % 2
                if bit:
                    cat_bin = <int>(category_bins[bit]) + offset
                    self.node_to_raw_bitset[node_idx].insert(cat_bin)
                val = val // 2
        self.node_to_binned_bitset[node_idx] = cat_bitset[0]

    cdef unsigned char raw_category_in_bitset(self, unsigned int node_idx, floating category) nogil:
        if self.node_to_raw_bitset.count(node_idx) == 0:
            return 0
        return self.node_to_raw_bitset[node_idx].count(<int>category)

    cdef unsigned char binned_category_in_bitset(self, unsigned int node_idx, X_BINNED_DTYPE_C category) nogil:
        if self.node_to_binned_bitset.count(node_idx) == 0:
            return 0
        return in_bitset(category, &self.node_to_binned_bitset[node_idx])
