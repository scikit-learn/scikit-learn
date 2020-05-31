# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# cython: nonecheck=False
# distutils: language=c++
from ._bitset cimport in_bitset
from .common cimport BITSET_INNER_DTYPE_C
from libc.limits cimport CHAR_BIT


cdef class PredictorBitSet:
    def insert_categories_bitset(self, unsigned int node_idx,
                                 X_DTYPE_C[:] category_bins,
                                 BITSET_INNER_DTYPE_C[:] cat_bitset):
        # get cateogries from cat_bitset
        cdef:
            BITSET_INNER_DTYPE_C val
            int k, offset, cat_bin
            int cardinality = category_bins.shape[0]
            int BITSET_SIZE = sizeof(BITSET_INNER_DTYPE_C) * CHAR_BIT

        self.node_to_binned_bitset[node_idx].resize(cat_bitset.shape[0])

        for k in range(cat_bitset.shape[0]):
            offset = BITSET_SIZE * k
            val = cat_bitset[k]
            self.node_to_binned_bitset[node_idx][k] = val
            while val and offset < cardinality:
                if val % 2:
                    cat_bin = <int>(category_bins[offset])
                    self.node_to_raw_bitset[node_idx].insert(cat_bin)
                val = val // 2
                offset += 1

    cdef unsigned char raw_category_in_bitset(self, unsigned int node_idx,
                                              X_DTYPE_C category) nogil:
        return self.node_to_raw_bitset[node_idx].count(<int>category)

    cdef unsigned char binned_category_in_bitset(self, unsigned int node_idx,
                                                 X_BINNED_DTYPE_C category) nogil:
        cdef:
            unsigned int i1 = category // 32
            unsigned int i2 = category % 32
            vector[BITSET_INNER_DTYPE_C] bitset = \
                self.node_to_binned_bitset[node_idx]
        return (bitset[i1] >> i2) & 1

    def get_raw_categories(self, unsigned int node_idx):
        """Used for testing"""
        return self.node_to_raw_bitset[node_idx]

    def get_binned_categories(self, unsigned int node_idx):
        """Used for testing"""
        return self.node_to_binned_bitset[node_idx]
