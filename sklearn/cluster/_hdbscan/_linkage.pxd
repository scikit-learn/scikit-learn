from ...utils._typedefs cimport float64_t, int64_t

# Packed shouldn't make a difference since they're all 8-byte quantities,
# but it's included just to be safe.
ctypedef packed struct MST_edge_t:
    int64_t current_node
    int64_t next_node
    float64_t distance
