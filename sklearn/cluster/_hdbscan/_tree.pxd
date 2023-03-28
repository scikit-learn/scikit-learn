cimport numpy as cnp

# This corresponds to the scipy.cluster.hierarchy format
ctypedef packed struct HIERARCHY_t:
    cnp.intp_t left_node
    cnp.intp_t right_node
    cnp.float64_t value
    cnp.intp_t cluster_size

# Effectively an edgelist encoding a parent/child pair, along with a value and
# the corresponding cluster_size in each row providing a tree structure.
ctypedef packed struct CONDENSED_t:
    cnp.intp_t parent
    cnp.intp_t child
    cnp.float64_t value
    cnp.intp_t cluster_size
