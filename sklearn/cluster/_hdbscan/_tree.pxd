cimport numpy as cnp
import numpy as np

# This corresponds to the scipy.cluster.hierarchy format
ctypedef packed struct HIERARCHY_t:
    cnp.intp_t left_node
    cnp.intp_t right_node
    cnp.float64_t value
    cnp.intp_t cluster_size

ctypedef packed struct CONDENSED_t:
    cnp.intp_t parent
    cnp.intp_t child
    cnp.float64_t value
    cnp.intp_t cluster_size
