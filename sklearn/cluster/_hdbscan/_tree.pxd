cimport numpy as cnp
import numpy as np


ctypedef packed struct HIERARCHY_t:
    cnp.intp_t left_node
    cnp.intp_t right_node
    cnp.float64_t value
    cnp.intp_t cluster_size
