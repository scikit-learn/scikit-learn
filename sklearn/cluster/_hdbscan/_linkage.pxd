cimport numpy as cnp
import numpy as np

# Numpy structured dtype representing a single ordered edge in Prim's algorithm
MST_edge_dtype = np.dtype([
    ("current_node", np.intp),
    ("next_node", np.intp),
    ("distance", np.float64),
])

ctypedef struct MST_edge_t:
    cnp.intp_t current_node
    cnp.intp_t next_node
    cnp.float64_t distance
