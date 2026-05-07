from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, uint8_t


ctypedef float32_t X_DTYPE_C  # C version of X_DTYPE
ctypedef float64_t Y_DTYPE_C  # C version of Y_DTYPE

cdef struct node_struct:
    # Base storage structure for the nodes in a Tree object
    intp_t left_child                    # id of the left child of the node
    intp_t right_child                   # id of the right child of the node
    intp_t feature                       # Feature used for splitting the node
    float64_t threshold                  # Threshold value at the node
    Y_DTYPE_C impurity                   # Impurity of the node (i.e., the value of the criterion)
    intp_t n_node_samples                # Number of samples at the node
    Y_DTYPE_C weighted_n_node_samples    # Weighted number of samples at the node
    uint8_t missing_go_to_left           # Whether features have missing values
