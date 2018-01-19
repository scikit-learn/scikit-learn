from cython cimport floating

cdef floating compute_dual_gap(int n_samples, int n_features, int n_tasks,
                               floating *W_ptr,  # C-order 2d array
                               floating reg, floating l2_reg,
                               floating *X_or_Gram_ptr,  # F-order 2d
                               floating *Y_or_Cov_ptr,  # F-order 2d
                               floating *R_or_RCov_ptr,  # F-order 2d
                               floating *XtA_ptr,  # F-order 2d
                               floating Y_norm2, bint precomputed,
                               int penalty_model, bint positive) nogil
