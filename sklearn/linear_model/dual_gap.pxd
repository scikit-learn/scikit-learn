from types cimport floating, complexing

cdef floating compute_dual_gap(int n_samples, int n_features, int n_tasks,
                               complexing *W_ptr,  # C-order 2d array
                               floating reg, floating l2_reg,
                               complexing *X_or_Gram_conj_ptr,  # F-order 2d
                               complexing *Y_or_Cov_ptr,  # F-order 2d
                               complexing *R_or_RCov_ptr,  # F-order 2d
                               complexing *XtA_ptr,  # F-order 2d
                               floating Y_norm2, bint precomputed,
                               int penalty_model, bint positive) nogil except *
