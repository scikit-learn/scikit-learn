"""
Binding for libdecisiontree
---------------------------

Notes
-----

Authors
-------
2011: Noel Dawe <end1@sfu.ca>
"""

import  numpy as np
cimport numpy as np

################################################################################
# Includes

cdef extern from "Object.h":
    cdef struct Object

cdef extern from "Node.h":
    cdef cppclass Node:

################################################################################
# Wrapper functions

def train(np.ndarray[np.float64_t, ndim=2, mode='c'] X,
          np.ndarray[np.float64_t, ndim=1, mode='c'] Y,
          np.ndarray[np.float64_t, ndim=1, mode='c'] sample_weight=np.empty(0),
          int min_leaf_size = 10,
          int nbins=0):

    """
    Parameters
    ----------
    X: array-like, dtype=float, size=[n_samples, n_features]

    Y: array, dtype=float, size=[n_samples]
        target vector

    Return
    ------
    root: root Node of the trained decision tree

    classification : resulting classification of the training set
    """

    cdef np.npy_intp nr

    if len(sample_weight) == 0:
        sample_weight = np.ones(X.shape[0], dtype=np.float64)
    else:
        assert sample_weight.shape[0] == X.shape[0], \
               "sample_weight and X have incompatible shapes: " + \
               "sample_weight has %s samples while X has %s" % \
               (sample_weight.shape[0], X.shape[0])

    # set libsvm problem
    problem = set_problem(X.data, Y.data, sample_weight.data, X.shape)

    # check parameters
    if problem == NULL:
        raise MemoryError("Seems we've run out of of memory")

    # this does the real work
    model = svm_train(problem, param)

    # from here until the end, we just copy the data returned by
    # svm_train
    SV_len  = get_l(model)
    n_class = get_nr(model)

    # copy model.sv_coef
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] sv_coef
    sv_coef = np.empty((n_class-1, SV_len), dtype=np.float64)
    copy_sv_coef (sv_coef.data, model)

    # copy model.rho into the intercept
    # the intercept is just model.rho but with sign changed
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] intercept
    intercept = np.empty(n_class*(n_class-1)/2, dtype=np.float64)
    copy_intercept (intercept.data, model, intercept.shape)

    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] support
    support = np.empty (SV_len, dtype=np.int32)
    copy_support (support.data, model)

    # copy model.SV
    cdef np.ndarray[np.float64_t, ndim=2, mode='c'] support_vectors
    if kernel_type == 4:
        support_vectors = np.empty((0, 0), dtype=np.float64)
    else:
        support_vectors = np.empty((SV_len, X.shape[1]), dtype=np.float64)
        copy_SV(support_vectors.data, model, support_vectors.shape)

    # copy model.nSV
    # TODO: do only in classification
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] n_class_SV
    n_class_SV = np.empty(n_class, dtype=np.int32)
    copy_nSV(n_class_SV.data, model)

    # copy label
    cdef np.ndarray[np.int32_t, ndim=1, mode='c'] label
    label = np.empty((n_class), dtype=np.int32)
    copy_label(label.data, model)

    # copy probabilities
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] probA
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] probB
    if probability != 0:
        if svm_type < 2: # SVC and NuSVC
            probA = np.empty(n_class*(n_class-1)/2, dtype=np.float64)
            probB = np.empty(n_class*(n_class-1)/2, dtype=np.float64)
            copy_probB(probB.data, model, probB.shape)
        else:
            probA = np.empty(1, dtype=np.float64)
            probB = np.empty(0, dtype=np.float64)
        copy_probA(probA.data, model, probA.shape)

    # memory deallocation
    svm_free_and_destroy_model(&model)
    free_problem(problem)
    free_param(param)

    return root, classification


def predict(np.ndarray[np.float64_t, ndim=2, mode='c'] X,
            Node root):
    """
    Predict target values of X given a model (low-level method)

    Parameters
    ----------
    X: array-like, dtype=float, size=[n_samples, n_features]

    root: decision tree root Node

    Return
    ------
    dec_values : array
        predicted values.
    """
    cdef np.ndarray[np.float64_t, ndim=1, mode='c'] dec_values
    
    #TODO: use check_model
    dec_values = np.empty(X.shape[0])
    if copy_predict(X.data, root, X.shape, dec_values.data) < 0:
        raise MemoryError("We've run out of of memory")
    return dec_values
