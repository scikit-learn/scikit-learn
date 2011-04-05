"""
Binding for libdecisiontree
---------------------------

Notes
-----

Authors
-------
2011: Noel Dawe <Noel.Dawe@cern.ch>
"""

import  numpy as np
cimport numpy as np

################################################################################
# Includes

cdef extern from "Node.h":
    cdef cppclass Node:
        void recursive_split(int, int)

cdef extern from "libdecisiontree_helper.cpp":
    void copy_predict (char *, Node *, np.npy_intp *, char *)
    Node* create_root(char *, char *, char *, np.npy_intp *)

################################################################################
# Wrapper functions

def fit(np.ndarray[np.float64_t, ndim=2, mode='c'] X,
        np.ndarray[np.float64_t, ndim=1, mode='c'] Y,
        np.ndarray[np.float64_t, ndim=1, mode='c'] sample_weight = np.empty(0),
        int minleafsize = 10,
        int nbins = 0):

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
    cdef Node* root

    if len(sample_weight) == 0:
        sample_weight = np.ones(X.shape[0], dtype=np.float64)
    else:
        assert sample_weight.shape[0] == X.shape[0], \
               "sample_weight and X have incompatible shapes: " + \
               "sample_weight has %s samples while X has %s" % \
               (sample_weight.shape[0], X.shape[0])

    root = create_root(X.data, Y.data, sample_weight.data, X.shape)
    root.recursive_split(minleafsize, nbins)
    return root

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
    
    dec_values = np.empty(X.shape[0])
    copy_predict(X.data, root, X.shape, dec_values.data)
    return dec_values
