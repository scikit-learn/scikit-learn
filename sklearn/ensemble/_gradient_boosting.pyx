# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer
#
# License: BSD Style.

cimport cython

import numpy as np
cimport numpy as np

# Define a datatype for the data array
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t


cdef void _predict_regression_tree_inplace(np.ndarray[DTYPE_t, ndim=2] X,
                                           np.ndarray[np.int32_t, ndim=2] children,
                                           np.ndarray[np.int32_t, ndim=1] feature,
                                           np.ndarray[np.float64_t, ndim=1] threshold,
                                           np.ndarray[np.float64_t, ndim=2] values,
                                           double scale,
                                           int k,
                                           np.ndarray[np.float64_t, ndim=2] pred):
    """Predicts output for regression tree and stores it in ``pred[i, k]`` """
    cdef int i = 0
    cdef int n = X.shape[0]
    cdef int node_id = 0
    cdef int K = values.shape[1]
    for i in range(n):
        node_id = 0
        # While node_id not a leaf
        while children[node_id, 0] != -1 and children[node_id, 1] != -1:
            if X[i, feature[node_id]] <= threshold[node_id]:
                node_id = children[node_id, 0]
            else:
                node_id = children[node_id, 1]
        pred[i, k] += scale * values[node_id, 0]


@cython.nonecheck(False)
def predict_stages(np.ndarray[object, ndim=2] estimators,
                   np.ndarray[DTYPE_t, ndim=2] X, double scale,
                   np.ndarray[np.float64_t, ndim=2] out):
    """Add predictions of ``estimators`` to ``out``.

    Each estimator is scaled by ``scale`` before its prediction
    is added to ``out``.
    """
    cdef int i
    cdef int k
    cdef int n_estimators = estimators.shape[0]
    cdef int K = estimators.shape[1]
    for i in range(n_estimators):
        for k in range(K):
            tree = estimators[i, k]
            _predict_regression_tree_inplace(X, tree.children, tree.feature,
                                             tree.threshold, tree.value,
                                             scale, k, out)


@cython.nonecheck(False)
def predict_stage(np.ndarray[object, ndim=2] estimators,
                  int stage,
                  np.ndarray[DTYPE_t, ndim=2] X, double scale,
                  np.ndarray[np.float64_t, ndim=2] pred):
    """Add predictions of ``estimators[stage]`` to ``out``.

    Each estimator in the stage is scaled by ``scale`` before
    its prediction is added to ``out``.
    """
    cdef int i
    cdef int k
    cdef int n_estimators = estimators.shape[0]
    cdef int K = estimators.shape[1]
    for k in range(K):
        tree = estimators[stage, k]
        _predict_regression_tree_inplace(X, tree.children, tree.feature,
                                         tree.threshold, tree.value,
                                         scale, k, pred)
