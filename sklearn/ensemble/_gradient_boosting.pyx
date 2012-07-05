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


cdef void _predict_regression_tree_inplace_fast(DTYPE_t *X,
                                                np.int32_t *children,
                                                np.int32_t *feature,
                                                np.float64_t *threshold,
                                                np.float64_t * value,
                                                double scale,
                                                Py_ssize_t k,
                                                Py_ssize_t K,
                                                Py_ssize_t n_samples,
                                                Py_ssize_t n_features,
                                                np.float64_t *out):
    """Predicts output for regression tree and stores it in ``out[i, k]``.

    This function operates directly on the data arrays of the tree
    data structures. This is 5x faster than the variant above because
    it allows us to avoid buffer validation.

    Parameters
    ----------
    X : DTYPE_t pointer
        The pointer to the data array of the input ``X``.
    children : np.int32_t pointer
        The pointer to the data array of the ``children`` array attribute
        of the :class:``sklearn.tree.Tree``.
    feature : np.int32_t pointer
        The pointer to the data array of the ``feature`` array attribute
        of the :class:``sklearn.tree.Tree``.
    threshold : np.float64_t pointer
        The pointer to the data array of the ``threshold`` array attribute
        of the :class:``sklearn.tree.Tree``.
    value : np.float64_t pointer
        The pointer to the data array of the ``value`` array attribute
        of the :class:``sklearn.tree.Tree``.
    scale : double
        A constant to scale the predictions.
    k : int
        The index of the tree output to be predicted. Must satisfy
        0 <= ``k`` < ``K``.
    K : int
        The number of regression tree outputs. For regression and
        binary classification ``K == 1``, for multi-class
        classification ``K == n_classes``.
    n_samples : int
        The number of samples in the input array ``X``;
        ``n_samples == X.shape[0]``.
    n_features : int
        The number of features; ``n_samples == X.shape[1]``.
    out : np.float64_t pointer
        The pointer to the data array where the predictions are stored.
        ``out`` is assumed to be a two-dimensional array of
        shape ``(n_samples, K)``.
    """
    cdef Py_ssize_t i
    cdef np.int32_t node_id
    cdef np.int32_t feature_idx
    cdef int stride = 2  # children.shape[1]
    for i in range(n_samples):
        node_id = 0
        # While node_id not a leaf
        while children[node_id * stride] != -1 and \
                  children[(node_id * stride) + 1] != -1:
            feature_idx = feature[node_id]
            if X[(i * n_features) + feature_idx] <= threshold[node_id]:
                node_id = children[node_id * stride]
            else:
                node_id = children[(node_id * stride) + 1]
        out[(i * K) + k] += scale * value[node_id]


@cython.nonecheck(False)
def predict_stages(np.ndarray[object, ndim=2] estimators,
                   np.ndarray[DTYPE_t, ndim=2, mode='c'] X, double scale,
                   np.ndarray[np.float64_t, ndim=2] out):
    """Add predictions of ``estimators`` to ``out``.

    Each estimator is scaled by ``scale`` before its prediction
    is added to ``out``.
    """
    cdef Py_ssize_t i
    cdef Py_ssize_t k
    cdef Py_ssize_t n_estimators = estimators.shape[0]
    cdef Py_ssize_t n_samples = X.shape[0]
    cdef Py_ssize_t n_features = X.shape[1]
    cdef Py_ssize_t K = estimators.shape[1]
    cdef object tree

    for i in range(n_estimators):
        for k in range(K):
            tree = estimators[i, k]

            # avoid buffer validation by casting to ndarray
            # and get data pointer
            # need brackets because of casting operator priority
            _predict_regression_tree_inplace_fast(
                <DTYPE_t*>(X.data),
                <np.int32_t*>((<np.ndarray>(tree.children)).data),
                <np.int32_t*>((<np.ndarray>(tree.feature)).data),
                <np.float64_t*>((<np.ndarray>(tree.threshold)).data),
                <np.float64_t*>((<np.ndarray>(tree.value)).data),
                scale, k, K, n_samples, n_features,
                <np.float64_t*>((<np.ndarray>out).data))


@cython.nonecheck(False)
def predict_stage(np.ndarray[object, ndim=2] estimators,
                  int stage,
                  np.ndarray[DTYPE_t, ndim=2] X, double scale,
                  np.ndarray[np.float64_t, ndim=2] out):
    """Add predictions of ``estimators[stage]`` to ``out``.

    Each estimator in the stage is scaled by ``scale`` before
    its prediction is added to ``out``.
    """
    cdef Py_ssize_t i
    cdef Py_ssize_t k
    cdef Py_ssize_t n_estimators = estimators.shape[0]
    cdef Py_ssize_t n_samples = X.shape[0]
    cdef Py_ssize_t n_features = X.shape[1]
    cdef Py_ssize_t K = estimators.shape[1]
    cdef object tree
    for k in range(K):
        tree = estimators[stage, k]

        _predict_regression_tree_inplace_fast(
                <DTYPE_t*>(X.data),
                <np.int32_t*>((<np.ndarray>(tree.children)).data),
                <np.int32_t*>((<np.ndarray>(tree.feature)).data),
                <np.float64_t*>((<np.ndarray>(tree.threshold)).data),
                <np.float64_t*>((<np.ndarray>(tree.value)).data),
                scale, k, K, n_samples, n_features,
                <np.float64_t*>((<np.ndarray>out).data))

