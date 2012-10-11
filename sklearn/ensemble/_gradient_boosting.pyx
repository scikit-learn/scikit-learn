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

from sklearn.tree._tree cimport Tree

# Define a datatype for the data array
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

# constant to mark tree leafs
cdef int LEAF = -1

cdef void _predict_regression_tree_inplace_fast(DTYPE_t *X,
                                                int *children_left,
                                                int *children_right,
                                                int *feature,
                                                double *threshold,
                                                double *value,
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

    The function assumes that the ndarray that wraps ``X`` is
    c-continuous.

    Parameters
    ----------
    X : DTYPE_t pointer
        The pointer to the data array of the input ``X``.
        Assumes that the array is c-continuous.
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
    for i in range(n_samples):
        node_id = 0
        # While node_id not a leaf
        while children_left[node_id] != -1 and \
                  children_right[node_id] != -1:
            feature_idx = feature[node_id]
            if X[(i * n_features) + feature_idx] <= threshold[node_id]:
                node_id = children_left[node_id]
            else:
                node_id = children_right[node_id]
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
    cdef Tree tree

    for i in range(n_estimators):
        for k in range(K):
            tree = estimators[i, k]

            # avoid buffer validation by casting to ndarray
            # and get data pointer
            # need brackets because of casting operator priority
            _predict_regression_tree_inplace_fast(
                <DTYPE_t*>(X.data),
                tree.children_left,
                tree.children_right,
                tree.feature,
                tree.threshold,
                tree.value,
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
    cdef Tree tree
    for k in range(K):
        tree = estimators[stage, k]

        _predict_regression_tree_inplace_fast(
                <DTYPE_t*>(X.data),
                tree.children_left,
                tree.children_right,
                tree.feature,
                tree.threshold,
                tree.value,
                scale, k, K, n_samples, n_features,
                <np.float64_t*>((<np.ndarray>out).data))


cdef inline int array_index(int val, int[::1] arr):
    """Find index of ``val`` in array ``arr``. """
    cdef int res = -1
    cdef int i = 0
    cdef int n = arr.shape[0]
    for i in range(n):
        if arr[i] == val:
            res = i
            break
    return res


cdef void _partial_dependency_tree(Tree tree, DTYPE_t[:, ::1] X,
                              int[::1] target_feature,
                              double[::1] out):
    """Partial dependency of ``target_feature`` set on the response.

    For each row in ``X`` a tree traversal is performed.
    Each traversal starts from the root with weight 1.0.
    At each non-terminal node that splits on a target variable either
    the left child or the right child is visited based on the feature
    value of the current sample and the weight is not modified.
    At each non-terminal node that splits on a complementary feature
    both children are visited and the weight is multiplied by the fraction
    of training samples which went to each child.
    At each terminal node the value of the node
    """
    cdef Py_ssize_t i = 0
    cdef Py_ssize_t n_features = X.shape[1]
    cdef int *children_left = tree.children_left
    cdef int *children_right = tree.children_right
    cdef int* feature = tree.feature
    cdef double* value = tree.value
    cdef double* threshold = tree.threshold
    cdef int* n_samples = tree.n_samples
    cdef int node_count = tree.node_count

    cdef int[::1] node_stack = np.zeros((node_count,), dtype=np.int)
    cdef double[::1] weight_stack = np.ones((node_count,), dtype=np.float64)
    cdef int stack_size = 1
    cdef double left_sample_frac
    cdef double current_weight

    for i in range(X.shape[0]):

        #for j in range(node_count):
        while stack_size > 0:
            stack_size -= 1

            # get top node on stack
            current_node = node_stack[stack_size]

            if children_left[current_node] == LEAF:
                out[i] += weight_stack[stack_size] * value[current_node]
            else:
                # non-terminal node
                feature_index = array_index(feature[current_node], target_feature)
                if feature_index != -1:
                    # split feature in target set
                    # push left or right child on stack
                    if X[i, feature_index] < threshold[current_node]:
                        # left
                        node_stack[stack_size] = children_left[current_node]
                    else:
                        # right
                        node_stack[stack_size] = children_right[current_node]
                    stack_size += 1
                else:
                    # split feature in complement set
                    # push both children onto stack

                    # push left child
                    node_stack[stack_size] = children_left[current_node]
                    current_weight = weight_stack[stack_size]
                    left_sample_frac = n_samples[children_left[current_node]] / n_samples[current_node]
                    assert 0.0 < left_sample_frac < 1.0
                    weight_stack[stack_size] = current_weight * left_sample_frac
                    stack_size +=1

                    # push right child
                    node_stack[stack_size] = children_right[current_node]
                    weight_stack[stack_size] = current_weight * (1.0 - left_sample_frac)
                    stack_size +=1
