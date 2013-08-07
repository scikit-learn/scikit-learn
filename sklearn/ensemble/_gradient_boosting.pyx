# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer
#
# Licence: BSD 3 clause

cimport cython

import numpy as np
cimport numpy as np
np.import_array()

from sklearn.tree._tree cimport Tree

ctypedef np.int32_t int32
ctypedef np.float64_t float64
ctypedef np.int8_t int8

from numpy import bool as np_bool

# no namespace lookup for numpy dtype and array creation
from numpy import zeros as np_zeros
from numpy import ones as np_ones
from numpy import bool as np_bool
from numpy import int8 as np_int8
from numpy import intp as np_intp
from numpy import float32 as np_float32
from numpy import float64 as np_float64

# Define a datatype for the data array
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t
ctypedef np.npy_intp SIZE_t


# constant to mark tree leafs
cdef int LEAF = -1

cdef void _predict_regression_tree_inplace_fast(DTYPE_t *X,
                                                SIZE_t *children_left,
                                                SIZE_t *children_right,
                                                SIZE_t *feature,
                                                double *threshold,
                                                double *value,
                                                double scale,
                                                Py_ssize_t k,
                                                Py_ssize_t K,
                                                Py_ssize_t n_samples,
                                                Py_ssize_t n_features,
                                                float64 *out):
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
    children : np.npy_intp pointer
        The pointer to the data array of the ``children`` array attribute
        of the :class:``sklearn.tree.Tree``.
    feature : np.npy_intp pointer
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
    cdef int32 node_id
    cdef SIZE_t feature_idx
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
                   np.ndarray[float64, ndim=2] out):
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
            tree = estimators[i, k].tree_

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
                <float64*>((<np.ndarray>out).data))
            ## out += scale * tree.predict(X).reshape((X.shape[0], 1))


@cython.nonecheck(False)
def predict_stage(np.ndarray[object, ndim=2] estimators,
                  int stage,
                  np.ndarray[DTYPE_t, ndim=2] X, double scale,
                  np.ndarray[float64, ndim=2] out):
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
        tree = estimators[stage, k].tree_

        _predict_regression_tree_inplace_fast(
                <DTYPE_t*>(X.data),
                tree.children_left,
                tree.children_right,
                tree.feature,
                tree.threshold,
                tree.value,
                scale, k, K, n_samples, n_features,
                <float64*>((<np.ndarray>out).data))
        ## out += scale * tree.predict(X).reshape((X.shape[0], 1))


cdef inline int array_index(int32 val, int32[::1] arr):
    """Find index of ``val`` in array ``arr``. """
    cdef int32 res = -1
    cdef int32 i = 0
    cdef int32 n = arr.shape[0]
    for i in range(n):
        if arr[i] == val:
            res = i
            break
    return res


cpdef _partial_dependence_tree(Tree tree, DTYPE_t[:, ::1] X,
                               int32[::1] target_feature,
                               double learn_rate,
                               double[::1] out):
    """Partial dependence of the response on the ``target_feature`` set.

    For each row in ``X`` a tree traversal is performed.
    Each traversal starts from the root with weight 1.0.

    At each non-terminal node that splits on a target variable either
    the left child or the right child is visited based on the feature
    value of the current sample and the weight is not modified.
    At each non-terminal node that splits on a complementary feature
    both children are visited and the weight is multiplied by the fraction
    of training samples which went to each child.

    At each terminal node the value of the node is multiplied by the
    current weight (weights sum to 1 for all visited terminal nodes).

    Parameters
    ----------
    tree : sklearn.tree.Tree
        A regression tree; tree.values.shape[1] == 1
    X : memory view on 2d ndarray
        The grid points on which the partial dependence
        should be evaluated. X.shape[1] == target_feature.shape[0].
    target_feature : memory view on 1d ndarray
        The set of target features for which the partial dependence
        should be evaluated. X.shape[1] == target_feature.shape[0].
    learn_rate : double
        Constant scaling factor for the leaf predictions.
    out : memory view on 1d ndarray
        The value of the partial dependence function on each grid
        point.
    """
    cdef Py_ssize_t i = 0
    cdef Py_ssize_t n_features = X.shape[1]
    cdef SIZE_t *children_left = tree.children_left
    cdef SIZE_t *children_right = tree.children_right
    cdef SIZE_t *feature = tree.feature
    cdef double *value = tree.value
    cdef double *threshold = tree.threshold
    cdef SIZE_t *n_node_samples = tree.n_node_samples
    cdef SIZE_t node_count = tree.node_count

    cdef SIZE_t stack_capacity = node_count * 2
    cdef SIZE_t[::1] node_stack = np_zeros((stack_capacity,), dtype=np_intp)
    cdef double[::1] weight_stack = np_ones((stack_capacity,), dtype=np_float64)
    cdef SIZE_t stack_size = 1
    cdef double left_sample_frac
    cdef double current_weight
    cdef double total_weight = 0.0
    cdef SIZE_t current_node

    for i in range(X.shape[0]):
        # init stacks for new example
        stack_size = 1
        node_stack[0] = 0
        weight_stack[0] = 1.0
        total_weight = 0.0

        while stack_size > 0:
            # get top node on stack
            stack_size -= 1
            current_node = node_stack[stack_size]

            if children_left[current_node] == LEAF:
                out[i] += weight_stack[stack_size] * value[current_node] * \
                          learn_rate
                total_weight += weight_stack[stack_size]
            else:
                # non-terminal node
                feature_index = array_index(feature[current_node], target_feature)
                if feature_index != -1:
                    # split feature in target set
                    # push left or right child on stack
                    if X[i, feature_index] <= threshold[current_node]:
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
                    left_sample_frac = n_node_samples[children_left[current_node]] / \
                                       <double>n_node_samples[current_node]
                    if left_sample_frac <= 0.0 or left_sample_frac >= 1.0:
                        raise ValueError("left_sample_frac:%f, "
                                         "n_samples current: %d, "
                                         "n_samples left: %d"
                                         % (left_sample_frac,
                                            n_node_samples[current_node],
                                            n_node_samples[children_left[current_node]]))
                    weight_stack[stack_size] = current_weight * left_sample_frac
                    stack_size +=1

                    # push right child
                    node_stack[stack_size] = children_right[current_node]
                    weight_stack[stack_size] = current_weight * \
                                               (1.0 - left_sample_frac)
                    stack_size +=1

        if not (0.999 < total_weight < 1.001):
            raise ValueError("Total weight should be 1.0 but was %.9f" %
                             total_weight)


def _random_sample_mask(int n_total_samples, int n_total_in_bag, random_state):
     """Create a random sample mask where ``n_total_in_bag`` elements are set.

     Parameters
     ----------
     n_total_samples : int
         The length of the resulting mask.

     n_total_in_bag : int
         The number of elements in the sample mask which are set to 1.

     random_state : np.RandomState
         A numpy ``RandomState`` object.

     Returns
     -------
     sample_mask : np.ndarray, shape=[n_total_samples]
         An ndarray where ``n_total_in_bag`` elements are set to ``True``
         the others are ``False``.
     """
     cdef np.ndarray[float64, ndim=1, mode="c"] rand = \
          random_state.rand(n_total_samples)
     cdef np.ndarray[int8, ndim=1, mode="c"] sample_mask = \
          np_zeros((n_total_samples,), dtype=np_int8)

     cdef int n_bagged = 0
     cdef int i = 0

     for i from 0 <= i < n_total_samples:
         if rand[i] * (n_total_samples - i) < (n_total_in_bag - n_bagged):
             sample_mask[i] = 1
             n_bagged += 1

     return sample_mask.astype(np_bool)
