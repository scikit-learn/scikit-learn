# Author: Peter Prettenhofer
#
# License: BSD 3 clause

from libc.stdlib cimport free
from libc.string cimport memset

import numpy as np
cimport numpy as cnp
cnp.import_array()

from scipy.sparse import issparse

from ..tree._tree cimport Node
from ..tree._tree cimport Tree
from ..tree._tree cimport DTYPE_t
from ..tree._tree cimport SIZE_t
from ..tree._tree cimport INT32_t
from ..tree._utils cimport safe_realloc


# no namespace lookup for numpy dtype and array creation
from numpy import zeros as np_zeros


# constant to mark tree leafs
cdef SIZE_t TREE_LEAF = -1

cdef void _predict_regression_tree_inplace_fast_dense(
    const DTYPE_t[:, ::1] X,
    Node* root_node,
    double *value,
    double scale,
    Py_ssize_t k,
    Py_ssize_t K,
    Py_ssize_t n_samples,
    Py_ssize_t n_features,
    cnp.float64_t *out
) nogil:
    """Predicts output for regression tree and stores it in ``out[i, k]``.

    This function operates directly on the data arrays of the tree
    data structures. This is 5x faster than the variant above because
    it allows us to avoid buffer validation.

    The function assumes that the ndarray that wraps ``X`` is
    c-continuous.

    Parameters
    ----------
    X : DTYPE_t 2d memory view
        The memory view on the data ndarray of the input ``X``.
        Assumes that the array is c-continuous.
    root_node : tree Node pointer
        Pointer to the main node array of the :class:``sklearn.tree.Tree``.
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
    cdef Node *node
    for i in range(n_samples):
        node = root_node
        # While node not a leaf
        while node.left_child != TREE_LEAF:
            if X[i, node.feature] <= node.threshold:
                node = root_node + node.left_child
            else:
                node = root_node + node.right_child
        out[i * K + k] += scale * value[node - root_node]

def _predict_regression_tree_stages_sparse(
    object[:, :] estimators,
    object X,
    double scale,
    cnp.float64_t[:, :] out
):
    """Predicts output for regression tree inplace and adds scaled value to ``out[i, k]``.

    The function assumes that the ndarray that wraps ``X`` is csr_matrix.
    """
    cdef DTYPE_t[::1] X_data = X.data
    cdef INT32_t[::1] X_indices = X.indices
    cdef INT32_t[::1] X_indptr = X.indptr

    cdef SIZE_t n_samples = X.shape[0]
    cdef SIZE_t n_features = X.shape[1]
    cdef SIZE_t n_stages = estimators.shape[0]
    cdef SIZE_t n_outputs = estimators.shape[1]

    # Indices and temporary variables
    cdef SIZE_t sample_i
    cdef SIZE_t feature_i
    cdef SIZE_t stage_i
    cdef SIZE_t output_i
    cdef Node *root_node = NULL
    cdef Node *node = NULL
    cdef double *value = NULL

    cdef Tree tree
    cdef Node** nodes = NULL
    cdef double** values = NULL
    safe_realloc(&nodes, n_stages * n_outputs)
    safe_realloc(&values, n_stages * n_outputs)
    for stage_i in range(n_stages):
        for output_i in range(n_outputs):
            tree = estimators[stage_i, output_i].tree_
            nodes[stage_i * n_outputs + output_i] = tree.nodes
            values[stage_i * n_outputs + output_i] = tree.value

    # Initialize auxiliary data-structure
    cdef DTYPE_t feature_value = 0.
    cdef DTYPE_t* X_sample = NULL

    # feature_to_sample as a data structure records the last seen sample
    # for each feature; functionally, it is an efficient way to identify
    # which features are nonzero in the present sample.
    cdef SIZE_t* feature_to_sample = NULL

    safe_realloc(&X_sample, n_features)
    safe_realloc(&feature_to_sample, n_features)

    memset(feature_to_sample, -1, n_features * sizeof(SIZE_t))

    # Cycle through all samples
    for sample_i in range(n_samples):
        for feature_i in range(X_indptr[sample_i], X_indptr[sample_i + 1]):
            feature_to_sample[X_indices[feature_i]] = sample_i
            X_sample[X_indices[feature_i]] = X_data[feature_i]

        # Cycle through all stages
        for stage_i in range(n_stages):
            # Cycle through all trees
            for output_i in range(n_outputs):
                root_node = nodes[stage_i * n_outputs + output_i]
                value = values[stage_i * n_outputs + output_i]
                node = root_node

                # While node not a leaf
                while node.left_child != TREE_LEAF:
                    # ... and node.right_child != TREE_LEAF:
                    if feature_to_sample[node.feature] == sample_i:
                        feature_value = X_sample[node.feature]
                    else:
                        feature_value = 0.

                    if feature_value <= node.threshold:
                        node = root_node + node.left_child
                    else:
                        node = root_node + node.right_child
                out[sample_i, output_i] += scale * value[node - root_node]

    # Free auxiliary arrays
    free(X_sample)
    free(feature_to_sample)
    free(nodes)
    free(values)


def predict_stages(
    object[:, :] estimators,
    object X,
    double scale,
    cnp.float64_t[:, :] out
):
    """Add predictions of ``estimators`` to ``out``.

    Each estimator is scaled by ``scale`` before its prediction
    is added to ``out``.
    """
    cdef Py_ssize_t i
    cdef Py_ssize_t k
    cdef Py_ssize_t n_estimators = estimators.shape[0]
    cdef Py_ssize_t K = estimators.shape[1]
    cdef Tree tree

    if issparse(X):
        if X.format != 'csr':
            raise ValueError("When X is a sparse matrix, a CSR format is"
                             " expected, got {!r}".format(type(X)))
        _predict_regression_tree_stages_sparse(
            estimators=estimators, X=X, scale=scale, out=out
        )
    else:
        if not isinstance(X, np.ndarray) or np.isfortran(X):
            raise ValueError(f"X should be C-ordered np.ndarray, got {type(X)}")

        for i in range(n_estimators):
            for k in range(K):
                tree = estimators[i, k].tree_

                # avoid buffer validation by casting to ndarray
                # and get data pointer
                # need brackets because of casting operator priority
                _predict_regression_tree_inplace_fast_dense(
                    X=X,
                    root_node=tree.nodes,
                    value=tree.value,
                    scale=scale,
                    k=k,
                    K=K,
                    n_samples=X.shape[0],
                    n_features=X.shape[1],
                    out=&out[0, 0]
                )
                ## out[:, k] += scale * tree.predict(X).ravel()


def predict_stage(
    object[:, :] estimators,
    int stage,
    object X,
    double scale,
    cnp.float64_t[:, :] out
):
    """Add predictions of ``estimators[stage]`` to ``out``.

    Each estimator in the stage is scaled by ``scale`` before
    its prediction is added to ``out``.
    """
    return predict_stages(
        estimators=estimators[stage:stage + 1], X=X, scale=scale, out=out
    )


def _random_sample_mask(
    cnp.npy_intp n_total_samples,
    cnp.npy_intp n_total_in_bag,
    random_state
):
     """Create a random sample mask where ``n_total_in_bag`` elements are set.

     Parameters
     ----------
     n_total_samples : int
         The length of the resulting mask.

     n_total_in_bag : int
         The number of elements in the sample mask which are set to 1.

     random_state : RandomState
         A numpy ``RandomState`` object.

     Returns
     -------
     sample_mask : np.ndarray, shape=[n_total_samples]
         An ndarray where ``n_total_in_bag`` elements are set to ``True``
         the others are ``False``.
     """
     cdef cnp.float64_t[::1] rand = random_state.uniform(size=n_total_samples)
     cdef cnp.uint8_t[::1] sample_mask = np_zeros((n_total_samples,), dtype=bool)

     cdef cnp.npy_intp n_bagged = 0
     cdef cnp.npy_intp i = 0

     for i in range(n_total_samples):
         if rand[i] * (n_total_samples - i) < (n_total_in_bag - n_bagged):
             sample_mask[i] = 1
             n_bagged += 1

     return sample_mask.base
