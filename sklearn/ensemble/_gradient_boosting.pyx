# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.stdlib cimport free
from libc.string cimport memset

import numpy as np
from scipy.sparse import issparse

from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, int32_t, uint8_t
# Note: _tree uses cimport numpy, cnp.import_array, so we need to include
# numpy headers in the build configuration of this extension
from sklearn.tree._tree cimport Node
from sklearn.tree._tree cimport Tree
from sklearn.tree._utils cimport safe_realloc


# no namespace lookup for numpy dtype and array creation
from numpy import zeros as np_zeros


# constant to mark tree leafs
cdef intp_t TREE_LEAF = -1

cdef void _predict_regression_tree_inplace_fast_dense(
    const float32_t[:, ::1] X,
    Node* root_node,
    double *value,
    double scale,
    Py_ssize_t k,
    float64_t[:, :] out
) noexcept nogil:
    """Predicts output for regression tree and stores it in ``out[i, k]``.

    This function operates directly on the data arrays of the tree
    data structures. This is 5x faster than the variant above because
    it allows us to avoid buffer validation.

    The function assumes that the ndarray that wraps ``X`` is
    c-continuous.

    Parameters
    ----------
    X : float32_t 2d memory view
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
    out : memory view on array of type np.float64_t
        The data array where the predictions are stored.
        ``out`` is assumed to be a two-dimensional array of
        shape ``(n_samples, K)``.
    """
    cdef intp_t n_samples = X.shape[0]
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
        out[i, k] += scale * value[node - root_node]


def _predict_regression_tree_stages_sparse(
    object[:, :] estimators,
    object X,
    double scale,
    float64_t[:, :] out
):
    """Predicts output for regression tree inplace and adds scaled value to ``out[i, k]``.

    The function assumes that the ndarray that wraps ``X`` is csr_matrix.
    """
    cdef const float32_t[::1] X_data = X.data
    cdef const int32_t[::1] X_indices = X.indices
    cdef const int32_t[::1] X_indptr = X.indptr

    cdef intp_t n_samples = X.shape[0]
    cdef intp_t n_features = X.shape[1]
    cdef intp_t n_stages = estimators.shape[0]
    cdef intp_t n_outputs = estimators.shape[1]

    # Indices and temporary variables
    cdef intp_t sample_i
    cdef intp_t feature_i
    cdef intp_t stage_i
    cdef intp_t output_i
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
    cdef float32_t feature_value = 0.
    cdef float32_t* X_sample = NULL

    # feature_to_sample as a data structure records the last seen sample
    # for each feature; functionally, it is an efficient way to identify
    # which features are nonzero in the present sample.
    cdef intp_t* feature_to_sample = NULL

    safe_realloc(&X_sample, n_features)
    safe_realloc(&feature_to_sample, n_features)

    memset(feature_to_sample, -1, n_features * sizeof(intp_t))

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
    float64_t[:, :] out
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
                    out=out
                )
                # out[:, k] += scale * tree.predict(X).ravel()


def predict_stage(
    object[:, :] estimators,
    int stage,
    object X,
    double scale,
    float64_t[:, :] out
):
    """Add predictions of ``estimators[stage]`` to ``out``.

    Each estimator in the stage is scaled by ``scale`` before
    its prediction is added to ``out``.
    """
    return predict_stages(
        estimators=estimators[stage:stage + 1], X=X, scale=scale, out=out
    )


def _random_sample_mask(
    intp_t n_total_samples,
    intp_t n_total_in_bag,
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
    cdef float64_t[::1] rand = random_state.uniform(size=n_total_samples)
    cdef uint8_t[::1] sample_mask = np_zeros((n_total_samples,), dtype=bool)

    cdef intp_t n_bagged = 0
    cdef intp_t i = 0

    for i in range(n_total_samples):
        if rand[i] * (n_total_samples - i) < (n_total_in_bag - n_bagged):
            sample_mask[i] = 1
            n_bagged += 1

    return sample_mask.base
