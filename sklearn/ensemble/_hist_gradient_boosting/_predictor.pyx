# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from cython.parallel import prange
from libc.stdlib cimport calloc, free
import numpy as np

from sklearn.utils._bitset cimport BITSET_INNER_DTYPE_C, in_bitset_2d_memoryview
from sklearn.utils._typedefs cimport intp_t, uint8_t, inlinable_isnan
from sklearn.ensemble._hist_gradient_boosting.common cimport X_DTYPE_C
from sklearn.ensemble._hist_gradient_boosting.common cimport Y_DTYPE_C
from sklearn.ensemble._hist_gradient_boosting.common import Y_DTYPE
from sklearn.ensemble._hist_gradient_boosting.common cimport X_BINNED_DTYPE_C
from sklearn.ensemble._hist_gradient_boosting.common cimport node_struct


def _predict_from_raw_data(  # raw data = non-binned data
        const node_struct [:] nodes,
        const X_DTYPE_C [:, :] numeric_data,
        const BITSET_INNER_DTYPE_C [:, ::1] raw_left_cat_bitsets,
        const BITSET_INNER_DTYPE_C [:, ::1] known_cat_bitsets,
        const unsigned int [::1] f_idx_map,
        int n_threads,
        Y_DTYPE_C [:] out):

    cdef:
        int i

    for i in prange(numeric_data.shape[0], schedule='static', nogil=True,
                    num_threads=n_threads):
        out[i] = _predict_one_from_raw_data(
            nodes, numeric_data, raw_left_cat_bitsets,
            known_cat_bitsets,
            f_idx_map, i)


cdef inline Y_DTYPE_C _predict_one_from_raw_data(
        const node_struct [:] nodes,
        const X_DTYPE_C [:, :] numeric_data,
        const BITSET_INNER_DTYPE_C [:, ::1] raw_left_cat_bitsets,
        const BITSET_INNER_DTYPE_C [:, ::1] known_cat_bitsets,
        const unsigned int [::1] f_idx_map,
        const int row) noexcept nogil:
    # Need to pass the whole array and the row index, else prange won't work.
    # See issue Cython #2798

    cdef:
        node_struct node = nodes[0]
        unsigned int node_idx = 0
        X_DTYPE_C data_val

    while True:
        if node.is_leaf:
            return node.value

        data_val = numeric_data[row, node.feature_idx]

        if inlinable_isnan(data_val):
            if node.missing_go_to_left:
                node_idx = node.left
            else:
                node_idx = node.right
        elif node.is_categorical:
            if in_bitset_2d_memoryview(
                    raw_left_cat_bitsets,
                    <X_BINNED_DTYPE_C>data_val,
                    node.bitset_idx):
                node_idx = node.left
            elif in_bitset_2d_memoryview(
                    known_cat_bitsets,
                    <X_BINNED_DTYPE_C>data_val,
                    f_idx_map[node.feature_idx]):
                node_idx = node.right
            else:
                # Treat unknown categories as missing.
                node_idx = node.left if node.missing_go_to_left else node.right
        else:
            if data_val <= node.num_threshold:
                node_idx = node.left
            else:
                node_idx = node.right
        node = nodes[node_idx]


def _predict_from_binned_data(
        node_struct [:] nodes,
        const X_BINNED_DTYPE_C [:, :] binned_data,
        BITSET_INNER_DTYPE_C [:, :] binned_left_cat_bitsets,
        const uint8_t missing_values_bin_idx,
        int n_threads,
        Y_DTYPE_C [:] out):

    cdef:
        int i

    for i in prange(binned_data.shape[0], schedule='static', nogil=True,
                    num_threads=n_threads):
        out[i] = _predict_one_from_binned_data(nodes,
                                               binned_data,
                                               binned_left_cat_bitsets, i,
                                               missing_values_bin_idx)


cdef inline Y_DTYPE_C _predict_one_from_binned_data(
        node_struct [:] nodes,
        const X_BINNED_DTYPE_C [:, :] binned_data,
        const BITSET_INNER_DTYPE_C [:, :] binned_left_cat_bitsets,
        const int row,
        const uint8_t missing_values_bin_idx) noexcept nogil:
    # Need to pass the whole array and the row index, else prange won't work.
    # See issue Cython #2798

    cdef:
        node_struct node = nodes[0]
        unsigned int node_idx = 0
        X_BINNED_DTYPE_C data_val

    while True:
        if node.is_leaf:
            return node.value

        data_val = binned_data[row, node.feature_idx]

        if data_val == missing_values_bin_idx:
            if node.missing_go_to_left:
                node_idx = node.left
            else:
                node_idx = node.right
        elif node.is_categorical:
            if in_bitset_2d_memoryview(
                    binned_left_cat_bitsets,
                    data_val,
                    node.bitset_idx):
                node_idx = node.left
            else:
                node_idx = node.right
        else:
            if data_val <= node.bin_threshold:
                node_idx = node.left
            else:
                node_idx = node.right
        node = nodes[node_idx]


def _compute_partial_dependence(
    node_struct [:] nodes,
    const X_DTYPE_C [:, ::1] X,
    const intp_t [:] target_features,
    Y_DTYPE_C [:] out
):
    """Partial dependence of the response on the ``target_features`` set.

    For each sample in ``X`` a tree traversal is performed.
    Each traversal starts from the root with weight 1.0.

    At each non-leaf node that splits on a target feature, either
    the left child or the right child is visited based on the feature
    value of the current sample, and the weight is not modified.
    At each non-leaf node that splits on a complementary feature,
    both children are visited and the weight is multiplied by the fraction
    of training samples which went to each child.

    At each leaf, the value of the node is multiplied by the current
    weight (weights sum to 1 for all visited terminal nodes).

    Parameters
    ----------
    nodes : view on array of PREDICTOR_RECORD_DTYPE, shape (n_nodes)
        The array representing the predictor tree.
    X : view on 2d ndarray, shape (n_samples, n_target_features)
        The grid points on which the partial dependence should be
        evaluated.
    target_features : view on 1d ndarray of intp_t, shape (n_target_features)
        The set of target features for which the partial dependence
        should be evaluated.
    out : view on 1d ndarray, shape (n_samples)
        The value of the partial dependence function on each grid
        point.
    """

    cdef:
        unsigned int current_node_idx
        unsigned int [:] node_idx_stack = np.zeros(shape=nodes.shape[0],
                                                   dtype=np.uint32)
        Y_DTYPE_C [::1] weight_stack = np.zeros(shape=nodes.shape[0],
                                                dtype=Y_DTYPE)
        node_struct * current_node  # pointer to avoid copying attributes

        unsigned int sample_idx
        intp_t feature_idx
        unsigned stack_size
        Y_DTYPE_C left_sample_frac
        Y_DTYPE_C current_weight
        Y_DTYPE_C total_weight  # used for sanity check only
        bint is_target_feature

    for sample_idx in range(X.shape[0]):
        # init stacks for current sample
        stack_size = 1
        node_idx_stack[0] = 0  # root node
        weight_stack[0] = 1  # all the samples are in the root node
        total_weight = 0

        while stack_size > 0:

            # pop the stack
            stack_size -= 1
            current_node_idx = node_idx_stack[stack_size]
            current_node = &nodes[current_node_idx]

            if current_node.is_leaf:
                out[sample_idx] += (weight_stack[stack_size] *
                                    current_node.value)
                total_weight += weight_stack[stack_size]
            else:
                # determine if the split feature is a target feature
                is_target_feature = False
                for feature_idx in range(target_features.shape[0]):
                    if target_features[feature_idx] == current_node.feature_idx:
                        is_target_feature = True
                        break

                if is_target_feature:
                    # In this case, we push left or right child on stack
                    if X[sample_idx, feature_idx] <= current_node.num_threshold:
                        node_idx_stack[stack_size] = current_node.left
                    else:
                        node_idx_stack[stack_size] = current_node.right
                    stack_size += 1
                else:
                    # In this case, we push both children onto the stack,
                    # and give a weight proportional to the number of
                    # samples going through each branch.

                    # push left child
                    node_idx_stack[stack_size] = current_node.left
                    left_sample_frac = (
                        <Y_DTYPE_C> nodes[current_node.left].count /
                        current_node.count)
                    current_weight = weight_stack[stack_size]
                    weight_stack[stack_size] = current_weight * left_sample_frac
                    stack_size += 1

                    # push right child
                    node_idx_stack[stack_size] = current_node.right
                    weight_stack[stack_size] = (
                        current_weight * (1 - left_sample_frac))
                    stack_size += 1

        # Sanity check. Should never happen.
        if not (0.999 < total_weight < 1.001):
            raise ValueError("Total weight should be 1.0 but was %.9f" %total_weight)


def _check_nodes(
    const node_struct [:] nodes,
    const intp_t n_features,
    const uint8_t [:] is_categorical,
):
    """Check the values of a node array in a single pass over the nodes.

    See `TreePredictor.check_state` for what is checked and why.

    Parameters
    ----------
    nodes : view on array of PREDICTOR_RECORD_DTYPE, shape (n_nodes)
        The array representing the predictor tree.
    n_features : int
        The number of columns of `X` at prediction time.
    is_categorical : view on 1d ndarray of uint8, shape (n_features)
        Whether each feature is categorical.

    Returns
    -------
    invalid : str or None
        What the first invalid split node gets wrong: the name of the field
        that is out of bounds, "shared child" if one of its children already
        is the child of a node, "numerical feature" if it is a categorical
        split on a feature that is not categorical. None if all nodes are
        valid.
    n_bitsets : int
        The number of rows the categorical splits need in the bitset arrays.
    """
    cdef:
        intp_t n_nodes = nodes.shape[0]
        intp_t i, left, right
        intp_t n_bitsets = 0
        uint8_t* has_parent
        str invalid = None

    if is_categorical.shape[0] != n_features:
        raise ValueError(
            f"is_categorical has {is_categorical.shape[0]} entries, expected "
            f"{n_features}."
        )
    if n_nodes == 0:
        return None, 0

    has_parent = <uint8_t*> calloc(n_nodes, sizeof(uint8_t))
    if has_parent == NULL:
        raise MemoryError()

    for i in range(n_nodes):
        if nodes[i].is_leaf:
            continue

        if nodes[i].feature_idx < 0 or nodes[i].feature_idx >= n_features:
            invalid = "feature_idx"
            break
        # `left` / `right` must be node ids strictly greater than the own index:
        # this keeps them in `[0, n_nodes)` and rules out cycles.
        left = nodes[i].left
        right = nodes[i].right
        if left <= i or left >= n_nodes:
            invalid = "left"
            break
        if right <= i or right >= n_nodes:
            invalid = "right"
            break
        if left == right or has_parent[left] or has_parent[right]:
            invalid = "shared child"
            break
        has_parent[left] = 1
        has_parent[right] = 1

        if nodes[i].is_categorical:
            if not is_categorical[nodes[i].feature_idx]:
                invalid = "numerical feature"
                break
            if nodes[i].bitset_idx >= n_bitsets:
                n_bitsets = <intp_t> nodes[i].bitset_idx + 1

    free(has_parent)
    return invalid, n_bitsets
