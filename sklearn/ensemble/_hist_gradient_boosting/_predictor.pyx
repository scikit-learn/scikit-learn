# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3

# Author: Nicolas Hug

cimport cython
from cython.parallel import prange
from libc.math cimport isnan
import numpy as np
cimport numpy as np
from numpy.math cimport INFINITY

from .common cimport X_DTYPE_C
from .common cimport Y_DTYPE_C
from .common import Y_DTYPE
from .common cimport X_BINNED_DTYPE_C
from .common cimport node_struct

np.import_array()


@cython.final
cdef class TreePredictor:
    """Tree class used for predictions.

    Parameters
    ----------
    n_nodes : int
        Number of nodes in tree.
    """
    def __init__(self, int n_nodes):
        self.nodes.resize(n_nodes)

    cdef node_struct* get(self, int node_idx) nogil:
        """Get reference to node."""
        return &self.nodes[node_idx]

    cdef int get_n_nodes(self) nogil:
        """Get number of nodes in tree."""
        return self.nodes.size()

    def get_n_leaf_nodes(self):
        """Return number of leaves."""
        cdef:
            int n_leaf_nodes = 0
            node_struct node

        for node in self.nodes:
            n_leaf_nodes += node.is_leaf
        return n_leaf_nodes

    def get_max_depth(self):
        """Return maximum depth among all leaves."""
        cdef:
            int max_depth = 0
            node_struct node

        for node in self.nodes:
            max_depth = max(node.depth, max_depth)
        return max_depth

    def predict(self, const X_DTYPE_C [:, :] X):
        """Predict raw values for non-binned data.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The raw predicted values.
        """
        cdef:
            int i
            Y_DTYPE_C [:] out = np.empty(X.shape[0], dtype=Y_DTYPE)

        for i in prange(X.shape[0], schedule='static', nogil=True):
            out[i] = _predict_one_from_numeric_data(self, X, i)
        return np.asarray(out)

    def predict_binned(self, const X_BINNED_DTYPE_C [:, :] X,
                       const unsigned char missing_values_bin_idx):
        """Predict raw values for binned data.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The input samples.
        missing_values_bin_idx : uint8
            Index of the bin that is used for missing values. This is the
            index of the last bin and is always equal to max_bins (as passed
            to the GBDT classes), or equivalently to n_bins - 1.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The raw predicted values.
        """
        cdef:
            int i
            Y_DTYPE_C [:] out = np.empty(X.shape[0], dtype=Y_DTYPE)

        for i in prange(X.shape[0], schedule='static', nogil=True):
            out[i] = _predict_one_from_binned_data(self, X, i,
                                                   missing_values_bin_idx)
        return np.asarray(out)

    def compute_partial_dependence(self, grid, target_features, out):
        """Fast partial dependence computation.

        Parameters
        ----------
        grid : ndarray, shape (n_samples, n_target_features)
            The grid points on which the partial dependence should be
            evaluated.
        target_features : ndarray, shape (n_target_features)
            The set of target features for which the partial dependence
            should be evaluated.
        out : ndarray, shape (n_samples)
            The value of the partial dependence function on each grid
            point.
        """
        _compute_partial_dependence(self, grid, target_features, out)

    def _set_node(self, int node_idx, unsigned int left=0,
                  unsigned int right=0, Y_DTYPE_C value=0.,
                  unsigned char is_leaf=False,
                  unsigned int feature_idx=0, X_DTYPE_C threshold=0.0):
        """Used for testing for setting the value of a node."""
        cdef:
            node_struct * node = &self.nodes[node_idx]
        node.left = left
        node.right = right
        node.value = value
        node.feature_idx = feature_idx
        node.threshold = threshold
        node.is_leaf = is_leaf


cdef inline Y_DTYPE_C _predict_one_from_numeric_data(
        TreePredictor tree,
        const X_DTYPE_C [:, :] numeric_data,
        const int row) nogil:
    # Need to pass the whole array and the row index, else prange won't work.
    # See issue Cython #2798

    cdef:
        node_struct * node = tree.get(0)

    while True:
        if node.is_leaf:
            return node.value

        if isnan(numeric_data[row, node.feature_idx]):
            if node.missing_go_to_left:
                node = tree.get(node.left)
            else:
                node = tree.get(node.right)
        else:
            if numeric_data[row, node.feature_idx] <= node.threshold:
                node = tree.get(node.left)
            else:
                node = tree.get(node.right)

cdef inline Y_DTYPE_C _predict_one_from_binned_data(
        TreePredictor tree,
        const X_BINNED_DTYPE_C [:, :] binned_data,
        const int row,
        const unsigned char missing_values_bin_idx) nogil:
    # Need to pass the whole array and the row index, else prange won't work.
    # See issue Cython #2798

    cdef:
        node_struct * node = tree.get(0)

    while True:
        if node.is_leaf:
            return node.value
        if binned_data[row, node.feature_idx] ==  missing_values_bin_idx:
            if node.missing_go_to_left:
                node = tree.get(node.left)
            else:
                node = tree.get(node.right)
        else:
            if binned_data[row, node.feature_idx] <= node.bin_threshold:
                node = tree.get(node.left)
            else:
                node = tree.get(node.right)

def _compute_partial_dependence(
    TreePredictor tree,
    const X_DTYPE_C [:, ::1] X,
    int [:] target_features,
    Y_DTYPE_C [:] out):
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
    weight (weights sum to 1 for all visited terminal tree).

    Parameters
    ----------
    tree : TreePredictor
        Object representing the predictor tree.
    X : view on 2d ndarray, shape (n_samples, n_target_features)
        The grid points on which the partial dependence should be
        evaluated.
    target_features : view on 1d ndarray, shape (n_target_features)
        The set of target features for which the partial dependence
        should be evaluated.
    out : view on 1d ndarray, shape (n_samples)
        The value of the partial dependence function on each grid
        point.
    """

    cdef:
        unsigned int current_node_idx
        unsigned int [:] node_idx_stack = np.zeros(shape=tree.get_n_nodes(),
                                                   dtype=np.uint32)
        Y_DTYPE_C [::1] weight_stack = np.zeros(shape=tree.get_n_nodes(),
                                                dtype=Y_DTYPE)
        node_struct * current_node  # pointer to avoid copying attributes

        unsigned int sample_idx
        unsigned feature_idx
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
            current_node = tree.get(current_node_idx)

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
                    if X[sample_idx, feature_idx] <= current_node.threshold:
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
                        <Y_DTYPE_C> tree.get(current_node.left).count /
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
            raise ValueError("Total weight should be 1.0 but was %.9f" %
                                total_weight)
