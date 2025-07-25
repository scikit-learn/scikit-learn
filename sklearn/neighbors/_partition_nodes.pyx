# BinaryTrees rely on partial sorts to partition their nodes during their
# initialisation.
#
# The C++ std library exposes nth_element, an efficient partial sort for this
# situation which has a linear time complexity as well as the best performances.
#
# To use std::algorithm::nth_element, a few fixture are defined using Cython:
# - partition_node_indices, a Cython function used in BinaryTrees, that calls
# - partition_node_indices_inner, a C++ function that wraps nth_element and uses
# - an IndexComparator to state how to compare KDTrees' indices
#
# IndexComparator has been defined so that partial sorts are stable with
# respect to the nodes initial indices.
#
# See for reference:
#  - https://en.cppreference.com/w/cpp/algorithm/nth_element.
#  - https://github.com/scikit-learn/scikit-learn/pull/11103
#  - https://github.com/scikit-learn/scikit-learn/pull/19473
from cython cimport floating


cdef extern from *:
    """
    #include <algorithm>

    template<class D, class I>
    class IndexComparator {
    private:
        const D *data;
        I split_dim, n_features;
    public:
        IndexComparator(const D *data, const I &split_dim, const I &n_features):
            data(data), split_dim(split_dim), n_features(n_features) {}

        bool operator()(const I &a, const I &b) const {
            D a_value = data[a * n_features + split_dim];
            D b_value = data[b * n_features + split_dim];
            return a_value == b_value ? a < b : a_value < b_value;
        }
    };

    template<class D, class I>
    void partition_node_indices_inner(
        const D *data,
        I *node_indices,
        const I &split_dim,
        const I &split_index,
        const I &n_features,
        const I &n_points) {
        IndexComparator<D, I> index_comparator(data, split_dim, n_features);
        std::nth_element(
            node_indices,
            node_indices + split_index,
            node_indices + n_points,
            index_comparator);
    }
    """
    void partition_node_indices_inner[D, I](
                const D *data,
                I *node_indices,
                I split_dim,
                I split_index,
                I n_features,
                I n_points) except +


cdef int partition_node_indices(
        const floating *data,
        intp_t *node_indices,
        intp_t split_dim,
        intp_t split_index,
        intp_t n_features,
        intp_t n_points) except -1:
    """Partition points in the node into two equal-sized groups.

    Upon return, the values in node_indices will be rearranged such that
    (assuming numpy-style indexing):

        data[node_indices[0:split_index], split_dim]
          <= data[node_indices[split_index], split_dim]

    and

        data[node_indices[split_index], split_dim]
          <= data[node_indices[split_index:n_points], split_dim]

    The algorithm is essentially a partial in-place quicksort around a
    set pivot.

    Parameters
    ----------
    data : double pointer
        Pointer to a 2D array of the training data, of shape [N, n_features].
        N must be greater than any of the values in node_indices.
    node_indices : int pointer
        Pointer to a 1D array of length n_points.  This lists the indices of
        each of the points within the current node.  This will be modified
        in-place.
    split_dim : int
        the dimension on which to split.  This will usually be computed via
        the routine ``find_node_split_dim``.
    split_index : int
        the index within node_indices around which to split the points.
    n_features: int
        the number of features (i.e columns) in the 2D array pointed by data.
    n_points : int
        the length of node_indices. This is also the number of points in
        the original dataset.
    Returns
    -------
    status : int
        integer exit status.  On return, the contents of node_indices are
        modified as noted above.
    """
    partition_node_indices_inner(
        data,
        node_indices,
        split_dim,
        split_index,
        n_features,
        n_points)
    return 0
