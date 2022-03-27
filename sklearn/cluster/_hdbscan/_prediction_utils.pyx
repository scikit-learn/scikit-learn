# Utility routines in cython for prediction in hdbscan
# Authors: Leland McInnes
# License: 3-clause BSD

import numpy as np
cimport numpy as np

from sklearn.metrics._dist_metrics cimport DistanceMetric

from libc.float cimport DBL_MAX
from libc.math cimport exp

cpdef get_tree_row_with_child(np.ndarray tree, np.intp_t child):

    cdef np.intp_t i
    cdef np.ndarray[np.intp_t, ndim = 1] child_array = tree['child']

    for i in range(tree.shape[0]):
        if child_array[i] == child:
            return tree[i]

    return tree[0]


cdef np.ndarray[np.float64_t, ndim=1] merge_height(
        np.intp_t point_cluster,
        np.float64_t point_lambda,
        np.ndarray[np.intp_t, ndim=1] clusters,
        np.ndarray cluster_tree):

    cdef np.intp_t i
    cdef np.intp_t j

    cdef np.intp_t left_cluster
    cdef np.intp_t right_cluster
    cdef int took_right_parent
    cdef int took_left_parent
    cdef np.intp_t cluster

    cdef np.ndarray[np.float64_t, ndim=1] result = np.empty(clusters.shape[0],
                                                            dtype=np.float64)
    cdef np.ndarray[np.intp_t, ndim=1] parents
    cdef np.ndarray[np.intp_t, ndim=1] children
    cdef np.ndarray[np.float64_t, ndim=1] lambdas

    # convert the cluster tree for fast direct access
    parents = cluster_tree['parent'].astype(np.intp)
    children = cluster_tree['child'].astype(np.intp)
    lambdas = cluster_tree['lambda_val'].astype(np.float64)


    for i in range(clusters.shape[0]):

        took_right_parent = False
        took_left_parent = False

        right_cluster = clusters[i]
        left_cluster = point_cluster

        while left_cluster != right_cluster:
            if left_cluster > right_cluster:
                took_left_parent = True
                last_cluster = left_cluster

                # Set left_cluster to be its parent
                for j in range(children.shape[0]):
                    if children[j] == left_cluster:
                        left_cluster = parents[j]
                        break
            else:
                took_right_parent = True
                last_cluster = right_cluster

                # Set right_cluster to be its parent
                for j in range(children.shape[0]):
                    if children[j] == right_cluster:
                        right_cluster = parents[j]
                        break

        if took_left_parent and took_right_parent:
            # Take the lambda value of last_cluster merging in
            for j in range(children.shape[0]):
                if children[j] == last_cluster:
                    result[i] = lambdas[j]
                    break
        else:
            result[i] = point_lambda

    return result


cpdef np.float64_t safe_always_positive_division(
        np.float64_t numerator,
        np.float64_t denominator):
    """ This is a helper function to divide numbers safely without getting a ZeroDivision error, the
    function handles zero division by assuming the denominator is always positive

    Parameters
    ----------
    numerator: floating
        any floating point type
    denominator: floating
        any floating point type

    Returns
    -------
    floating
    """
    if denominator <= 0:
        # prevent zero division or negative result
        denominator = 1e-8
    return numerator / denominator


cpdef np.ndarray[np.float64_t, ndim=1] per_cluster_scores(
        np.intp_t neighbor,
        np.float32_t lambda_,
        np.ndarray[np.intp_t, ndim=1] clusters,
        np.ndarray tree,
        dict max_lambda_dict,
        np.ndarray cluster_tree):

    cdef np.intp_t point_cluster
    cdef np.float64_t point_lambda
    cdef np.float64_t max_lambda

    cdef np.intp_t i

    cdef np.ndarray[np.float64_t, ndim=1] result

    point_row = get_tree_row_with_child(tree, neighbor)
    point_cluster = point_row['parent']
    point_lambda = lambda_
    max_lambda = max_lambda_dict[point_cluster]

    # Save an allocation by assigning and reusing result ...
    # height = merge_height(point_cluster, point_lambda,
    #                       clusters, cluster_tree)
    result = merge_height(point_cluster, point_lambda,
                          clusters, cluster_tree)

    # Cythonize: result = np.exp(-(max_lambda / height))
    for i in range(result.shape[0]):
        # result[i] = exp(-(max_lambda / result[i]))
        result[i] = safe_always_positive_division(max_lambda, (max_lambda - result[i]))

    return result


cpdef all_points_prob_in_some_cluster(
        np.ndarray[np.intp_t, ndim=1] clusters,
        np.ndarray tree,
        dict max_lambda_dict,
        np.ndarray cluster_tree):

    cdef np.ndarray[np.float64_t, ndim=1] heights
    cdef np.intp_t num_points = tree['parent'].min()
    cdef np.ndarray[np.float64_t, ndim=1] result
    cdef np.intp_t point
    cdef np.intp_t point_cluster
    cdef np.float64_t point_lambda
    cdef np.float64_t max_lambda

    cdef np.intp_t i

    result = np.empty(num_points, dtype=np.float64)

    point_tree = tree[tree['child_size'] == 1]

    for i in range(point_tree.shape[0]):
        point_row = point_tree[i]
        point = point_row['child']
        point_cluster = point_row['parent']
        point_lambda = point_row['lambda_val']

        # Can we not do a faster merge height operation here?
        heights = merge_height(point_cluster, point_lambda,
                               clusters, cluster_tree)
        max_lambda = max(max_lambda_dict[clusters[heights.argmax()]],
                         point_lambda)
        result[point] = (heights.max() / max_lambda)

    return result
