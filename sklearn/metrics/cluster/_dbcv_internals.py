import numpy as np
from scipy.spatial.distance import cdist

from ..pairwise import pairwise_distances
from ._dbcv_linkage import mst_linkage_core


def density_separation(
    X,
    labels,
    cluster_id1,
    cluster_id2,
    internal_nodes1,
    internal_nodes2,
    core_distances1,
    core_distances2,
    metric="euclidean",
    no_coredist=False,
    **kwd_args,
):
    """
    Compute the density separation between two clusters. This is the minimum
    distance between pairs of points, one from internal nodes of MSTs of each cluster.

    Parameters
    ----------
    X : array (n_samples, n_features) or (n_samples, n_samples)
        The input data of the clustering. This can be the data, or, if
        metric is set to `precomputed` the pairwise distance matrix used
        for the clustering.

    labels : array (n_samples)
        The label array output by the clustering, providing an integral
        cluster label to each data point, with -1 for noise points.

    cluster_id1 : integer
        The first cluster label to compute separation between.

    cluster_id2 : integer
        The second cluster label to compute separation between.

    internal_nodes1 : array
        The vertices of the MST for `cluster_id1` that were internal vertices.

    internal_nodes2 : array
        The vertices of the MST for `cluster_id2` that were internal vertices.

    core_distances1 : array (size of cluster_id1,)
        The all-points-core_distances of all points in the cluster
        specified by cluster_id1.

    core_distances2 : array (size of cluster_id2,)
        The all-points-core_distances of all points in the cluster
        specified by cluster_id2.

    metric : string
        The metric used to compute distances for the clustering (and
        to be re-used in computing distances for mr distance). If
        set to `precomputed` then X is assumed to be the precomputed
        distance matrix between samples.

    **kwd_args :
        Extra arguments to pass to the distance computation for other
        metrics, such as minkowski, Mahanalobis etc.

    Returns
    -------
    The 'density separation' between the clusters specified by
    `cluster_id1` and `cluster_id2`.

    References
    ----------
    Moulavi, D., Jaskowiak, P.A., Campello, R.J., Zimek, A. and Sander, J.,
    2014. Density-Based Clustering Validation. In SDM (pp. 839-847).
    """
    if metric == "precomputed":
        sub_select = X[labels == cluster_id1, :][:, labels == cluster_id2]
        distance_matrix = sub_select[internal_nodes1, :][:, internal_nodes2]
    else:
        cluster1 = X[labels == cluster_id1][internal_nodes1]
        cluster2 = X[labels == cluster_id2][internal_nodes2]
        distance_matrix = cdist(cluster1, cluster2, metric, **kwd_args)

    if no_coredist:
        return distance_matrix.min()

    else:
        core_dist_matrix1 = np.tile(
            core_distances1[internal_nodes1], (distance_matrix.shape[1], 1)
        ).T
        core_dist_matrix2 = np.tile(
            core_distances2[internal_nodes2], (distance_matrix.shape[0], 1)
        )

        mr_dist_matrix = np.dstack(
            [distance_matrix, core_dist_matrix1, core_dist_matrix2]
        ).max(axis=-1)

        return mr_dist_matrix.min()


def internal_minimum_spanning_tree(mr_distances):
    """
    Compute the 'internal' minimum spanning tree given a matrix of mutual
    reachability distances. Given a minimum spanning tree the 'internal'
    graph is the subgraph induced by vertices of degree greater than one.

    Parameters
    ----------
    mr_distances : array (cluster_size, cluster_size)
        The pairwise mutual reachability distances, inferred to be the edge
        weights of a complete graph. Since MSTs are computed per cluster
        this is the all-points-mutual-reacability for points within a single
        cluster.

    Returns
    -------
    internal_nodes : array
        An array listing the indices of the internal nodes of the MST

    internal_edges : array (?, 3)
        An array of internal edges in weighted edge list format; that is
        an edge is an array of length three listing the two vertices
        forming the edge and weight of the edge.

    References
    ----------
    Moulavi, D., Jaskowiak, P.A., Campello, R.J., Zimek, A. and Sander, J.,
    2014. Density-Based Clustering Validation. In SDM (pp. 839-847).
    """
    single_linkage_data = mst_linkage_core(mr_distances)
    min_span_tree = single_linkage_data.copy()
    for index, row in enumerate(min_span_tree[1:], 1):
        candidates = np.where(np.isclose(mr_distances[int(row[1])], row[2]))[0]
        candidates = np.intersect1d(
            candidates, single_linkage_data[:index, :2].astype(int)
        )
        candidates = candidates[candidates != row[1]]
        assert len(candidates) > 0
        row[0] = candidates[0]

    vertices = np.arange(mr_distances.shape[0])[
        np.bincount(min_span_tree.T[:2].flatten().astype(np.intp)) > 1
    ]
    if not len(vertices):
        vertices = [0]
    # A little "fancy" we select from the flattened array reshape back
    # (Fortran format to get indexing right) and take the product to do an and
    # then convert back to boolean type.
    edge_selection = np.prod(
        np.in1d(min_span_tree.T[:2], vertices).reshape(
            (min_span_tree.shape[0], 2), order="F"
        ),
        axis=1,
    ).astype(bool)

    # Density sparseness is not well defined if there are no
    # internal edges (as per the referenced paper). However
    # MATLAB code from the original authors simply selects the
    # largest of *all* the edges in the case that there are
    # no internal edges, so we do the same here
    if np.any(edge_selection):
        # If there are any internal edges, then subselect them out
        edges = min_span_tree[edge_selection]
    else:
        # If there are no internal edges then we want to take the
        # max over all the edges that exist in the MST, so we simply
        # do nothing and return all the edges in the MST.
        edges = min_span_tree.copy()

    return vertices, edges


def distances_between_points(
    X,
    labels,
    cluster_id,
    metric="euclidean",
    d=None,
    no_coredist=False,
    print_max_raw_to_coredist_ratio=False,
    **kwd_args,
):
    """
    Compute pairwise distances for all the points of a cluster.

    If metric is 'precomputed' then assume X is a distance matrix for the full
    dataset. Note that in this case you must pass in 'd' the dimension of the
    dataset.

    Parameters
    ----------
    X : array (n_samples, n_features) or (n_samples, n_samples)
        The input data of the clustering. This can be the data, or, if
        metric is set to `precomputed` the pairwise distance matrix used
        for the clustering.

    labels : array (n_samples)
        The label array output by the clustering, providing an integral
        cluster label to each data point, with -1 for noise points.

    cluster_id : integer
        The cluster label for which to compute the distances

    metric : string
        The metric used to compute distances for the clustering (and
        to be re-used in computing distances for mr distance). If
        set to `precomputed` then X is assumed to be the precomputed
        distance matrix between samples.

    d : integer (or None)
        The number of features (dimension) of the dataset. This need only
        be set in the case of metric being set to `precomputed`, where
        the ambient dimension of the data is unknown to the function.

    **kwd_args :
        Extra arguments to pass to the distance computation for other
        metrics, such as minkowski, Mahanalobis etc.

    Returns
    -------

    distances : array (n_samples, n_samples)
        The distances between all points in `X` with `label` equal to `cluster_id`.

    core_distances : array (n_samples,)
        The all-points-core_distance of all points in `X` with `label` equal
        to `cluster_id`.

    References
    ----------
    Moulavi, D., Jaskowiak, P.A., Campello, R.J., Zimek, A. and Sander, J.,
    2014. Density-Based Clustering Validation. In SDM (pp. 839-847).
    """
    if metric == "precomputed":
        if d is None:
            raise ValueError("If metric is precomputed a d value must be provided!")
        distance_matrix = X[labels == cluster_id, :][:, labels == cluster_id]
    else:
        subset_X = X[labels == cluster_id, :]
        distance_matrix = pairwise_distances(subset_X, metric=metric, **kwd_args)
        d = X.shape[1]

    if no_coredist:
        return distance_matrix, None

    else:
        core_distances = all_points_core_distance(distance_matrix.copy(), d=d)
        core_dist_matrix = np.tile(core_distances, (core_distances.shape[0], 1))
        stacked_distances = np.dstack(
            [distance_matrix, core_dist_matrix, core_dist_matrix.T]
        )

        if print_max_raw_to_coredist_ratio:
            print(
                "Max raw distance to coredistance ratio: "
                + str(max_ratio(stacked_distances))
            )

        return stacked_distances.max(axis=-1), core_distances


def all_points_core_distance(distance_matrix, d=2.0):
    """
    Compute the all-points-core-distance for all the points of a cluster.

    Parameters
    ----------
    distance_matrix : array (cluster_size, cluster_size)
        The pairwise distance matrix between points in the cluster.

    d : integer
        The dimension of the data set, which is used in the computation
        of the all-point-core-distance as per the paper.

    Returns
    -------
    core_distances : array (cluster_size,)
        The all-points-core-distance of each point in the cluster

    References
    ----------
    Moulavi, D., Jaskowiak, P.A., Campello, R.J., Zimek, A. and Sander, J.,
    2014. Density-Based Clustering Validation. In SDM (pp. 839-847).
    """
    distance_matrix[distance_matrix != 0] = (
        1.0 / distance_matrix[distance_matrix != 0]
    ) ** d
    result = distance_matrix.sum(axis=1)
    result /= distance_matrix.shape[0] - 1

    if result.sum() == 0:
        result = np.zeros(len(distance_matrix))
    else:
        result **= -1.0 / d

    return result


def max_ratio(stacked_distances):
    max_ratio = 0
    for i in range(stacked_distances.shape[0]):
        for j in range(stacked_distances.shape[1]):
            dist = stacked_distances[i][j][0]
            coredist = stacked_distances[i][j][1]
            if dist == 0 or coredist / dist <= max_ratio:
                continue
            max_ratio = coredist / dist

    return max_ratio
