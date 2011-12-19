"""
Several basic tests for hierarchical clustering procedures

Author : Jan Hendrik Metzen, 2011
"""

import random
import numpy as np

from sklearn.neighbors import kneighbors_graph
from sklearn.cluster.linkage import WardsLinkage, CompleteLinkage


def ward_distance(cluster1, cluster2):
    """ Compute distance of two cluster according to Ward's criterion.

    According to this criterion, at each step during agglomerative clustering,
    the pair of clusters with minimum cluster distance are merged, where the
    cluster distance is the increase of within-cluster variance by merging
    two clusters.
    """
    centroid1 = np.mean(cluster1, axis=0)
    centroid2 = np.mean(cluster2, axis=0)
    dist = np.linalg.norm(centroid1 - centroid2) ** 2 / \
                        (1.0 / cluster1.shape[0] + 1.0 / cluster2.shape[0])
    return dist


def complete_linkage_distance(cluster1, cluster2, metric):
    """ Compute complete linkage distance based on brute force search. """
    if metric == "euclidean":
        metric = lambda x, y: np.linalg.norm(x - y)
    nodes = np.vstack((cluster1, cluster2))
    dist = 0.0
    for i in range(nodes.shape[0]):
        for j in range(i + 1, nodes.shape[0]):
            dist = max(dist, metric(nodes[i, :], nodes[j, :]))
    return dist


def test_wards_linkage():
    """
    Check that linkage based on ward's criterion computes distance of
    two clusterings correctly.
    """
    n_features = 2
    n_samples = 30
    n_nodes = 2 * n_samples - 1

    np.random.seed(0)
    X = np.random.randn(n_samples, n_features)
    connectivity = kneighbors_graph(X, n_neighbors=5).tolil()
    wards_linkage = WardsLinkage(X, connectivity)
    parent = np.arange(n_nodes, dtype=np.int)
    children = []
    open_nodes = np.ones(n_nodes, dtype=bool)

    k = n_samples
    while wards_linkage.has_more_candidates():
        linkage_distance, i, j = wards_linkage.fetch_candidate()
        if not open_nodes[i] or not open_nodes[j]:
            continue
        # Check that merge distance returned by wards_linkage matches the one
        # obtained by explicitly computing it on the two clusters
        cluster1 = X[wards_linkage.get_nodes_of_cluster(i), :]
        cluster2 = X[wards_linkage.get_nodes_of_cluster(j), :]
        true_dist = ward_distance(cluster1, cluster2)

        assert np.allclose([linkage_distance], [true_dist])

        # Continue with merging
        wards_linkage.update(i, j, k, parent)
        parent[i] = parent[j] = k
        children.append([i, j])
        open_nodes[i], open_nodes[j] = False, False
        k += 1


def test_complete_linkage(metric="euclidean"):
    """
    Check that complete linkage computes distance of two clusterings correctly.
    """
    n_features = 2
    n_samples = 30
    n_nodes = 2 * n_samples - 1

    np.random.seed(0)
    X = np.random.randn(n_samples, n_features)
    connectivity = kneighbors_graph(X, n_neighbors=5).tolil()
    complete_linkage = CompleteLinkage(X, connectivity,
                                       precompute_distances=True,
                                       metric=metric)
    parent = np.arange(n_nodes, dtype=np.int)
    children = []
    open_nodes = np.ones(n_nodes, dtype=bool)

    k = n_samples
    while complete_linkage.has_more_candidates():
        linkage_distance, i, j = complete_linkage.fetch_candidate()
        if not open_nodes[i] or not open_nodes[j]:
            continue
        # Check that merge distance returned by wards_linkage matches the one
        # obtained by explicitly computing it on the two clusters
        cluster1 = X[complete_linkage.get_nodes_of_cluster(i), :]
        cluster2 = X[complete_linkage.get_nodes_of_cluster(j), :]
        true_dist = complete_linkage_distance(cluster1, cluster2,
                                              metric)

        assert np.allclose([linkage_distance], [true_dist])

        # Continue with merging
        complete_linkage.update(i, j, k, parent)
        parent[i] = parent[j] = k
        children.append([i, j])
        open_nodes[i], open_nodes[j] = False, False
        k += 1


def test_complete_linkage_noneuclidean_distance():
    """
    Check that complete linkage computes distance of two clusterings correctly\
    when the pairwise distance is non-euclidean
    """
    metric = lambda x, y: random.Random(sum(x + y)).random()

    test_complete_linkage(metric=metric)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
