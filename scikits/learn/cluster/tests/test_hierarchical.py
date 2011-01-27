"""
Several basic tests for hierarchical clustering procedures

Author : Vincent Michel, 2010
"""

import numpy as np
from scikits.learn.cluster import Ward, WardAgglomeration, ward_tree
from scikits.learn.feature_extraction.image import img_to_graph


def test_structured_ward_tree():
    """
    Check that we obtain the correct solution for structured ward tree.
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(50, 100)
    connectivity = img_to_graph(mask, mask)
    children, n_comp, n_leaves = ward_tree(X.T, connectivity)
    n_nodes = 2 * X.shape[1] - 1
    assert(len(children) + n_leaves == n_nodes)


def test_unstructured_ward_tree():
    """
    Check that we obtain the correct solution for unstructured ward tree.
    """
    np.random.seed(0)
    X = np.random.randn(50, 100)
    children, n_nodes, n_leaves = ward_tree(X.T)
    n_nodes = 2 * X.shape[1] - 1
    assert(len(children) + n_leaves == n_nodes)


def test_height_ward_tree():
    """
    Check that the height of ward tree is sorted.
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(50, 100)
    connectivity = img_to_graph(mask, mask)
    children, n_nodes, n_leaves = ward_tree(X.T, connectivity)
    n_nodes = 2 * X.shape[1] - 1
    assert(len(children) + n_leaves == n_nodes)


def test_ward_clustering():
    """
    Check that we obtain the correct number of clusters with Ward clustering.
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(100, 50)
    connectivity = img_to_graph(mask, mask)
    clustering = Ward(n_clusters=10, connectivity=connectivity)
    clustering.fit(X)
    assert(np.size(np.unique(clustering.labels_)) == 10)


def test_ward_agglomeration():
    """
    Check that we obtain the correct solution in a simplistic case
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(50, 100)
    connectivity = img_to_graph(mask, mask)
    ward = WardAgglomeration(n_clusters=5, connectivity=connectivity)
    ward.fit(X)
    assert(np.size(np.unique(ward.labels_)) == 5)

    Xred = ward.transform(X)
    assert(Xred.shape[1] == 5)
    Xfull = ward.inverse_transform(Xred)
    assert(np.unique(Xfull[0]).size == 5)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
