"""
Several basic tests for hierarchical clustering procedures

Author : Vincent Michel, 2010
"""

import numpy as np
from scikits.learn.cluster import Ward, ward_tree
from scikits.learn.feature_extraction.image import img_to_graph


def test_structured_ward_tree():
    """
    Check that we obtain the correct solution for structured ward tree.
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(50, 100)
    adjacency_matrix = img_to_graph(mask, mask)
    parent, children, height, adjacency_matrix = ward_tree(X.T,
                                                        adjacency_matrix)
    n_nodes = 2 * X.shape[1] - 1
    assert(len(parent) == n_nodes)
    assert(len(children) == n_nodes)
    assert(len(height) == n_nodes)


def test_unstructured_ward_tree():
    """
    Check that we obtain the correct solution for unstructured ward tree.
    """
    np.random.seed(0)
    X = np.random.randn(50, 100)
    parent, children, height, adjacency_matrix = ward_tree(X.T)
    n_nodes = 2 * X.shape[1] - 1
    assert(len(parent) == n_nodes)
    assert(len(children) == n_nodes)
    assert(len(height) == n_nodes)


def test_height_ward_tree():
    """
    Check that the height of ward tree is sorted.
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(50, 100)
    adjacency_matrix = img_to_graph(mask, mask)
    parent, children, height, adjacency_matrix = ward_tree(X.T,\
                                                adjacency_matrix)
    n_nodes = 2 * X.shape[1] - 1
    assert(np.sum(np.argsort(height)[100:] - np.arange(n_nodes)[100:]) == 0)


def test_ward_clustering():
    """
    Check that we obtain the correct number of clusters with Ward clustering.
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(50, 100)
    adjacency_matrix = img_to_graph(mask, mask)
    clustering = Ward(10)
    clustering.fit(X.T, adjacency_matrix)
    assert(np.size(np.unique(clustering.labels_)) == 10)

    Xred = clustering.transform(X.T)
    Xfull = clustering.inverse_transform(Xred)
    assert(np.unique(Xfull[0]).size == 10)
    Xfull = clustering.inverse_transform(Xred[:,0])
    assert(np.unique(Xfull).size == 10)

if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
