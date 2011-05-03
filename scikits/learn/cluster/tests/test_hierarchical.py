"""
Several basic tests for hierarchical clustering procedures

Author : Vincent Michel, 2010
"""

import numpy as np
from scipy.cluster import hierarchy

from scikits.learn.cluster import Ward, WardAgglomeration, ward_tree
from scikits.learn.cluster.hierarchical import _hc_cut
from scikits.learn.feature_extraction.image import grid_to_graph


def test_structured_ward_tree():
    """
    Check that we obtain the correct solution for structured ward tree.
    """
    np.random.seed(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = np.random.randn(50, 100)
    connectivity = grid_to_graph(*mask.shape)
    children, n_components, n_leaves = ward_tree(X.T, connectivity)
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
    connectivity = grid_to_graph(*mask.shape)
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
    connectivity = grid_to_graph(*mask.shape)
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
    connectivity = grid_to_graph(*mask.shape)
    ward = WardAgglomeration(n_clusters=5, connectivity=connectivity)
    ward.fit(X)
    assert(np.size(np.unique(ward.labels_)) == 5)

    Xred = ward.transform(X)
    assert(Xred.shape[1] == 5)
    Xfull = ward.inverse_transform(Xred)
    assert(np.unique(Xfull[0]).size == 5)


def assess_same_labelling(cut1, cut2):
    """Util for comparison with scipy"""
    co_clust = []
    for cut in [cut1, cut2]:
        n = len(cut)
        k = cut.max() + 1
        ecut = np.zeros((n, k))
        ecut[np.arange(n), cut] = 1
        co_clust.append(np.dot(ecut, ecut.T))
    assert((co_clust[0] == co_clust[1]).all())


def test_scikit_vs_scipy():
    """Test scikit ward with full connectivity (i.e. unstructured) vs scipy
    """
    from scipy.sparse import lil_matrix
    n, p, k = 10, 5, 3

    connectivity = lil_matrix(np.ones((n, n)))
    for i in range(5):
        X = .1*np.random.normal(size=(n, p))
        X -= 4*np.arange(n)[:, np.newaxis]
        X -= X.mean(axis=1)[:, np.newaxis]

        out = hierarchy.ward(X)

        children_ = out[:, :2].astype(np.int)
        children, _, n_leaves = ward_tree(X, connectivity)

        cut = _hc_cut(k, children, n_leaves)
        cut_ = _hc_cut(k, children_, n_leaves)
        assess_same_labelling(cut, cut_)

if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
