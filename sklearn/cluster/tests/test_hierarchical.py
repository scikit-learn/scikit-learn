"""
Several basic tests for hierarchical clustering procedures

Author : Vincent Michel, 2010
"""

import numpy as np
from scipy.cluster import hierarchy
from nose.tools import assert_true

from sklearn.cluster import Ward, WardAgglomeration, ward_tree
from sklearn.cluster.hierarchical import _hc_cut
from sklearn.feature_extraction.image import grid_to_graph


def test_structured_ward_tree():
    """
    Check that we obtain the correct solution for structured ward tree.
    """
    rnd = np.random.RandomState(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = rnd.randn(50, 100)
    connectivity = grid_to_graph(*mask.shape)
    children, n_components, n_leaves = ward_tree(X.T, connectivity)
    n_nodes = 2 * X.shape[1] - 1
    assert_true(len(children) + n_leaves == n_nodes)


def test_unstructured_ward_tree():
    """
    Check that we obtain the correct solution for unstructured ward tree.
    """
    rnd = np.random.RandomState(0)
    X = rnd.randn(50, 100)
    children, n_nodes, n_leaves = ward_tree(X.T)
    n_nodes = 2 * X.shape[1] - 1
    assert_true(len(children) + n_leaves == n_nodes)


def test_height_ward_tree():
    """
    Check that the height of ward tree is sorted.
    """
    rnd = np.random.RandomState(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = rnd.randn(50, 100)
    connectivity = grid_to_graph(*mask.shape)
    children, n_nodes, n_leaves = ward_tree(X.T, connectivity)
    n_nodes = 2 * X.shape[1] - 1
    assert_true(len(children) + n_leaves == n_nodes)


def test_ward_clustering():
    """
    Check that we obtain the correct number of clusters with Ward clustering.
    """
    rnd = np.random.RandomState(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = rnd.randn(100, 50)
    connectivity = grid_to_graph(*mask.shape)
    clustering = Ward(n_clusters=10, connectivity=connectivity)
    clustering.fit(X)
    assert_true(np.size(np.unique(clustering.labels_)) == 10)


def test_ward_agglomeration():
    """
    Check that we obtain the correct solution in a simplistic case
    """
    rnd = np.random.RandomState(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = rnd.randn(50, 100)
    connectivity = grid_to_graph(*mask.shape)
    ward = WardAgglomeration(n_clusters=5, connectivity=connectivity)
    ward.fit(X)
    assert_true(np.size(np.unique(ward.labels_)) == 5)

    Xred = ward.transform(X)
    assert_true(Xred.shape[1] == 5)
    Xfull = ward.inverse_transform(Xred)
    assert_true(np.unique(Xfull[0]).size == 5)


def assess_same_labelling(cut1, cut2):
    """Util for comparison with scipy"""
    co_clust = []
    for cut in [cut1, cut2]:
        n = len(cut)
        k = cut.max() + 1
        ecut = np.zeros((n, k))
        ecut[np.arange(n), cut] = 1
        co_clust.append(np.dot(ecut, ecut.T))
    assert_true((co_clust[0] == co_clust[1]).all())


def test_scikit_vs_scipy():
    """Test scikit ward with full connectivity (i.e. unstructured) vs scipy
    """
    from scipy.sparse import lil_matrix
    n, p, k = 10, 5, 3
    rnd = np.random.RandomState(0)

    connectivity = lil_matrix(np.ones((n, n)))
    for i in range(5):
        X = .1 * rnd.normal(size=(n, p))
        X -= 4 * np.arange(n)[:, np.newaxis]
        X -= X.mean(axis=1)[:, np.newaxis]

        out = hierarchy.ward(X)

        children_ = out[:, :2].astype(np.int)
        children, _, n_leaves = ward_tree(X, connectivity)

        cut = _hc_cut(k, children, n_leaves)
        cut_ = _hc_cut(k, children_, n_leaves)
        assess_same_labelling(cut, cut_)


def test_connectivity_popagation():
    """
    Check that connectivity in the ward tree is propagated correctly during
    merging.
    """
    from sklearn.neighbors import NearestNeighbors

    X = np.array([(.014, .120), (.014, .099), (.014, .097),
                  (.017, .153), (.017, .153), (.018, .153),
                  (.018, .153), (.018, .153), (.018, .153),
                  (.018, .153), (.018, .153), (.018, .153),
                  (.018, .152), (.018, .149), (.018, .144),
                 ])
    nn = NearestNeighbors(n_neighbors=10, warn_on_equidistant=False).fit(X)
    connectivity = nn.kneighbors_graph(X)
    ward = Ward(n_clusters=4, connectivity=connectivity)
    # If changes are not propagated correctly, fit crashes with an
    # IndexError
    ward.fit(X)


def test_connectivity_fixing_non_lil():
    """
    Check non regression of a bug if a non item assignable connectivity is
    provided with more than one component.
    """
    # create dummy data
    x = np.array([[0, 0], [1, 1]])
    # create a mask with several components to force connectivity fixing
    m = np.array([[True, False], [False, True]])
    c = grid_to_graph(n_x=2, n_y=2, mask=m)
    w = Ward(connectivity=c)
    w.fit(x)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
