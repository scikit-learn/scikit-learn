"""
Several basic tests for hierarchical clustering procedures

Author : Vincent Michel, 2010
"""

import numpy as np
from scipy.cluster import hierarchy

from sklearn.cluster import Ward, WardAgglomeration, ward_tree
from sklearn.cluster.hierarchical import _hc_cut
from sklearn.feature_extraction.image import grid_to_graph


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

def test_connectivity_popagation():
    """
    Check that connectivity in the ward tree is propagated correctly during merging.  
    """
    from sklearn.neighbors import kneighbors_graph
    
    X = np.array([(0.0140, 0.1020), (0.0141, 0.0997), (0.0142, 0.0967),
                  (0.0142, 0.0951), (0.0141, 0.0954), (0.0139, 0.0961),
                  (0.0134, 0.0974), (0.0124, 0.0991), (0.0111, 0.1012),
                  (0.0094, 0.1035), (0.0074, 0.1052), (0.0057, 0.1044),
                  (0.0039, 0.1033), (0.0018, 0.1017), (-0.0005, 0.0996),
                  (-0.0030, 0.0973), (-0.0055, 0.0951), (-0.0076, 0.0952),
                  (-0.0094, 0.0959), (-0.0106, 0.0975), (-0.0111, 0.0999), 
                  (-0.0107, 0.1034), (-0.0094, 0.1048), (-0.0076, 0.1047),
                  (-0.0053, 0.1053), (-0.0028, 0.1053), (-0.0001, 0.1053),
                  (0.0026, 0.1053), (0.0055, 0.1053), (0.0086, 0.1053), 
                  (0.0121, 0.1053), (0.0161, 0.1053), (0.0164, 0.1051), 
                  (0.0158, 0.1053), (0.0161, 0.1053), (0.0168, 0.1053), 
                  (0.0167, 0.1053), (0.0167, 0.1053), (0.0167, 0.1053), 
                  (0.0167, 0.1053), (0.0167, 0.1053), (0.0167, 0.1053), 
                  (0.0167, 0.1053), (0.0167, 0.1053), (0.0168, 0.1053), 
                  (0.0168, 0.1053), (0.0168, 0.1053), (0.0168, 0.1053), 
                  (0.0168, 0.1053), (0.0168, 0.1053), (0.0168, 0.1053), 
                  (0.0168, 0.1052), (0.0168, 0.1049), (0.0168, 0.1044), 
                  (0.0168, 0.1037), (0.0168, 0.1027), (0.0168, 0.1012), 
                  (0.0168, 0.0994), (0.0168, 0.0971), (0.0168, 0.0951), 
                  (0.0168, 0.0951), (0.0168, 0.0952), (0.0168, 0.0955),
                  (0.0168, 0.0964), (0.0168, 0.0979), (0.0165, 0.1002),
                  (0.0159, 0.1032), (0.0154, 0.1049), (0.0147, 0.1046), 
                  (0.0138, 0.1051), (0.0126, 0.1053), (0.0112, 0.1054), 
                  (0.0097, 0.1053), (0.0079, 0.1053), (0.0060, 0.1053), 
                  (0.0040, 0.1053), (0.0020, 0.1053), (-0.0002, 0.1053),
                  (-0.0026, 0.1053), (-0.0051, 0.1053), (-0.0077, 0.1053),
                  (-0.0104, 0.1053), (-0.0132, 0.1053), (-0.0160, 0.1053),
                  (-0.0187, 0.1053), (-0.0212, 0.1054), (-0.0239, 0.1053),
                  (-0.0267, 0.1054), (-0.0296, 0.1054), (-0.0327, 0.1054),
                  (-0.0360, 0.1054), (-0.0398, 0.1054), (-0.0441, 0.1054),
                  (-0.0486, 0.1054), (-0.0533, 0.1054), (-0.0584, 0.1054), 
                  (-0.0601, 0.1053), (-0.0602, 0.1054), (-0.0600, 0.1053),
                  (-0.0596, 0.1054)])
    
    connectivity = kneighbors_graph(X, n_neighbors=10)
    ward = Ward(n_clusters=4, connectivity=connectivity)
    ward.fit(X) # If changes are not propagated correctly, fit crashes with an Exception

if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
