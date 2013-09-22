"""
Several basic tests for hierarchical clustering procedures

"""
# Authors: Vincent Michel, 2010, Gael Varoquaux 2012
# License: BSD 3 clause
import warnings
from tempfile import mkdtemp

import numpy as np
from scipy import sparse
from scipy.cluster import hierarchy

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_almost_equal

from sklearn.cluster import Ward, WardAgglomeration, ward_tree
from sklearn.cluster import AgglomerativeClustering, FeatureAgglomeration
from sklearn.cluster.hierarchical import (_hc_cut, _TREE_BUILDERS,
                                          linkage_tree)
from sklearn.feature_extraction.image import grid_to_graph
from sklearn.metrics.pairwise import PAIRED_DISTANCES
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.cluster._hierarchical import average_merge, max_merge
from sklearn.utils.fast_dict import IntFloatDict


def test_linkage_misc():
    # Misc tests on linkage
    X = np.ones((5, 5))
    assert_raises(ValueError,
                  AgglomerativeClustering(linkage='foobar').fit,
                  X)
    assert_raises(ValueError, linkage_tree, X, linkage='foobar')
    assert_raises(ValueError, linkage_tree, X, connectivity=np.ones((4, 4)))

    # Smoke test FeatureAgglomeration
    FeatureAgglomeration().fit(X)

    with warnings.catch_warnings(record=True) as warning_list:
        warnings.simplefilter("always", UserWarning)
        # Use the copy argument, to raise a warning
        Ward(copy=True).fit(X)
    # We should be getting 2 warnings: one for using Ward that is
    # deprecated, one for using the copy argument
    assert_equal(len(warning_list), 2)

    with warnings.catch_warnings(record=True) as warning_list:
        warnings.simplefilter("always", UserWarning)
        # Use the copy argument, to raise a warning
        ward_tree(X, copy=True)
    # We should be getting 1 warnings: for using the copy argument
    assert_equal(len(warning_list), 1)


def test_structured_linkage_tree():
    """
    Check that we obtain the correct solution for structured linkage trees.
    """
    rnd = np.random.RandomState(0)
    mask = np.ones([10, 10], dtype=np.bool)
    # Avoiding a mask with only 'True' entries
    mask[4:7, 4:7] = 0
    X = rnd.randn(50, 100)
    connectivity = grid_to_graph(*mask.shape)
    for tree_builder in _TREE_BUILDERS.values():
        children, n_components, n_leaves, parent = \
            tree_builder(X.T, connectivity)
        n_nodes = 2 * X.shape[1] - 1
        assert_true(len(children) + n_leaves == n_nodes)
        # Check that ward_tree raises a ValueError with a connectivity matrix
        # of the wrong shape
        assert_raises(ValueError,
                      tree_builder, X.T, np.ones((4, 4)))


def test_unstructured_linkage_tree():
    """
    Check that we obtain the correct solution for unstructured linkage trees.
    """
    rnd = np.random.RandomState(0)
    X = rnd.randn(50, 100)

    for tree_builder in _TREE_BUILDERS.values():
        for this_X in (X, X[0]):
            with warnings.catch_warnings(record=True) as warning_list:
                warnings.simplefilter("always", UserWarning)
                warnings.simplefilter("ignore", DeprecationWarning)

                # With specified a number of clusters just for the sake of
                # raising a warning and testing the warning code
                children, n_nodes, n_leaves, parent = tree_builder(
                    this_X.T, n_clusters=10)
            assert_equal(len(warning_list), 1)
            n_nodes = 2 * X.shape[1] - 1
            assert_equal(len(children) + n_leaves, n_nodes)


def test_height_linkage_tree():
    """
    Check that the height of the results of linkage tree is sorted.
    """
    rnd = np.random.RandomState(0)
    mask = np.ones([10, 10], dtype=np.bool)
    X = rnd.randn(50, 100)
    connectivity = grid_to_graph(*mask.shape)
    for linkage_func in _TREE_BUILDERS.values():
        children, n_nodes, n_leaves, parent = linkage_func(X.T, connectivity)
        n_nodes = 2 * X.shape[1] - 1
        assert_true(len(children) + n_leaves == n_nodes)


def test_agglomerative_clustering():
    """
    Check that we obtain the correct number of clusters with
    agglomerative clustering.
    """
    rnd = np.random.RandomState(0)
    mask = np.ones([10, 10], dtype=np.bool)
    n_samples = 100
    X = rnd.randn(n_samples, 50)
    connectivity = grid_to_graph(*mask.shape)
    for linkage in ("ward", "complete", "average"):
        clustering = AgglomerativeClustering(n_clusters=10,
                                             connectivity=connectivity,
                                             linkage=linkage)
        clustering.fit(X)
        # test caching
        clustering = AgglomerativeClustering(
            n_clusters=10, connectivity=connectivity,
            memory=mkdtemp(),
            linkage=linkage)
        clustering.fit(X)
        labels = clustering.labels_
        assert_true(np.size(np.unique(labels)) == 10)
        # Turn caching off now
        clustering = AgglomerativeClustering(n_clusters=10,
                            connectivity=connectivity, linkage=linkage)
        # Check that we obtain the same solution with early-stopping of the
        # tree building
        clustering.compute_full_tree = False
        clustering.fit(X)
        np.testing.assert_array_equal(clustering.labels_, labels)
        clustering.connectivity = None
        clustering.fit(X)
        assert_true(np.size(np.unique(clustering.labels_)) == 10)
        # Check that we raise a TypeError on dense matrices
        clustering = AgglomerativeClustering(
            n_clusters=10,
            connectivity=sparse.lil_matrix(
                connectivity.todense()[:10, :10]),
            linkage=linkage)
        assert_raises(ValueError, clustering.fit, X)

    # Test that using ward with another metric than euclidean raises an
    # exception
    clustering = AgglomerativeClustering(
        n_clusters=10,
        connectivity=connectivity.todense(),
        affinity="manhattan",
        linkage="ward")
    assert_raises(ValueError, clustering.fit, X)

    # Test using another metric than euclidean works with linkage complete
    for affinity in PAIRED_DISTANCES.keys():
        # Compare our (structured) implementation to scipy
        clustering = AgglomerativeClustering(
            n_clusters=10,
            connectivity=np.ones((n_samples, n_samples)),
            affinity=affinity,
            linkage="complete")
        clustering.fit(X)
        clustering2 = AgglomerativeClustering(
            n_clusters=10,
            connectivity=None,
            affinity=affinity,
            linkage="complete")
        clustering2.fit(X)
        assert_equal(normalized_mutual_info_score(
                        clustering2.labels_,
                        clustering.labels_), 1)



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
    assert_array_almost_equal(ward.transform(Xfull), Xred)


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
    """Test scikit linkage with full connectivity (i.e. unstructured) vs scipy
    """
    n, p, k = 10, 5, 3
    rnd = np.random.RandomState(0)

    # Not using a lil_matrix here, just to check that non sparse
    # matrices are well handled
    connectivity = np.ones((n, n))
    for linkage in _TREE_BUILDERS.keys():
        for i in range(5):
            X = .1 * rnd.normal(size=(n, p))
            X -= 4 * np.arange(n)[:, np.newaxis]
            X -= X.mean(axis=1)[:, np.newaxis]

            out = hierarchy.linkage(X, method=linkage)

            children_ = out[:, :2].astype(np.int)
            children, _, n_leaves, _ = _TREE_BUILDERS[linkage](X, connectivity)

            cut = _hc_cut(k, children, n_leaves)
            cut_ = _hc_cut(k, children_, n_leaves)
            assess_same_labelling(cut, cut_)

    # Test error management in _hc_cut
    assert_raises(ValueError, _hc_cut, n_leaves + 1, children, n_leaves)


def test_connectivity_popagation():
    """
    Check that connectivity in the ward tree is propagated correctly during
    merging.
    """
    from sklearn.neighbors import kneighbors_graph

    X = np.array([(.014, .120), (.014, .099), (.014, .097),
                  (.017, .153), (.017, .153), (.018, .153),
                  (.018, .153), (.018, .153), (.018, .153),
                  (.018, .153), (.018, .153), (.018, .153),
                  (.018, .152), (.018, .149), (.018, .144),
                  ])
    connectivity = kneighbors_graph(X, 10)
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
    with warnings.catch_warnings(record=True):
        w.fit(x)


def test_int_float_dict():
    rng = np.random.RandomState(0)
    keys = np.unique(rng.randint(100, size=10).astype(np.intp))
    values = rng.rand(len(keys))

    d = IntFloatDict(keys, values)
    for key, value in zip(keys, values):
        assert d[key] == value

    other_keys = np.arange(50).astype(np.intp)[::2]
    other_values = 0.5*np.ones(50)[::2]
    other = IntFloatDict(other_keys, other_values)
    # Complete smoke test
    max_merge(d, other, mask=np.ones(100, dtype=np.intp), n_a=1, n_b=1)
    average_merge(d, other, mask=np.ones(100, dtype=np.intp), n_a=1, n_b=1)



if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
