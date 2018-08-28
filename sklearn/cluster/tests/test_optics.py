# Authors: Shane Grigsby <refuge@rocktalus.com>
#          Amy X. Zhang <axz@mit.edu>
# License: BSD 3 clause

import numpy as np
import pytest

from sklearn.datasets.samples_generator import make_blobs
from sklearn.cluster.optics_ import OPTICS
from sklearn.cluster.optics_ import _TreeNode, _cluster_tree
from sklearn.cluster.optics_ import _find_local_maxima
from sklearn.metrics.cluster import contingency_matrix
from sklearn.cluster.dbscan_ import DBSCAN
from sklearn.utils.testing import assert_equal, assert_warns
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_allclose
from sklearn.utils import _IS_32BIT

from sklearn.cluster.tests.common import generate_clustered_data


def test_correct_number_of_clusters():
    # in 'auto' mode

    n_clusters = 3
    X = generate_clustered_data(n_clusters=n_clusters)
    # Parameters chosen specifically for this task.
    # Compute OPTICS
    clust = OPTICS(max_bound=5.0 * 6.0, min_samples=4, metric='euclidean')
    clust.fit(X)
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(clust.labels_)) - int(-1 in clust.labels_)
    assert_equal(n_clusters_1, n_clusters)


def test_minimum_number_of_sample_check():
    # test that we check a minimum number of samples
    msg = ("Number of training samples (n_samples=1) must be greater than "
           "min_samples (min_samples=10) used for clustering.")

    # Compute OPTICS
    X = [[1, 1]]
    clust = OPTICS(max_bound=5.0 * 0.3, min_samples=10)

    # Run the fit
    assert_raise_message(ValueError, msg, clust.fit, X)


def test_empty_extract():
    # Test extract where fit() has not yet been run.
    msg = ("This OPTICS instance is not fitted yet. Call 'fit' with "
           "appropriate arguments before using this method.")
    clust = OPTICS(max_bound=5.0 * 0.3, min_samples=10)
    assert_raise_message(ValueError, msg, clust.extract_dbscan, 0.01)


def test_bad_extract():
    # Test an extraction of eps too close to original eps
    msg = "Specify an epsilon smaller than 0.015. Got 0.3."
    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers,
                                cluster_std=0.4, random_state=0)

    # Compute OPTICS
    clust = OPTICS(max_bound=5.0 * 0.003, min_samples=10)
    clust2 = clust.fit(X)
    assert_raise_message(ValueError, msg, clust2.extract_dbscan, 0.3)


def test_close_extract():
    # Test extract where extraction eps is close to scaled epsPrime

    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers,
                                cluster_std=0.4, random_state=0)

    # Compute OPTICS
    clust = OPTICS(max_bound=1.0, min_samples=10)
    clust3 = clust.fit(X)
    # check warning when centers are passed
    assert_warns(RuntimeWarning, clust3.extract_dbscan, .3)
    # Cluster ordering starts at 0; max cluster label = 2 is 3 clusters
    assert_equal(max(clust3.extract_dbscan(.3)[1]), 2)


@pytest.mark.parametrize('eps', [0.1, .3, .5])
@pytest.mark.parametrize('min_samples', [3, 10, 20])
def test_dbscan_optics_parity(eps, min_samples):
    # Test that OPTICS clustering labels are <= 5% difference of DBSCAN

    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers,
                                cluster_std=0.4, random_state=0)

    # calculate optics with dbscan extract at 0.3 epsilon
    op = OPTICS(min_samples=min_samples).fit(X)
    core_optics, labels_optics = op.extract_dbscan(eps)

    # calculate dbscan labels
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)

    contingency = contingency_matrix(db.labels_, labels_optics)
    agree = min(np.sum(np.max(contingency, axis=0)),
                np.sum(np.max(contingency, axis=1)))
    disagree = X.shape[0] - agree

    # verify core_labels match
    assert_array_equal(core_optics, db.core_sample_indices_)

    non_core_count = len(labels_optics) - len(core_optics)
    percent_mismatch = np.round((disagree - 1) / non_core_count, 2)

    # verify label mismatch is <= 5% labels
    assert percent_mismatch <= 0.05


def test_auto_extract_hier():
    # Tests auto extraction gets correct # of clusters with varying density

    # Generate sample data
    rng = np.random.RandomState(0)
    n_points_per_cluster = 250

    C1 = [-5, -2] + .8 * rng.randn(n_points_per_cluster, 2)
    C2 = [4, -1] + .1 * rng.randn(n_points_per_cluster, 2)
    C3 = [1, -2] + .2 * rng.randn(n_points_per_cluster, 2)
    C4 = [-2, 3] + .3 * rng.randn(n_points_per_cluster, 2)
    C5 = [3, -2] + 1.6 * rng.randn(n_points_per_cluster, 2)
    C6 = [5, 6] + 2 * rng.randn(n_points_per_cluster, 2)
    X = np.vstack((C1, C2, C3, C4, C5, C6))

    # Compute OPTICS

    clust = OPTICS(min_samples=9)

    # Run the fit
    clust.fit(X)

    assert_equal(len(set(clust.labels_)), 6)


@pytest.mark.parametrize("reach, n_child, members", [
    (np.array([np.inf, 0.9, 0.9, 1.0, 0.89, 0.88, 10, .9, .9, .9, 10, 0.9,
               0.9, 0.89, 0.88, 10, .9, .9, .9, .9]), 2, np.r_[0:6]),
    (np.array([np.inf, 0.9, 0.9, 0.9, 0.89, 0.88, 10, .9, .9, .9, 10, 0.9,
               0.9, 0.89, 0.88, 100, .9, .9, .9, .9]), 1, np.r_[0:15])])
def test_cluster_sigmin_pruning(reach, n_child, members):
    # Tests pruning left and right, insignificant splitpoints, empty nodelists
    # Parameters chosen specifically for this task

    # Case 1: Three pseudo clusters, 2 of which are too small
    # Case 2: Two pseudo clusters, 1 of which are too small
    # Normalize
    reach = reach / np.max(reach[1:])

    ordering = np.r_[0:20]
    cluster_boundaries = _find_local_maxima(reach, 5)
    root = _TreeNode(ordering, 0, 20, None)

    # Build cluster tree inplace on root node
    _cluster_tree(root, None, cluster_boundaries, reach, ordering,
                  5, .75, .7, .4, .3)
    assert_equal(root.split_point, cluster_boundaries[0])
    assert_equal(n_child, len(root.children))
    assert_array_equal(members, root.children[0].points)


def test_reach_dists():
    # Tests against known extraction array

    rng = np.random.RandomState(0)
    n_points_per_cluster = 50

    C1 = [-5, -2] + 8 * rng.randn(n_points_per_cluster, 2)
    C2 = [4, -1] + 1 * rng.randn(n_points_per_cluster, 2)
    C3 = [1, -2] + 2 * rng.randn(n_points_per_cluster, 2)
    C4 = [-2, 3] + 3 * rng.randn(n_points_per_cluster, 2)
    C5 = [3, -2] + 16 * rng.randn(n_points_per_cluster, 2)
    C6 = [5, 6] + 20 * rng.randn(n_points_per_cluster, 2)
    X = np.vstack((C1, C2, C3, C4, C5, C6))

    # Compute OPTICS

    clust = OPTICS(min_samples=10, metric='minkowski')

    # Run the fit
    clust.fit(X)

    # Expected values, matches 'RD' results from:
    # http://chemometria.us.edu.pl/download/optics.py

    v = [np.inf,  7.59885069,  6.17890089,  1.05945707,  1.30847269,
         3.69580892,  0.79365587,  1.48862487,  1.67837779,  3.14568804,
         7.53854141,  2.72790669,  8.24245162,  2.0378668,  6.23092657,
         1.41644975,  6.82449182,  2.3740666,  3.82029265,  3.1659918,
         5.41312893,  7.41189078,  3.6499356,  4.65838153,  4.7790213,
         3.77505395,  5.41312893,  1.24740706,  1.24740706,  3.76568821,
         3.76568821,  5.94047022,  2.50040119,  5.44241104,  3.76568821,
         0.79365587,  4.51358931,  2.79521045,  3.76568821,  2.3740666,
         4.36094583,  5.3200597,  6.23092657,  1.1851792,  4.36094583,
         3.56833046,  1.95885905,  1.60343068,  4.35506869,  1.30644367,
         1.03067612,  0.58566183,  0.8509913,  0.78701264,  1.75672123,
         1.647896,  0.76893526,  0.59517811,  0.67697899,  0.62607998,
         0.79112254,  0.80706656,  0.76240524,  0.73744965,  0.67697899,
         0.66665535,  0.6786622,  0.59517811,  0.67697899,  0.79112254,
         0.58566183,  0.62607998,  4.12425707,  0.59517811,  0.76240524,
         0.83682819,  0.78701264,  0.81416946,  0.79112254,  0.76240524,
         0.79516239,  0.67697899,  0.59517811,  0.76240524,  1.29508353,
         0.83412449,  0.75708599,  0.59517811,  0.76240524,  0.66665535,
         0.76240524,  0.89652267,  0.90300374,  0.59517811,  0.73744965,
         0.58566183,  0.78384875,  0.67697899,  0.62607998,  1.79315246,
         0.88200422,  0.66665535,  2.0482758,  1.05562986,  0.97287579,
         0.76240524,  1.05945707,  0.88200422,  0.79365587,  0.79112254,
         0.88200422,  0.76303359,  1.30876907,  0.7629389,  0.88200422,
         1.19273656,  0.88200422,  1.42757935,  2.58860441,  0.72822051,
         0.88200422,  0.88200422,  1.96044868,  0.58566183,  1.99021141,
         1.13371735,  1.03469971,  0.61358788,  1.81518036,  1.30876907,
         1.0976403,  0.88200422,  1.30876907,  0.80936808,  0.91629931,
         2.64046804,  0.67697899,  0.88200422,  0.97287579,  1.25072302,
         0.79365587,  1.36943228,  0.82840271,  1.48862487,  0.72439367,
         0.88200422,  1.68717255,  0.79365587,  0.79365587,  0.7629389,
         2.80341016,  1.5853267,  1.48862487,  1.95007681,  1.52146026,
         1.83591272,  1.2600627,  1.25072302,  1.25072302,  1.25072302,
         1.24740706,  1.48862487,  1.24740706,  3.53766973,  1.95885905,
         0.7629389,  0.7629389,  3.56833046,  2.34116433,  2.53639304,
         1.282933,  1.2600627,  1.24740706,  1.33097173,  1.56632179,
         1.66211533,  0.91844338,  1.25852424,  1.83912019,  1.60343068,
         1.56632179,  1.66011024,  1.60343068,  1.25072302,  1.5592943,
         1.48862487,  1.48862487,  1.64533938,  0.88200422,  1.83591272,
         2.34116433,  2.79897124,  1.21087988,  1.69485371,  2.07750488,
         1.48862487,  1.30876907,  1.41644975,  2.41362097,  2.22219201,
         6.82449182,  5.86214331,  7.52179844,  6.20734688,  5.3200597,
         7.41189078,  1.99021141, 10.12563499,  4.5948987, 12.493863,
         6.20734688, 11.08001958,  0.96195202, 18.62694837,  1.30847269,
         7.0559566,  6.20734688,  4.36094583,  6.20734688,  6.20734688,
         6.82449182,  0.93105834,  6.27263296,  8.24245162,  3.1659918,
         5.52143272,  4.12425707, 10.79984612, 12.55260431,  1.14228938,
         6.35782964,  1.30644367, 13.86060527,  5.94047022, 11.08001958,
         9.96891342,  4.65113924, 12.3547993, 11.84255723,  6.20734688,
         8.53664632,  4.1253444, 16.82829994,  8.1786272,  0.91844338,
         4.1253444,  9.57370251, 17.84963742,  4.96798516,  3.80240412,
         5.42876951,  8.1786272,  1.99021141, 14.81182123, 15.9525836,
         12.493863, 18.22386404,  3.92572449,  9.57370251, 14.54063695,
         4.53537202, 10.06701936,  5.90471515, 12.493863, 21.99435576,
         15.14654311,  8.54895489,  7.36169303,  1.48862487,  1.56632179,
         8.13436918, 13.86060527,  1.30876907, 19.29062527,  8.71319171,
         5.86214331,  7.17085401, 13.86060527, 13.41032703,  7.49252559,
         15.00240034, 11.46569783,  4.15139784,  8.27349636,  7.59885069,
         13.32055744,  8.27349636,  6.20734688,  5.94047022, 10.79984612,
         12.24758945, 11.14464297,  2.29458767,  6.01813273, 20.81401082,
         8.56541443,  5.3200597,  6.9362395,  1.25852424,  4.36094583]

    # FIXME: known failure in 32bit Linux; numerical imprecision results in
    # different ordering in quick_scan
    if _IS_32BIT:  # pragma: no cover
        assert_allclose(clust.reachability_, np.array(v), rtol=1e-2)
    else:
        # we compare to truncated decimals, so use atol
        assert_allclose(clust.reachability_, np.array(v), atol=1e-5)
