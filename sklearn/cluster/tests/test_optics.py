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

    v = [   np.inf, 7.5988507, 6.1789009, 1.0594571, 1.3084727, 3.6958089,
            0.7936559, 1.4886249, 1.6783778, 3.1456880, 7.5385414, 2.7279067,
            8.2424516, 2.0378668, 6.2309266, 1.4164497, 6.8244918, 2.3740666,
            3.8202927, 3.1659918, 5.4131289, 7.4118908, 3.6499356, 4.6583815,
            4.7790213, 3.7750540, 5.4131289, 1.2474071, 1.2474071, 3.7656882,
            3.7656882, 5.9404702, 2.5004012, 5.4424110, 3.7656882, 0.7936559,
            4.5135893, 2.7952104, 3.7656882, 2.3740666, 4.3609458, 5.3200597,
            6.2309266, 1.1851792, 4.3609458, 3.5683305, 1.9588590, 1.6034307,
            4.3550687, 1.3064437, 1.0306761, 0.5856618, 0.8509913, 0.7870126,
            1.7567212, 1.6478960, 0.7689353, 0.5951781, 0.6769790, 0.6260800,
            0.7911225, 0.8070666, 0.7624052, 0.7374497, 0.6769790, 0.6666553,
            0.6786622, 0.5951781, 0.6769790, 0.7911225, 0.5856618, 0.6260800,
            4.1242571, 0.5951781, 0.7624052, 0.8368282, 0.7870126, 0.8141695,
            0.7911225, 0.7624052, 0.7951624, 0.6769790, 0.5951781, 0.7624052,
            1.2950835, 0.8341245, 0.7570860, 0.5951781, 0.7624052, 0.6666553,
            0.7624052, 0.8965227, 0.9030037, 0.5951781, 0.7374497, 0.5856618,
            0.7838487, 0.6769790, 0.6260800, 1.7931525, 0.8820042, 0.6666553,
            2.0482758, 1.0556299, 0.9728758, 0.7624052, 1.0594571, 0.8820042,
            0.7936559, 0.7911225, 0.8820042, 0.7630336, 1.3087691, 0.7629389,
            0.8820042, 1.1927366, 0.8820042, 1.4275793, 2.5886044, 0.7282205,
            0.8820042, 0.8820042, 1.9604487, 0.5856618, 1.9902114, 1.1337173,
            1.0346997, 0.6135879, 1.8151804, 1.3087691, 1.0976403, 0.8820042,
            1.3087691, 0.8093681, 0.9162993, 2.6404680, 0.6769790, 0.8820042,
            0.9728758, 1.2507230, 0.7936559, 1.3694323, 0.8284027, 1.4886249,
            0.7243937, 0.8820042, 1.6871725, 0.7936559, 0.7936559, 0.7629389,
            2.8034102, 1.5853267, 1.4886249, 1.9500768, 1.5214603, 1.8359127,
            1.2600627, 1.2507230, 1.2507230, 1.2507230, 1.2474071, 1.4886249,
            1.2474071, 3.5376697, 1.9588590, 0.7629389, 0.7629389, 3.5683305,
            2.3411643, 2.5363930, 1.2829330, 1.2600627, 1.2474071, 1.3309717,
            1.5663218, 1.6621153, 0.9184434, 1.2585242, 1.8391202, 1.6034307,
            1.5663218, 1.6601102, 1.6034307, 1.2507230, 1.5592943, 1.4886249,
            1.4886249, 1.6453394, 0.8820042, 1.8359127, 2.3411643, 2.7989712,
            1.2108799, 1.6948537, 2.0775049, 1.4886249, 1.3087691, 1.4164497,
            2.4136210, 2.2221920, 6.8244918, 5.8621433, 7.5217984, 6.2073469,
            5.3200597, 7.4118908, 1.9902114,10.1256350, 4.5948987,12.4938630,
            6.2073469,11.0800196, 0.9619520,18.6269484, 1.3084727, 7.0559566,
            6.2073469, 4.3609458, 6.2073469, 6.2073469, 6.8244918, 0.9310583,
            6.2726330, 8.2424516, 3.1659918, 5.5214327, 4.1242571, 0.7998461,
            12.5526043, 1.1422894, 6.3578296, 1.3064437,13.8606053, 5.9404702,
            11.0800196, 9.9689134, 4.6511392,12.3547993,11.8425572, 6.2073469,
            8.5366463, 4.1253444,16.8282999, 8.1786272, 0.9184434, 4.1253444,
            9.5737025,17.8496374, 4.9679852, 3.8024041, 5.4287695, 8.1786272,
            1.9902114,14.8118212,15.9525836,12.4938630,18.2238640, 3.9257245,
            9.5737025,14.5406370, 4.5353720,10.0670194, 5.9047151,12.4938630,
            21.9943558,15.1465431, 8.5489549, 7.3616930, 1.4886249, 1.5663218,
            8.1343692,13.8606053, 1.3087691,19.2906253, 8.7131917, 5.8621433,
            7.1708540,13.8606053,13.4103270, 7.4925256,15.0024003,11.4656978,
            4.1513978, 8.2734964, 7.5988507,13.3205574, 8.2734964, 6.2073469,
            5.9404702,10.7998461,12.2475895,11.1446430, 2.2945877, 6.0181327,
            20.8140108, 8.5654144, 5.3200597, 6.9362395, 1.2585242, 4.3609458]

    # FIXME: known failure in 32bit Linux; numerical imprecision results in
    # different ordering in quick_scan
    if _IS_32BIT:  # pragma: no cover
        assert_allclose(clust.reachability_, np.array(v), rtol=1e-2)
    else:
        # we compare to truncated decimals, so use atol
        assert_allclose(clust.reachability_, np.array(v), atol=1e-5)
