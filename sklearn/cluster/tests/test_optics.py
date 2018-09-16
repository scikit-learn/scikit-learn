# Authors: Shane Grigsby <refuge@rocktalus.com>
#          Amy X. Zhang <axz@mit.edu>
# License: BSD 3 clause

from __future__ import print_function, division
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


rng = np.random.RandomState(0)
n_points_per_cluster = 250
C1 = [-5, -2] + .8 * rng.randn(n_points_per_cluster, 2)
C2 = [4, -1] + .1 * rng.randn(n_points_per_cluster, 2)
C3 = [1, -2] + .2 * rng.randn(n_points_per_cluster, 2)
C4 = [-2, 3] + .3 * rng.randn(n_points_per_cluster, 2)
C5 = [3, -2] + 1.6 * rng.randn(n_points_per_cluster, 2)
C6 = [5, 6] + 2 * rng.randn(n_points_per_cluster, 2)
X = np.vstack((C1, C2, C3, C4, C5, C6))


def test_correct_number_of_clusters():
    # in 'auto' mode

    n_clusters = 3
    X = generate_clustered_data(n_clusters=n_clusters)
    # Parameters chosen specifically for this task.
    # Compute OPTICS
    clust = OPTICS(max_eps=5.0 * 6.0, min_samples=4, metric='euclidean')
    clust.fit(X)
    # number of clusters, ignoring noise if present
    n_clusters_1 = len(set(clust.labels_)) - int(-1 in clust.labels_)
    assert_equal(n_clusters_1, n_clusters)

    # check attribute types and sizes
    assert clust.core_sample_indices_.ndim == 1
    assert clust.core_sample_indices_.size > 0
    assert clust.core_sample_indices_.dtype.kind == 'i'

    assert clust.labels_.shape == (len(X),)
    assert clust.labels_.dtype.kind == 'i'

    assert clust.reachability_.shape == (len(X),)
    assert clust.reachability_.dtype.kind == 'f'

    assert clust.core_distances_.shape == (len(X),)
    assert clust.core_distances_.dtype.kind == 'f'

    assert clust.ordering_.shape == (len(X),)
    assert clust.ordering_.dtype.kind == 'i'
    assert set(clust.ordering_) == set(range(len(X)))


def test_minimum_number_of_sample_check():
    # test that we check a minimum number of samples
    msg = ("Number of training samples (n_samples=1) must be greater than "
           "min_samples (min_samples=10) used for clustering.")

    # Compute OPTICS
    X = [[1, 1]]
    clust = OPTICS(max_eps=5.0 * 0.3, min_samples=10)

    # Run the fit
    assert_raise_message(ValueError, msg, clust.fit, X)


def test_empty_extract():
    # Test extract where fit() has not yet been run.
    msg = ("This OPTICS instance is not fitted yet. Call 'fit' with "
           "appropriate arguments before using this method.")
    clust = OPTICS(max_eps=5.0 * 0.3, min_samples=10)
    assert_raise_message(ValueError, msg, clust.extract_dbscan, 0.01)


def test_bad_extract():
    # Test an extraction of eps too close to original eps
    msg = "Specify an epsilon smaller than 0.015. Got 0.3."
    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers,
                                cluster_std=0.4, random_state=0)

    # Compute OPTICS
    clust = OPTICS(max_eps=5.0 * 0.003, min_samples=10)
    clust2 = clust.fit(X)
    assert_raise_message(ValueError, msg, clust2.extract_dbscan, 0.3)


def test_close_extract():
    # Test extract where extraction eps is close to scaled epsPrime

    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers,
                                cluster_std=0.4, random_state=0)

    # Compute OPTICS
    clust = OPTICS(max_eps=1.0, min_samples=10)
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
    clust = OPTICS(min_samples=9).fit(X)
    assert_equal(len(set(clust.labels_)), 6)


# try arbitrary minimum sizes
@pytest.mark.parametrize('min_cluster_size', range(2, X.shape[0] // 10, 23))
def test_min_cluster_size(min_cluster_size):
    redX = X[::10]  # reduce for speed
    clust = OPTICS(min_samples=9, min_cluster_size=min_cluster_size).fit(redX)
    cluster_sizes = np.bincount(clust.labels_[clust.labels_ != -1])
    if cluster_sizes.size:
        assert min(cluster_sizes) >= min_cluster_size
    # check behaviour is the same when min_cluster_size is a fraction
    clust_frac = OPTICS(min_samples=9,
                        min_cluster_size=min_cluster_size / redX.shape[0])
    clust_frac.fit(redX)
    assert_array_equal(clust.labels_, clust_frac.labels_)


@pytest.mark.parametrize('min_cluster_size', [0, -1, 1.1, 2.2])
def test_min_cluster_size_invalid(min_cluster_size):
    clust = OPTICS(min_cluster_size=min_cluster_size)
    with pytest.raises(ValueError, match="must be a positive integer or a "):
        clust.fit(X)


def test_min_cluster_size_invalid2():
    clust = OPTICS(min_cluster_size=len(X) + 1)
    with pytest.raises(ValueError, match="must be no greater than the "):
        clust.fit(X)


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

    rng = np.random.RandomState(12)
    n_points_per_cluster = 50

    C1 = [-5, -2] + .8 * rng.randn(n_points_per_cluster, 2)
    C2 = [4, -1] + .1 * rng.randn(n_points_per_cluster, 2)
    C3 = [1, -2] + .2 * rng.randn(n_points_per_cluster, 2)
    C4 = [-2, 3] + .3 * rng.randn(n_points_per_cluster, 2)
    C5 = [3, -2] + 1.6 * rng.randn(n_points_per_cluster, 2)
    C6 = [5, 6] + 2 * rng.randn(n_points_per_cluster, 2)
    X = np.vstack((C1, C2, C3, C4, C5, C6))

    # Compute OPTICS

    clust = OPTICS(min_samples=2, metric='minkowski')

    # Run the fit
    clust.fit(X)

    # Expected values, matches 'RD' results from:
    # http://chemometria.us.edu.pl/download/optics.py

    v = [np.inf, 0.42960718, 0.11846289, 0.04075161, 0.66462963,
       0.14038844, 0.25011207, 0.27893129, 0.46997676, 0.4246995 ,
       0.07988421, 0.42518697, 0.44432568, 0.28499101, 0.64324688,
       0.18981377, 0.26916942, 0.25865072, 0.35161532, 0.10411231,
       0.06366346, 0.31535031, 0.30748184, 0.74869561, 0.8999724 ,
       0.2948563 , 0.72984042, 0.33044984, 0.25210377, 0.23342478,
       0.27769411, 0.13377085, 0.25293594, 0.42487955, 0.27249511,
       0.21687232, 0.45289033, 1.05978849, 0.1168313 , 0.18793158,
       0.56142578, 0.15393804, 0.61081356, 0.83740986, 0.22767581,
       0.54131684, 0.27399932, 0.22147788, 0.25752016, 0.13851306,
       0.04889065, 0.01675676, 0.03160531, 0.02890523, 0.11911402,
       0.08550431, 0.07215453, 0.03832064, 0.02563518, 0.01143894,
       0.02783358, 0.04300404, 0.01424015, 0.04820423, 0.01539277,
       0.02495634, 0.06118386, 0.0196823 , 0.06387437, 0.04539848,
       0.03416283, 0.26393258, 0.04899627, 0.03451833, 0.04782168,
       0.03494979, 0.04468399, 0.067822  , 0.03939023, 0.03901752,
       0.07799087, 0.05325595, 0.01695318, 0.02019025, 0.00842129,
       0.08857675, 0.04814282, 0.02333937, 0.0194488 , 0.00636206,
       0.02428171, 0.01232994, 0.04752247, 0.06648731, 0.01923359,
       0.06478506, 0.02050194, 0.02570519, 0.00626935, 0.04566003,
       0.07229984, 0.05972662, 0.12402862, 0.05786566, 0.12902458,
       0.19397117, 0.06481248, 0.19744853, 0.13050219, 0.07368615,
       0.10293114, 0.02518936, 0.04737045, 0.07272921, 0.10337971,
       0.05448018, 0.04560843, 0.37185355, 0.17882154, 0.06383366,
       0.09666484, 0.15841938, 0.02802419, 0.07680612, 0.05505145,
       0.04726204, 0.04298729, 0.14588457, 0.05204779, 0.08758364,
       0.03986742, 0.0404145 , 0.02376484, 0.04621751, 0.15894636,
       0.06739961, 0.14745807, 0.07519989, 0.10413549, 0.03118367,
       0.27884321, 0.12577852, 0.01946909, 0.07165066, 0.03484924,
       0.07072021, 0.06056985, 0.0289722 , 0.04860705, 0.27330833,
       0.19253633, 0.07473675, 0.09600461, 0.03133039, 0.06652879,
       0.06213692, 0.12765036, 0.06305507, 0.03230761, 0.07761318,
       0.10812441, 0.12604482, 3.38166132, 0.07954152, 0.06355475,
       0.0718342 , 0.09100107, 0.06351387, 0.03578337, 0.09761135,
       0.08706804, 0.1112927 , 0.15610375, 0.03129849, 0.14103905,
       0.05477851, 0.08592214, 0.0985648 , 0.38047783, 0.07056325,
       0.08912907, 0.07893207, 0.07048013, 0.12146806, 0.18014014,
       0.07106212, 0.04583092, 0.05484359, 0.13149053, 0.08295075,
       0.05190986, 0.04277681, 0.09795172, 0.05191173, 0.06240285,
       0.08663622, 0.0387884 , 0.0314522 , 0.08132232, 0.10009446,
       0.40594025, 0.6958039 , 1.02272736, 1.01009808, 0.36132809,
       1.40000558, 0.50876207, 0.35528016, 0.41655233, 1.54185878,
       0.58283302, 1.08161647, 0.79083769, 0.21787959, 0.42869318,
       0.43233582, 0.77078842, 0.42731442, 0.43857058, 0.43947967,
       0.44450538, 0.21250175, 1.3021191 , 0.07902971, 0.29901135,
       0.43341179, 0.48832866, 1.50014159, 0.40381285, 0.71428016,
       1.67998055, 0.56486474, 0.29504376, 0.99582401, 2.73126254,
       0.63858302, 0.49901512, 0.1625456 , 0.29009947, 0.41233532,
       0.28210301, 0.45252974, 0.98694788, 0.84519324, 0.7537702 ,
       0.4332758 , 0.31252018, 1.44771384, 1.69060132, 0.36163728,
       2.42023675, 0.71811649, 1.46964434, 0.27110825, 0.58160667,
       0.30436792, 0.29073531, 0.89860137, 0.36691359, 0.93946811,
       1.31842005, 0.68243415, 0.1456141 , 1.0761098 , 0.41250304,
       0.44101933, 0.88379282, 0.29976315, 0.48026722, 0.44065701,
       0.3953454 , 0.20985811, 0.11029291, 1.10077171, 0.85881447,
       0.35109906, 0.59462279, 0.79243895, 0.97558335, 1.46257997,
       0.54302544, 0.7990609 , 0.73341026, 0.54179151, 1.22279402,
       0.60442441, 0.99768341, 0.59304996, 0.21894289, 0.62705249,
       0.39311448, 0.70272104, 0.81815493, 0.51196661, 0.63019314,
       0.62327168, 1.09937516, 0.64112772, 0.6290978 , 0.51004568]
    
    # we compare to truncated decimals, so use atol
    assert_allclose(clust.reachability_, np.array(v), atol=1e-5)
