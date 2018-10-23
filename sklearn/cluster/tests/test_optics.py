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
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster.dbscan_ import DBSCAN
from sklearn.utils.testing import assert_equal, assert_warns
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_allclose

from sklearn.cluster.tests.common import generate_clustered_data


rng = np.random.RandomState(0)
n_points_per_cluster = 50
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
    clust = OPTICS(max_eps=5.0 * 6.0, min_samples=4)
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
    msg = "Specify an epsilon smaller than 0.15. Got 0.3."
    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers,
                                cluster_std=0.4, random_state=0)

    # Compute OPTICS
    clust = OPTICS(max_eps=5.0 * 0.03, min_samples=10)
    clust2 = clust.fit(X)
    assert_raise_message(ValueError, msg, clust2.extract_dbscan, 0.3)


def test_bad_reachability():
    msg = "All reachability values are inf. Set a larger max_eps."
    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers,
                                cluster_std=0.4, random_state=0)

    clust = OPTICS(max_eps=5.0 * 0.003, min_samples=10)
    assert_raise_message(ValueError, msg, clust.fit, X)


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


def test_compare_to_ELKI():
    # Expected values, computed with (future) ELKI 0.7.5 using:
    # java -jar elki.jar cli -dbc.in csv -dbc.filter FixedDBIDsFilter
    #   -algorithm clustering.optics.OPTICSHeap -optics.minpts 5
    # where the FixedDBIDsFilter gives 0-indexed ids.
    r = [np.inf, 0.7865694338710508, 0.4373157299595305, 0.4121908069391695,
         0.302907091394212, 0.20815674060999778, 0.20815674060999778,
         0.15190193459676368, 0.15190193459676368, 0.28229645104833345,
         0.302907091394212, 0.30507239477026865, 0.30820580778767087,
         0.3289019667317037, 0.3458462228589966, 0.3458462228589966,
         0.2931114364132193, 0.2931114364132193, 0.2562790168458507,
         0.23654635530592025, 0.37903448688824876, 0.3920764620583683,
         0.4121908069391695, 0.4364542226186831, 0.45523658462146793,
         0.458757846268185, 0.458757846268185, 0.4752907412198826,
         0.42350366820623375, 0.42350366820623375, 0.42350366820623375,
         0.47758738570352993, 0.47758738570352993, 0.4776963110272057,
         0.5272079288923731, 0.5591861752070968, 0.5592057084987357,
         0.5609913790596295, 0.5909117211348757, 0.5940470220777727,
         0.5940470220777727, 0.6861627576116127, 0.687795873252133,
         0.7538541412862811, 0.7865694338710508, 0.8038180561910464,
         0.8038180561910464, 0.8242451615289921, 0.8548361202185057,
         0.8790098789921685, 2.9281214555815764, 1.3256656984284734,
         0.19590944671099267, 0.1339924636672767, 0.1137384200258616,
         0.061455005237474075, 0.061455005237474075, 0.061455005237474075,
         0.045627777293497276, 0.045627777293497276, 0.045627777293497276,
         0.04900902556283447, 0.061455005237474075, 0.06225461602815799,
         0.06835750467748272, 0.07882900172724974, 0.07882900172724974,
         0.07650735397943846, 0.07650735397943846, 0.07650735397943846,
         0.07650735397943846, 0.07650735397943846, 0.07113275489288699,
         0.07890196345324527, 0.07052683707634783, 0.07052683707634783,
         0.07052683707634783, 0.08284027053523288, 0.08725436842020087,
         0.08725436842020087, 0.09010229261951723, 0.09128578974358925,
         0.09154172670176584, 0.0968576383038391, 0.12007572768323092,
         0.12024155806196564, 0.12141990481584404, 0.1339924636672767,
         0.13694322786307633, 0.14275793459246572, 0.15093125027309579,
         0.17927454395170142, 0.18151803569400365, 0.1906028449191095,
         0.1906028449191095, 0.19604486784973194, 0.2096539172540186,
         0.2096539172540186, 0.21614333983312325, 0.22036454909290296,
         0.23610322103910933, 0.26028003932256766, 0.2607126030060721,
         0.2891824876072483, 0.3258089271514364, 0.35968687619960743,
         0.4512973330510512, 0.4746141313843085, 0.5958585488429471,
         0.6468718886525733, 0.6878453052524358, 0.6911582799500199,
         0.7172169499815705, 0.7209874999572031, 0.6326884657912096,
         0.5755681293026617, 0.5755681293026617, 0.5755681293026617,
         0.6015042225447333, 0.6756244556376542, 0.4722384908959966,
         0.08775739179493615, 0.06665303472021758, 0.056308477780164796,
         0.056308477780164796, 0.05507767260835565, 0.05368146914586802,
         0.05163427719303039, 0.05163427719303039, 0.05163427719303039,
         0.04918757627098621, 0.04918757627098621, 0.05368146914586802,
         0.05473720349424546, 0.05473720349424546, 0.048442038421760626,
         0.048442038421760626, 0.04598840269934622, 0.03984301937835033,
         0.04598840269934622, 0.04598840269934622, 0.04303884892957088,
         0.04303884892957088, 0.04303884892957088, 0.0431802780806032,
         0.0520412490141781, 0.056308477780164796, 0.05080724020124642,
         0.05080724020124642, 0.05080724020124642, 0.06385565101399236,
         0.05840878369200427, 0.0474472391259039, 0.0474472391259039,
         0.04232512684465669, 0.04232512684465669, 0.04232512684465669,
         0.0474472391259039, 0.051802632822946656, 0.051802632822946656,
         0.05316405104684577, 0.05316405104684577, 0.05840878369200427,
         0.06385565101399236, 0.08025248922898705, 0.08775739179493615,
         0.08993337040710143, 0.08993337040710143, 0.08993337040710143,
         0.08993337040710143, 0.297457175321605, 0.29763608186278934,
         0.3415255849656254, 0.34713336941105105, 0.44108940848708167,
         0.35942962652965604, 0.35942962652965604, 0.33609522256535296,
         0.5008111387107295, 0.5333587622018111, 0.6223243743872802,
         0.6793840035409552, 0.7445032492109848, 0.7445032492109848,
         0.6556432627279256, 0.6556432627279256, 0.6556432627279256,
         0.8196566935960162, 0.8724089149982351, 0.9352758042365477,
         0.9352758042365477, 1.0581847953137133, 1.0684332509194163,
         1.0887817699873303, 1.2552604310322708, 1.3993856001769436,
         1.4869615658197606, 1.6588098267326852, 1.679969559453028,
         1.679969559453028, 1.6860509219163458, 1.6860509219163458,
         1.1465697826627317, 0.992866533434785, 0.7691908270707519,
         0.578131499171622, 0.578131499171622, 0.578131499171622,
         0.5754243919945694, 0.8416199360035114, 0.8722493727270406,
         0.9156549976203665, 0.9156549976203665, 0.7472322844356064,
         0.715219324518981, 0.715219324518981, 0.715219324518981,
         0.7472322844356064, 0.820988298336316, 0.908958489674247,
         0.9234036745782839, 0.9519521817942455, 0.992866533434785,
         0.992866533434785, 0.9995692674695029, 1.0727415198904493,
         1.1395519941203158, 1.1395519941203158, 1.1741737271442092,
         1.212860115632712, 0.8724097897372123, 0.8724097897372123,
         0.8724097897372123, 1.2439272570611581, 1.2439272570611581,
         1.3524538390109015, 1.3524538390109015, 1.2982303284415664,
         1.3610655849680207, 1.3802783392089437, 1.3802783392089437,
         1.4540636953090629, 1.5879329500533819, 1.5909193228826986,
         1.72931779186001, 1.9619075944592093, 2.1994355761906257,
         2.2508672067362165, 2.274436122235927, 2.417635732260135,
         3.014235905390584, 0.30616929141177107, 0.16449675872754976,
         0.09071681523805683, 0.09071681523805683, 0.09071681523805683,
         0.08727060912039632, 0.09151721189581336, 0.12277953408786725,
         0.14285575406641507, 0.16449675872754976, 0.16321992344119793,
         0.1330971730344373, 0.11429891993167259, 0.11429891993167259,
         0.11429891993167259, 0.11429891993167259, 0.11429891993167259,
         0.0945498340011516, 0.11410457435712089, 0.1196414019798306,
         0.12925682285016715, 0.12925682285016715, 0.12925682285016715,
         0.12864887158869853, 0.12864887158869853, 0.12864887158869853,
         0.13369634918690246, 0.14330826543275352, 0.14877705862323184,
         0.15203263952428328, 0.15696350160889708, 0.1585326700393211,
         0.1585326700393211, 0.16034306786654595, 0.16034306786654595,
         0.15053328296567992, 0.16396729418886688, 0.16763548009617293,
         0.1732029325454474, 0.21163390061029352, 0.21497664171864372,
         0.22125889949299, 0.240251070192081, 0.240251070192081,
         0.2413620965310808, 0.26319419022234064, 0.26319419022234064,
         0.27989712380504483, 0.2909782800714374]
    o = [0, 3, 6, 7, 15, 4, 27, 28, 49, 17, 35, 47, 46, 39, 13, 19,
         22, 29, 30, 38, 34, 32, 43, 8, 25, 9, 37, 23, 33, 40, 44, 11, 36, 5,
         45, 48, 41, 26, 24, 20, 31, 2, 16, 10, 18, 14, 42, 12, 1, 21, 234,
         132, 112, 115, 107, 110, 120, 114, 100, 131, 137, 145, 130, 121, 134,
         116, 149, 108, 111, 113, 142, 148, 119, 104, 126, 133, 138, 127, 101,
         105, 103, 106, 125, 140, 123, 147, 144, 129, 141, 117, 143, 136, 128,
         122, 124, 102, 109, 249, 146, 118, 135, 245, 139, 224, 241, 217, 202,
         248, 233, 214, 236, 211, 206, 231, 212, 221, 229, 244, 208, 226, 83,
         76, 53, 77, 88, 62, 66, 65, 89, 93, 79, 95, 74, 70, 82, 51, 73, 87,
         67, 94, 56, 52, 63, 80, 75, 57, 96, 60, 69, 90, 86, 58, 68, 81, 64,
         84, 85, 97, 59, 98, 61, 71, 78, 92, 50, 91, 55, 54, 72, 99, 210, 201,
         216, 239, 203, 218, 219, 222, 240, 294, 243, 246, 204, 220, 200, 215,
         230, 225, 205, 207, 237, 223, 235, 209, 228, 238, 227, 285, 232, 256,
         281, 270, 260, 252, 272, 268, 292, 298, 269, 275, 257, 250, 284, 283,
         286, 295, 297, 293, 289, 258, 299, 282, 262, 296, 287, 267, 255, 263,
         288, 276, 251, 266, 274, 271, 277, 261, 279, 290, 253, 254, 291, 259,
         280, 278, 273, 247, 265, 242, 264, 213, 199, 174, 154, 152, 180, 186,
         195, 170, 181, 176, 187, 173, 157, 159, 158, 172, 182, 183, 151, 197,
         177, 160, 156, 171, 175, 184, 193, 161, 179, 196, 185, 192, 165, 166,
         164, 189, 155, 162, 188, 153, 178, 169, 194, 150, 163, 198, 190, 191,
         168, 167]
    p = [-1, 0, 3, 6, 7, 15, 15, 27, 27, 4, 7, 49, 47, 4, 39, 39,
         19, 19, 29, 30, 30, 13, 6, 43, 34, 32, 32, 25, 23, 23, 23, 3, 11, 46,
         46, 45, 9, 38, 33, 26, 26, 8, 20, 33, 0, 18, 18, 2, 18, 44, 0, 234,
         132, 112, 115, 107, 107, 120, 114, 114, 114, 114, 107, 100, 100, 134,
         134, 149, 149, 108, 108, 108, 148, 113, 104, 104, 104, 142, 127, 127,
         126, 138, 126, 148, 127, 148, 127, 112, 147, 116, 117, 101, 145, 128,
         128, 122, 136, 136, 249, 102, 102, 118, 143, 146, 245, 123, 139, 241,
         241, 217, 248, 202, 248, 224, 231, 212, 212, 212, 229, 229, 226, 83,
         76, 53, 53, 88, 62, 66, 66, 66, 93, 93, 79, 93, 70, 82, 82, 73, 87,
         73, 94, 56, 56, 56, 63, 67, 53, 96, 96, 96, 69, 86, 58, 58, 81, 81,
         81, 58, 64, 64, 59, 59, 86, 69, 78, 83, 84, 55, 55, 55, 72, 50, 201,
         210, 216, 203, 203, 219, 54, 240, 239, 240, 236, 236, 220, 220, 220,
         217, 139, 243, 243, 204, 211, 246, 215, 223, 294, 209, 227, 227, 209,
         281, 270, 260, 252, 272, 272, 272, 298, 269, 298, 275, 275, 284, 283,
         283, 283, 284, 286, 283, 298, 299, 260, 260, 250, 299, 258, 258, 296,
         250, 276, 276, 276, 289, 289, 267, 267, 279, 261, 277, 277, 258, 266,
         290, 209, 207, 290, 228, 278, 228, 290, 199, 174, 154, 154, 154, 186,
         186, 180, 170, 174, 187, 173, 157, 159, 157, 157, 159, 183, 183, 172,
         197, 160, 160, 171, 171, 171, 151, 173, 184, 151, 196, 185, 185, 179,
         179, 189, 177, 165, 175, 162, 164, 181, 169, 169, 181, 178, 178, 178,
         168]

    # Tests against known extraction array
    # Does NOT work with metric='euclidean', because sklearn euclidean has
    # worse numeric precision. 'minkowski' is slower but more accurate.
    clust = OPTICS(min_samples=5).fit(X)

    assert_array_equal(clust.ordering_, np.array(o))
    assert_array_equal(clust.predecessor_[clust.ordering_], np.array(p))
    assert_allclose(clust.reachability_[clust.ordering_], np.array(r))
    # ELKI currently does not print the core distances (which are not used much
    # in literature, but we can at least ensure to have this consistency:
    for i in clust.ordering_[1:]:
        assert (clust.reachability_[i] >=
                clust.core_distances_[clust.predecessor_[i]])


def test_precomputed_dists():
    redX = X[::10]
    dists = pairwise_distances(redX, metric='euclidean')
    clust1 = OPTICS(min_samples=10, algorithm='brute',
                    metric='precomputed').fit(dists)
    clust2 = OPTICS(min_samples=10, algorithm='brute',
                    metric='euclidean').fit(redX)

    assert_allclose(clust1.reachability_, clust2.reachability_)
    assert_array_equal(clust1.labels_, clust2.labels_)


def test_processing_order():
    """Early dev version of OPTICS would not consider all unprocessed points,
    but only direct neighbors. This tests against this mistake."""
    Y = [[0], [10], [-10], [25]]
    clust = OPTICS(min_samples=3, max_eps=15).fit(Y)
    assert_array_equal(clust.reachability_, [np.inf, 10, 10, 15])
    assert_array_equal(clust.core_distances_, [10, 15, 20, 25])
    assert_array_equal(clust.ordering_, [0, 1, 2, 3])
