import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix

from sklearn import datasets
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_greater
from sklearn.metrics.cluster import silhouette_score
from sklearn.metrics.cluster import silhouette_samples
from sklearn.metrics import pairwise_distances
from sklearn.metrics.cluster import calinski_harabaz_score


def test_silhouette():
    # Tests the Silhouette Coefficient.
    dataset = datasets.load_iris()
    X_dense = dataset.data
    X_csr = csr_matrix(X_dense)
    X_dok = sp.dok_matrix(X_dense)
    X_lil = sp.lil_matrix(X_dense)
    y = dataset.target

    for X in [X_dense, X_csr, X_dok, X_lil]:
        D = pairwise_distances(X, metric='euclidean')
        # Given that the actual labels are used, we can assume that S would be
        # positive.
        score_precomputed = silhouette_score(D, y, metric='precomputed')
        assert_greater(score_precomputed, 0)
        # Test without calculating D
        score_euclidean = silhouette_score(X, y, metric='euclidean')
        assert_almost_equal(score_precomputed, score_euclidean)

        if X is X_dense:
            score_dense_without_sampling = score_precomputed
        else:
            assert_almost_equal(score_euclidean,
                                score_dense_without_sampling)

        # Test with sampling
        score_precomputed = silhouette_score(D, y, metric='precomputed',
                                             sample_size=int(X.shape[0] / 2),
                                             random_state=0)
        score_euclidean = silhouette_score(X, y, metric='euclidean',
                                           sample_size=int(X.shape[0] / 2),
                                           random_state=0)
        assert_greater(score_precomputed, 0)
        assert_greater(score_euclidean, 0)
        assert_almost_equal(score_euclidean, score_precomputed)

        if X is X_dense:
            score_dense_with_sampling = score_precomputed
        else:
            assert_almost_equal(score_euclidean, score_dense_with_sampling)


def test_cluster_size_1():
    # Assert Silhouette Coefficient == 0 when there is 1 sample in a cluster
    # (cluster 0). We also test the case where there are identical samples
    # as the only members of a cluster (cluster 2). To our knowledge, this case
    # is not discussed in reference material, and we choose for it a sample
    # score of 1.
    X = [[0.], [1.], [1.], [2.], [3.], [3.]]
    labels = np.array([0, 1, 1, 1, 2, 2])

    # Cluster 0: 1 sample -> score of 0 by Rousseeuw's convention
    # Cluster 1: intra-cluster = [.5, .5, 1]
    #            inter-cluster = [1, 1, 1]
    #            silhouette    = [.5, .5, 0]
    # Cluster 2: intra-cluster = [0, 0]
    #            inter-cluster = [arbitrary, arbitrary]
    #            silhouette    = [1., 1.]

    silhouette = silhouette_score(X, labels)
    assert_false(np.isnan(silhouette))
    ss = silhouette_samples(X, labels)
    assert_array_equal(ss, [0, .5, .5, 0, 1, 1])


def test_correct_labelsize():
    # Assert 1 < n_labels < n_samples
    dataset = datasets.load_iris()
    X = dataset.data

    # n_labels = n_samples
    y = np.arange(X.shape[0])
    assert_raises_regexp(ValueError,
                         'Number of labels is %d\. Valid values are 2 '
                         'to n_samples - 1 \(inclusive\)' % len(np.unique(y)),
                         silhouette_score, X, y)

    # n_labels = 1
    y = np.zeros(X.shape[0])
    assert_raises_regexp(ValueError,
                         'Number of labels is %d\. Valid values are 2 '
                         'to n_samples - 1 \(inclusive\)' % len(np.unique(y)),
                         silhouette_score, X, y)


def test_non_encoded_labels():
    dataset = datasets.load_iris()
    X = dataset.data
    labels = dataset.target
    assert_equal(
        silhouette_score(X, labels * 2 + 10), silhouette_score(X, labels))
    assert_array_equal(
        silhouette_samples(X, labels * 2 + 10), silhouette_samples(X, labels))


def test_non_numpy_labels():
    dataset = datasets.load_iris()
    X = dataset.data
    y = dataset.target
    assert_equal(
        silhouette_score(list(X), list(y)), silhouette_score(X, y))


def test_calinski_harabaz_score():
    rng = np.random.RandomState(seed=0)

    # Assert message when there is only one label
    assert_raise_message(ValueError, "Number of labels is",
                         calinski_harabaz_score,
                         rng.rand(10, 2), np.zeros(10))

    # Assert message when all point are in different clusters
    assert_raise_message(ValueError, "Number of labels is",
                         calinski_harabaz_score,
                         rng.rand(10, 2), np.arange(10))

    # Assert the value is 1. when all samples are equals
    assert_equal(1., calinski_harabaz_score(np.ones((10, 2)),
                                            [0] * 5 + [1] * 5))

    # Assert the value is 0. when all the mean cluster are equal
    assert_equal(0., calinski_harabaz_score([[-1, -1], [1, 1]] * 10,
                                            [0] * 10 + [1] * 10))

    # General case (with non numpy arrays)
    X = ([[0, 0], [1, 1]] * 5 + [[3, 3], [4, 4]] * 5 +
         [[0, 4], [1, 3]] * 5 + [[3, 1], [4, 0]] * 5)
    labels = [0] * 10 + [1] * 10 + [2] * 10 + [3] * 10
    assert_almost_equal(calinski_harabaz_score(X, labels),
                        45 * (40 - 4) / (5 * (4 - 1)))
