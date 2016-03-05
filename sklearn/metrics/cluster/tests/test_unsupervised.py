import numpy as np
from scipy.sparse import csr_matrix

from sklearn import datasets
from sklearn.metrics.cluster.unsupervised import silhouette_score
from sklearn.metrics import pairwise_distances
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regexp


def test_silhouette():
    # Tests the Silhouette Coefficient.
    dataset = datasets.load_iris()
    X = dataset.data
    y = dataset.target
    D = pairwise_distances(X, metric='euclidean')
    # Given that the actual labels are used, we can assume that S would be
    # positive.
    silhouette = silhouette_score(D, y, metric='precomputed')
    assert(silhouette > 0)
    # Test without calculating D
    silhouette_metric = silhouette_score(X, y, metric='euclidean')
    assert_almost_equal(silhouette, silhouette_metric)
    # Test with sampling
    silhouette = silhouette_score(D, y, metric='precomputed',
                                  sample_size=int(X.shape[0] / 2),
                                  random_state=0)
    silhouette_metric = silhouette_score(X, y, metric='euclidean',
                                         sample_size=int(X.shape[0] / 2),
                                         random_state=0)
    assert(silhouette > 0)
    assert(silhouette_metric > 0)
    assert_almost_equal(silhouette_metric, silhouette)
    # Test with sparse X
    X_sparse = csr_matrix(X)
    D = pairwise_distances(X_sparse, metric='euclidean')
    silhouette = silhouette_score(D, y, metric='precomputed')
    assert(silhouette > 0)


def test_no_nan():
    # Assert Silhouette Coefficient != nan when there is 1 sample in a class.
    # This tests for the condition that caused issue 960.
    # Note that there is only one sample in cluster 0. This used to cause the
    # silhouette_score to return nan (see bug #960).
    labels = np.array([1, 0, 1, 1, 1])
    # The distance matrix doesn't actually matter.
    D = np.random.RandomState(0).rand(len(labels), len(labels))
    silhouette = silhouette_score(D, labels, metric='precomputed')
    assert_false(np.isnan(silhouette))


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
        silhouette_score(X, labels + 10), silhouette_score(X, labels))


def test_non_numpy_labels():
    dataset = datasets.load_iris()
    X = dataset.data
    y = dataset.target
    assert_equal(
        silhouette_score(list(X), list(y)), silhouette_score(X, y))
