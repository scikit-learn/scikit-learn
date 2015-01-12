import numpy as np
from scipy.sparse import csr_matrix

from .... import datasets
from ..unsupervised import silhouette_score, silhouette_samples
from ... import pairwise_distances
from sklearn.utils.testing import assert_false, assert_almost_equal
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
    # The distance matrix needs only be a valid distance matrix
    D = pairwise_distances(np.random.RandomState(0).rand(len(labels), 2))
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


def test_paper_example():
    lower = [5.58,
             7.00, 6.50,
             7.08, 7.00, 3.83,
             4.83, 5.08, 8.17, 5.83,
             2.17, 5.75, 6.67, 6.92, 4.92,
             6.42, 5.00, 5.58, 6.00, 4.67, 6.42,
             3.42, 5.50, 6.42, 6.42, 5.00, 3.92, 6.17,
             2.50, 4.92, 6.25, 7.33, 4.50, 2.25, 6.33, 2.75,
             6.08, 6.67, 4.25, 2.67, 6.00, 6.17, 6.17, 6.92, 6.17,
             5.25, 6.83, 4.50, 3.75, 5.75, 5.42, 6.08, 5.83, 6.67, 3.67,
             4.75, 3.00, 6.08, 6.67, 5.00, 5.58, 4.83, 6.17, 5.67, 6.50, 6.92]
    D = np.zeros((12, 12))
    D[np.tril_indices(12, -1)] = lower
    D += D.T

    names = ['BEL', 'BRA', 'CHI', 'CUB', 'EGY', 'FRA', 'IND', 'ISR', 'USA',
             'USS', 'YUG', 'ZAI']

    labels1 = [1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1]
    labels2 = [1, 2, 3, 3, 1, 1, 2, 1, 1, 3, 3, 2]

    expected1 = {'USA': .43, 'BEL': .39, 'FRA': .35, 'ISR': .30, 'BRA': .22,
                 'EGY': .20, 'ZAI': .19, 'CUB': .40, 'USS': .34, 'CHI': .33,
                 'YUG': .26, 'IND': -.04}
    score1 = .28
    expected2 = {'USA': .47, 'FRA': .44, 'BEL': .42, 'ISR': .37, 'EGY': .02,
                 'ZAI': .28, 'BRA': .25, 'IND': .17, 'CUB': .48, 'USS': .44,
                 'YUG': .31, 'CHI': .31}
    score2 = .33

    for labels, expected, score in [(labels1, expected1, score1),
                                    (labels2, expected2, score2)]:
        expected = [expected[name] for name in names]
        assert_almost_equal(expected, silhouette_samples(D, np.array(labels),
                                                         metric='precomputed'),
                            decimal=2)
        assert_almost_equal(score, silhouette_score(D, np.array(labels),
                                                    metric='precomputed'),
                            decimal=2)


def test_sample_weight():
    n_samples = 10
    rng = np.random.RandomState(0)
    X = rng.rand(n_samples, 2)
    sample_weight = rng.randint(3, size=n_samples)
    # Ensure both the 0 and multiple cases are handled
    assert np.any(sample_weight == 2)
    assert np.any(sample_weight == 0)
    y = rng.randint(2, size=n_samples)
    X_repeated = np.repeat(X, sample_weight, axis=0)
    y_repeated = np.repeat(y, sample_weight, axis=0)
    repeated_input = silhouette_samples(X_repeated, y_repeated)
    weighed_s = silhouette_samples(X, y, sample_weight=sample_weight)
    repeated_output = np.repeat(weighed_s, sample_weight, axis=0)
    assert_almost_equal(repeated_input, repeated_output)
    assert_almost_equal(silhouette_score(X_repeated, y_repeated),
                        silhouette_score(X, y, sample_weight=sample_weight))
