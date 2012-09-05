import numpy as np
from scipy.sparse import csr_matrix

from .... import datasets
from ..unsupervised import silhouette_score
from ... import pairwise_distances
from nose.tools import assert_false


def test_silhouette():
    """Tests the Silhouette Coefficient. """
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
    assert(silhouette == silhouette_metric)
    # Test with sampling
    silhouette = silhouette_score(D, y, metric='precomputed',
                                  sample_size=int(X.shape[0] / 2))
    assert(silhouette > 0)
    # Test with sparse X
    X_sparse = csr_matrix(X)
    D = pairwise_distances(X_sparse, metric='euclidean')
    silhouette = silhouette_score(D, y, metric='precomputed')
    assert(silhouette > 0)


def test_no_nan():
    """Assert Silhouette Coefficient != nan when there is 1 sample in a class.

        This tests for the condition that caused issue 960.
    """
    # Note that there is only one sample in cluster 0. This used to cause the
    # silhouette_score to return nan (see bug #960).
    labels = np.array([1, 0, 1, 1, 1])
    # The distance matrix doesn't actually matter.
    D = np.random.RandomState(0).rand(len(labels), len(labels))
    silhouette = silhouette_score(D, labels, metric='precomputed')
    assert_false(np.isnan(silhouette))
