from scipy.sparse import csr_matrix

from ... import datasets
from ..unsupervised import silhouette_coefficient
from .. import pairwise_distances


def test_silhouette():
    """ Tests the Silhouette Coefficient. """
    dataset = datasets.load_iris()
    X = dataset.data
    y = dataset.target
    D = pairwise_distances(X, metric='euclidean')
    # Given that the actual labels are used, we can assume that S would be
    # positive.
    silhouette = silhouette_coefficient(y, D)
    assert(silhouette > 0)
    # Test with sparse X
    X_sparse = csr_matrix(X)
    D = pairwise_distances(X, metric='euclidean')
    silhouette = silhouette_coefficient(y, D)
    assert(silhouette > 0)
