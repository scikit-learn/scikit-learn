from scipy.sparse import csr_matrix

from .... import datasets
from ..unsupervised import silhouette_score
from ... import pairwise_distances


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
