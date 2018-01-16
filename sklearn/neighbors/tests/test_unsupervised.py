import numpy as np

from sklearn.metrics import euclidean_distances
from sklearn.neighbors.unsupervised import NearestNeighborsTransformer

from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_warns_message


def test_transformer_query_parameters():
    # Test the check that one of radius and n_neighbors is given, but not both
    X = np.random.randn(10, 3)

    msg = 'Please specify either radius or n_neighbors.'
    assert_warns_message(UserWarning, msg,
                         NearestNeighborsTransformer().fit_transform, X)
    msg = 'Please do not specify both radius and n_neighbors.'
    assert_raise_message(ValueError, msg,
                         NearestNeighborsTransformer(
                             n_neighbors=5, radius=1).fit_transform, X)


def test_transformer_nonzero():
    # Test the number of neighbors returned
    n_neighbors = 5
    X = np.random.randn(20, 10)
    radius = np.percentile(euclidean_distances(X), 10)

    # with n_neighbors specified
    for mode in ('distance', 'connectivity'):
        Xt = NearestNeighborsTransformer(n_neighbors=n_neighbors,
                                         mode=mode).fit_transform(X)
        assert Xt.data.shape == (X.shape[0] * n_neighbors, )

    # with radius specified
    for mode in ('distance', 'connectivity'):
        Xt = NearestNeighborsTransformer(radius=radius,
                                         mode=mode).fit_transform(X)
        assert not Xt.data.shape == (X.shape[0] * n_neighbors, )


def test_transformer_shape():
    # Test the shape of the returned array
    n_samples_fit, n_samples_transform, n_features = 100, 80, 10
    X_fit = np.random.randn(n_samples_fit, n_features)
    X_transform = np.random.randn(n_samples_transform, n_features)
    est = NearestNeighborsTransformer(n_neighbors=5)

    Xt = est.fit_transform(X_fit)
    assert Xt.shape == (n_samples_fit, n_samples_fit)

    Xt = est.transform(X_transform)
    assert Xt.shape == (n_samples_transform, n_samples_fit)
