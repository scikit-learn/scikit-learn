import numpy as np
from numpy.random import randn
from nose.tools import assert_true
from nose.tools import assert_equal

from scipy.sparse import csr_matrix
from numpy.testing import assert_almost_equal

from .. import datasets
from ..pca import PCA
from ..pca import ProbabilisticPCA
from ..pca import RandomizedPCA
from ..pca import _assess_dimension_
from ..pca import _infer_dimension_

iris = datasets.load_iris()


def test_pca():
    """PCA on dense arrays"""
    pca = PCA(n_components=2)
    X = iris.data
    X_r = pca.fit(X).transform(X)
    np.testing.assert_equal(X_r.shape[1], 2)

    pca = PCA()
    pca.fit(X)
    assert_almost_equal(pca.explained_variance_ratio_.sum(), 1.0, 3)


def test_whitening():
    """Check that PCA output has unit-variance"""
    np.random.seed(0)
    n_samples = 100
    n_features = 80
    n_components = 30
    rank = 50

    # some low rank data with correlated features
    X = np.dot(randn(n_samples, rank),
               np.dot(np.diag(np.linspace(10.0, 1.0, rank)),
                      randn(rank, n_features)))
    # the component-wise variance of the first 50 features is 3 times the
    # mean component-wise variance of the remaingin 30 features
    X[:, :50] *= 3

    assert_equal(X.shape, (n_samples, n_features))

    # the component-wise variance is thus highly varying:
    assert_almost_equal(X.std(axis=0).std(), 43.9, 1)

    # whiten the data while projecting to the lower dim subspace
    pca = PCA(n_components=n_components, whiten=True).fit(X)
    X_whitened = pca.transform(X)
    assert_equal(X_whitened.shape, (n_samples, n_components))

    # all output component have unit variances
    assert_almost_equal(X_whitened.std(axis=0), np.ones(n_components))

    # is possible to project on the low dim space without scaling by the
    # singular values
    pca = PCA(n_components=n_components, whiten=False).fit(X)
    X_unwhitened = pca.transform(X)
    assert_equal(X_unwhitened.shape, (n_samples, n_components))

    # in that case the output components still have varying variances
    assert_almost_equal(X_unwhitened.std(axis=0).std(), 74.1, 1)


def test_pca_check_projection():
    """Test that the projection of data is correct"""
    np.random.seed(0)
    n, p = 100, 3
    X = randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    Xt = 0.1 * randn(1, p) + np.array([3, 4, 5])

    Yt = PCA(n_components=2).fit(X).transform(Xt)
    Yt /= np.sqrt((Yt ** 2).sum())

    assert_almost_equal(np.abs(Yt[0][0]), 1., 1)


def test_pca_inverse():
    """Test that the projection of data can be inverted"""
    np.random.seed(0)
    n, p = 50, 3
    X = randn(n, p) # spherical data
    X[:, 1] *= .00001 # make middle component relatively small
    X += [5, 4, 3] # make a large mean

    # same check that we can find the original data from the transformed
    # signal (since the data is almost of rank n_components)
    pca = PCA(n_components=2).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    assert_almost_equal(X, Y_inverse, decimal=3)

    # same as above with whitening (approximate reconstruction)
    pca = PCA(n_components=2, whiten=True)
    pca.fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    relative_max_delta = (np.abs(X - Y_inverse) / np.abs(X).mean()).max()
    assert_almost_equal(relative_max_delta, 0.11, decimal=2)


def test_randomized_pca_check_projection():
    """Test that the projection by RandomizedPCA on dense data is correct"""
    n, p = 100, 3
    X = randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    Xt = 0.1 * randn(1, p) + np.array([3, 4, 5])

    Yt = RandomizedPCA(n_components=2).fit(X).transform(Xt)
    Yt /= np.sqrt((Yt ** 2).sum())

    assert_almost_equal(np.abs(Yt[0][0]), 1., 1)


def test_randomized_pca_check_list():
    """Test that the projection by RandomizedPCA on list data is correct"""
    X = [[1.0, 0.0], [0.0, 1.0]]
    X_transformed = RandomizedPCA(n_components=1).fit(X).transform(X)
    assert_equal(X_transformed.shape, (2, 1))
    assert_almost_equal(X_transformed.mean(), 0.00, 2)
    assert_almost_equal(X_transformed.std(), 0.71, 2)


def test_randomized_pca_inverse():
    """Test that RandomizedPCA is inversible on dense data"""
    np.random.seed(0)
    n, p = 50, 3
    X = randn(n, p) # spherical data
    X[:, 1] *= .00001 # make middle component relatively small
    X += [5, 4, 3] # make a large mean

    # same check that we can find the original data from the transformed signal
    # (since the data is almost of rank n_components)
    pca = RandomizedPCA(n_components=2).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    assert_almost_equal(X, Y_inverse, decimal=2)

    # same as above with whitening (approximate reconstruction)
    pca = RandomizedPCA(n_components=2, whiten=True).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    relative_max_delta = (np.abs(X - Y_inverse) / np.abs(X).mean()).max()
    assert_almost_equal(relative_max_delta, 0.11, decimal=2)


def test_sparse_randomized_pca_check_projection():
    """Test that the projection by RandomizedPCA on sparse data is correct"""
    np.random.seed(0)
    n, p = 100, 3
    X = randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    X = csr_matrix(X)
    Xt = 0.1 * randn(1, p) + np.array([3, 4, 5])
    Xt = csr_matrix(Xt)

    Yt = RandomizedPCA(n_components=2).fit(X).transform(Xt)
    Yt /= np.sqrt((Yt ** 2).sum())

    np.testing.assert_almost_equal(np.abs(Yt[0][0]), 1., 1)


def test_sparse_randomized_pca_inverse():
    """Test that RandomizedPCA is inversible on sparse data"""
    np.random.seed(0)
    n, p = 50, 3
    X = randn(n, p) # spherical data
    X[:, 1] *= .00001 # make middle component relatively small
    # no large means because the sparse version of randomized pca does not do
    # centering to avoid breaking the sparsity
    X = csr_matrix(X)

    # same check that we can find the original data from the transformed signal
    # (since the data is almost of rank n_components)
    pca = RandomizedPCA(n_components=2).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    assert_almost_equal(X.todense(), Y_inverse, decimal=2)

    # same as above with whitening (approximate reconstruction)
    pca = RandomizedPCA(n_components=2, whiten=True).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    relative_max_delta = (np.abs(X.todense() - Y_inverse)
                          / np.abs(X).mean()).max()
    # XXX: this does not seam to work as expected:
    assert_almost_equal(relative_max_delta, 0.91, decimal=2)


def test_pca_dim():
    """Check automated dimensionality setting"""
    np.random.seed(0)
    n, p = 100, 5
    X = randn(n, p) * .1
    X[:10] += np.array([3, 4, 5, 1, 2])
    pca = PCA(n_components='mle').fit(X)
    assert_equal(pca.n_components, 1)


def test_infer_dim_1():
    """TODO: explain what this is testing

    Or at least use explicit variable names...
    """
    n, p = 1000, 5
    X = randn(n, p) * .1 + randn(n, 1) * np.array([3, 4, 5, 1, 2]) \
            + np.array([1, 0, 7, 4, 6])
    pca = PCA(n_components=p)
    pca.fit(X)
    spect = pca.explained_variance_
    ll = []
    for k in range(p):
         ll.append(_assess_dimension_(spect, k, n, p))
    ll = np.array(ll)
    assert_true(ll[1] > ll.max() - .01 * n)


def test_infer_dim_2():
    """TODO: explain what this is testing

    Or at least use explicit variable names...
    """
    n, p = 1000, 5
    X = randn(n, p) * .1
    X[:10] += np.array([3, 4, 5, 1, 2])
    X[10:20] += np.array([6, 0, 7, 2, -1])
    pca = PCA(n_components=p)
    pca.fit(X)
    spect = pca.explained_variance_
    assert_true(_infer_dimension_(spect, n, p) > 1)


def test_infer_dim_3():
    """
    """
    n, p = 100, 5
    X = randn(n, p)*.1
    X[:10] += np.array([3, 4, 5, 1, 2])
    X[10:20] += np.array([6, 0, 7, 2, -1])
    X[30:40] += 2*np.array([-1, 1, -1, 1, -1])
    pca = PCA(n_components=p)
    pca.fit(X)
    spect = pca.explained_variance_
    assert_true(_infer_dimension_(spect, n, p) > 2)

def test_infer_dim_by_explained_variance():
    X = iris.data
    pca = PCA(n_components=0.95)
    pca.fit(X)
    assert_equal(pca.n_components, 2)

    pca = PCA(n_components=0.01)
    pca.fit(X)
    assert_equal(pca.n_components, 1)


def test_probabilistic_pca_1():
    """Test that probabilistic PCA yields a reasonable score"""
    n, p = 1000, 3
    X = randn(n, p)*.1 + np.array([3, 4, 5])
    ppca = ProbabilisticPCA(n_components=2)
    ppca.fit(X)
    ll1 = ppca.score(X)
    h = 0.5 * np.log(2 * np.pi * np.exp(1) / 0.1 ** 2) * p
    np.testing.assert_almost_equal(ll1.mean() / h, 1, 0)


def test_probabilistic_pca_2():
    """Test that probabilistic PCA correctly separated different datasets"""
    n, p = 100, 3
    X = randn(n, p) * .1 + np.array([3, 4, 5])
    ppca = ProbabilisticPCA(n_components=2)
    ppca.fit(X)
    ll1 = ppca.score(X)
    ll2 = ppca.score(randn(n, p) * .2 + np.array([3, 4, 5]))
    assert_true(ll1.mean() > ll2.mean())


def test_probabilistic_pca_3():
    """The homoscedastic model should work slightly worth
    than the heteroscedastic one in over-fitting condition
    """
    n, p = 100, 3
    X = randn(n, p)*.1 + np.array([3, 4, 5])
    ppca = ProbabilisticPCA(n_components=2)
    ppca.fit(X)
    ll1 = ppca.score(X)
    ppca.fit(X, homoscedastic=False)
    ll2 = ppca.score(X)
    assert_true(ll1.mean() < ll2.mean())


def test_probabilistic_pca_4():
    """Check that ppca select the right model"""
    n, p = 200, 3
    Xl = randn(n, p) + randn(n, 1) * np.array([3, 4, 5]) + np.array([1, 0, 7])
    Xt = randn(n, p) + randn(n, 1) * np.array([3, 4, 5]) + np.array([1, 0, 7])
    ll = np.zeros(p)
    for k in range(p):
        ppca = ProbabilisticPCA(n_components=k)
        ppca.fit(Xl)
        ll[k] = ppca.score(Xt).mean()

    assert_true(ll.argmax() == 1)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
