import numpy as np

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.decomposition import RandomizedPCA
from sklearn.decomposition.pca import _assess_dimension
from sklearn.decomposition.pca import _infer_dimension

iris = datasets.load_iris()


def test_pca():
    # PCA on dense arrays
    pca = PCA(n_components=2)
    X = iris.data
    X_r = pca.fit(X).transform(X)
    np.testing.assert_equal(X_r.shape[1], 2)

    X_r2 = pca.fit_transform(X)
    assert_array_almost_equal(X_r, X_r2)

    pca = PCA()
    pca.fit(X)
    assert_almost_equal(pca.explained_variance_ratio_.sum(), 1.0, 3)

    X_r = pca.transform(X)
    X_r2 = pca.fit_transform(X)

    assert_array_almost_equal(X_r, X_r2)

    # Test get_covariance and get_precision with n_components == n_features
    # with n_components < n_features and with n_components == 0
    for n_components in [0, 2, X.shape[1]]:
        pca.n_components = n_components
        pca.fit(X)
        cov = pca.get_covariance()
        precision = pca.get_precision()
        assert_array_almost_equal(np.dot(cov, precision),
                                  np.eye(X.shape[1]), 12)


def test_whitening():
    # Check that PCA output has unit-variance
    rng = np.random.RandomState(0)
    n_samples = 100
    n_features = 80
    n_components = 30
    rank = 50

    # some low rank data with correlated features
    X = np.dot(rng.randn(n_samples, rank),
               np.dot(np.diag(np.linspace(10.0, 1.0, rank)),
                      rng.randn(rank, n_features)))
    # the component-wise variance of the first 50 features is 3 times the
    # mean component-wise variance of the remaingin 30 features
    X[:, :50] *= 3

    assert_equal(X.shape, (n_samples, n_features))

    # the component-wise variance is thus highly varying:
    assert_almost_equal(X.std(axis=0).std(), 43.9, 1)

    for this_PCA, copy in [(x, y) for x in (PCA, RandomizedPCA)
                           for y in (True, False)]:
        # whiten the data while projecting to the lower dim subspace
        X_ = X.copy()  # make sure we keep an original across iterations.
        pca = this_PCA(n_components=n_components, whiten=True, copy=copy)
        # test fit_transform
        X_whitened = pca.fit_transform(X_.copy())
        assert_equal(X_whitened.shape, (n_samples, n_components))
        X_whitened2 = pca.transform(X_)
        assert_array_almost_equal(X_whitened, X_whitened2)

        assert_almost_equal(X_whitened.std(axis=0), np.ones(n_components))
        assert_almost_equal(X_whitened.mean(axis=0), np.zeros(n_components))

        X_ = X.copy()
        pca = this_PCA(n_components=n_components, whiten=False,
                       copy=copy).fit(X_)
        X_unwhitened = pca.transform(X_)
        assert_equal(X_unwhitened.shape, (n_samples, n_components))

        # in that case the output components still have varying variances
        assert_almost_equal(X_unwhitened.std(axis=0).std(), 74.1, 1)
        # we always center, so no test for non-centering.


def test_explained_variance():
    # Check that PCA output has unit-variance
    rng = np.random.RandomState(0)
    n_samples = 100
    n_features = 80

    X = rng.randn(n_samples, n_features)

    pca = PCA(n_components=2).fit(X)
    rpca = RandomizedPCA(n_components=2, random_state=42).fit(X)
    assert_array_almost_equal(pca.explained_variance_,
                              rpca.explained_variance_, 1)
    assert_array_almost_equal(pca.explained_variance_ratio_,
                              rpca.explained_variance_ratio_, 3)

    # compare to empirical variances
    X_pca = pca.transform(X)
    assert_array_almost_equal(pca.explained_variance_,
                              np.var(X_pca, axis=0))

    X_rpca = rpca.transform(X)
    assert_array_almost_equal(rpca.explained_variance_,
                              np.var(X_rpca, axis=0))


def test_pca_check_projection():
    # Test that the projection of data is correct
    rng = np.random.RandomState(0)
    n, p = 100, 3
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    Xt = 0.1 * rng.randn(1, p) + np.array([3, 4, 5])

    Yt = PCA(n_components=2).fit(X).transform(Xt)
    Yt /= np.sqrt((Yt ** 2).sum())

    assert_almost_equal(np.abs(Yt[0][0]), 1., 1)


def test_pca_inverse():
    # Test that the projection of data can be inverted
    rng = np.random.RandomState(0)
    n, p = 50, 3
    X = rng.randn(n, p)  # spherical data
    X[:, 1] *= .00001  # make middle component relatively small
    X += [5, 4, 3]  # make a large mean

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
    assert_almost_equal(X, Y_inverse, decimal=3)


def test_pca_validation():
    X = [[0, 1], [1, 0]]
    for n_components in [-1, 3]:
        assert_raises(ValueError, PCA(n_components).fit, X)


def test_randomized_pca_check_projection():
    # Test that the projection by RandomizedPCA on dense data is correct
    rng = np.random.RandomState(0)
    n, p = 100, 3
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    Xt = 0.1 * rng.randn(1, p) + np.array([3, 4, 5])

    Yt = RandomizedPCA(n_components=2, random_state=0).fit(X).transform(Xt)
    Yt /= np.sqrt((Yt ** 2).sum())

    assert_almost_equal(np.abs(Yt[0][0]), 1., 1)


def test_randomized_pca_check_list():
    # Test that the projection by RandomizedPCA on list data is correct
    X = [[1.0, 0.0], [0.0, 1.0]]
    X_transformed = RandomizedPCA(n_components=1,
                                  random_state=0).fit(X).transform(X)
    assert_equal(X_transformed.shape, (2, 1))
    assert_almost_equal(X_transformed.mean(), 0.00, 2)
    assert_almost_equal(X_transformed.std(), 0.71, 2)


def test_randomized_pca_inverse():
    # Test that RandomizedPCA is invertible on dense data
    rng = np.random.RandomState(0)
    n, p = 50, 3
    X = rng.randn(n, p)  # spherical data
    X[:, 1] *= .00001  # make middle component relatively small
    X += [5, 4, 3]  # make a large mean

    # same check that we can find the original data from the transformed signal
    # (since the data is almost of rank n_components)
    pca = RandomizedPCA(n_components=2, random_state=0).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    assert_almost_equal(X, Y_inverse, decimal=2)

    # same as above with whitening (approximate reconstruction)
    pca = RandomizedPCA(n_components=2, whiten=True,
                        random_state=0).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    relative_max_delta = (np.abs(X - Y_inverse) / np.abs(X).mean()).max()
    assert_almost_equal(relative_max_delta, 0.11, decimal=2)


def test_pca_dim():
    # Check automated dimensionality setting
    rng = np.random.RandomState(0)
    n, p = 100, 5
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5, 1, 2])
    pca = PCA(n_components='mle').fit(X)
    assert_equal(pca.n_components, 'mle')
    assert_equal(pca.n_components_, 1)


def test_infer_dim_1():
    # TODO: explain what this is testing
    # Or at least use explicit variable names...
    n, p = 1000, 5
    rng = np.random.RandomState(0)
    X = (rng.randn(n, p) * .1 + rng.randn(n, 1) * np.array([3, 4, 5, 1, 2])
         + np.array([1, 0, 7, 4, 6]))
    pca = PCA(n_components=p)
    pca.fit(X)
    spect = pca.explained_variance_
    ll = []
    for k in range(p):
        ll.append(_assess_dimension(spect, k, n, p))
    ll = np.array(ll)
    assert_greater(ll[1], ll.max() - .01 * n)


def test_infer_dim_2():
    # TODO: explain what this is testing
    # Or at least use explicit variable names...
    n, p = 1000, 5
    rng = np.random.RandomState(0)
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5, 1, 2])
    X[10:20] += np.array([6, 0, 7, 2, -1])
    pca = PCA(n_components=p)
    pca.fit(X)
    spect = pca.explained_variance_
    assert_greater(_infer_dimension(spect, n, p), 1)


def test_infer_dim_3():
    n, p = 100, 5
    rng = np.random.RandomState(0)
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5, 1, 2])
    X[10:20] += np.array([6, 0, 7, 2, -1])
    X[30:40] += 2 * np.array([-1, 1, -1, 1, -1])
    pca = PCA(n_components=p)
    pca.fit(X)
    spect = pca.explained_variance_
    assert_greater(_infer_dimension(spect, n, p), 2)


def test_infer_dim_bad_spec():
    """
    Test that we deal with a spectrum that drops to near zero
    see issue https://github.com/scikit-learn/scikit-learn/issues/4441
    """
    spectrum = np.array([1, 1e-30, 1e-30, 1e-30])
    n_samples = 10
    n_features = 5
    ret = _infer_dimension(spectrum, n_samples, n_features)
    assert_equal(ret, 0)


def test_assess_dimension_tiny_eigenvals():
    """
    Test that we deal with tiny eigenvalues appropriately when
    `mle` inferring `n_components`
    see issue https://github.com/scikit-learn/scikit-learn/issues/4441
    """
    spectrum = np.array([1, 1e-30, 1e-30, 1e-30])
    n_samples = 10
    n_features = 5
    rank = 4
    ret = _assess_dimension(spectrum, rank, n_samples, n_features)
    assert_equal(ret, -np.inf)


def test_infer_dim_mle():
    """
    Test that we deal with tiny eigenvalues appropriately when
    `mle` inferring `n_components` with a pathalogical `X` dastaset
    see issue https://github.com/scikit-learn/scikit-learn/issues/4441
    """
    X, _ = datasets.make_classification(n_informative=1, n_repeated=18,
                                        n_redundant=1, n_clusters_per_class=1,
                                        random_state=20150609)
    pca = PCA(n_components='mle').fit(X)
    assert_equal(pca.n_components_, 0)


def test_infer_dim_by_explained_variance():
    X = iris.data
    pca = PCA(n_components=0.95)
    pca.fit(X)
    assert_equal(pca.n_components, 0.95)
    assert_equal(pca.n_components_, 2)

    pca = PCA(n_components=0.01)
    pca.fit(X)
    assert_equal(pca.n_components, 0.01)
    assert_equal(pca.n_components_, 1)

    rng = np.random.RandomState(0)
    # more features than samples
    X = rng.rand(5, 20)
    pca = PCA(n_components=.5).fit(X)
    assert_equal(pca.n_components, 0.5)
    assert_equal(pca.n_components_, 2)


def test_pca_score():
    # Test that probabilistic PCA scoring yields a reasonable score
    n, p = 1000, 3
    rng = np.random.RandomState(0)
    X = rng.randn(n, p) * .1 + np.array([3, 4, 5])
    pca = PCA(n_components=2)
    pca.fit(X)
    ll1 = pca.score(X)
    h = -0.5 * np.log(2 * np.pi * np.exp(1) * 0.1 ** 2) * p
    np.testing.assert_almost_equal(ll1 / h, 1, 0)


def test_pca_score2():
    # Test that probabilistic PCA correctly separated different datasets
    n, p = 100, 3
    rng = np.random.RandomState(0)
    X = rng.randn(n, p) * .1 + np.array([3, 4, 5])
    pca = PCA(n_components=2)
    pca.fit(X)
    ll1 = pca.score(X)
    ll2 = pca.score(rng.randn(n, p) * .2 + np.array([3, 4, 5]))
    assert_greater(ll1, ll2)

    # Test that it gives the same scores if whiten=True
    pca = PCA(n_components=2, whiten=True)
    pca.fit(X)
    ll2 = pca.score(X)
    assert_almost_equal(ll1, ll2)


def test_pca_score3():
    # Check that probabilistic PCA selects the right model
    n, p = 200, 3
    rng = np.random.RandomState(0)
    Xl = (rng.randn(n, p) + rng.randn(n, 1) * np.array([3, 4, 5])
          + np.array([1, 0, 7]))
    Xt = (rng.randn(n, p) + rng.randn(n, 1) * np.array([3, 4, 5])
          + np.array([1, 0, 7]))
    ll = np.zeros(p)
    for k in range(p):
        pca = PCA(n_components=k)
        pca.fit(Xl)
        ll[k] = pca.score(Xt)

    assert_true(ll.argmax() == 1)
