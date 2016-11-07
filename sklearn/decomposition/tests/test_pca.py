import numpy as np
import scipy as sp
from itertools import product

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_no_warnings
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_less

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.decomposition import RandomizedPCA
from sklearn.decomposition.pca import _assess_dimension_
from sklearn.decomposition.pca import _infer_dimension_

iris = datasets.load_iris()
solver_list = ['full', 'arpack', 'randomized', 'auto']


def test_pca():
    # PCA on dense arrays
    X = iris.data

    for n_comp in np.arange(X.shape[1]):
        pca = PCA(n_components=n_comp, svd_solver='full')

        X_r = pca.fit(X).transform(X)
        np.testing.assert_equal(X_r.shape[1], n_comp)

        X_r2 = pca.fit_transform(X)
        assert_array_almost_equal(X_r, X_r2)

        X_r = pca.transform(X)
        X_r2 = pca.fit_transform(X)
        assert_array_almost_equal(X_r, X_r2)

        # Test get_covariance and get_precision
        cov = pca.get_covariance()
        precision = pca.get_precision()
        assert_array_almost_equal(np.dot(cov, precision),
                                  np.eye(X.shape[1]), 12)

    # test explained_variance_ratio_ == 1 with all components
    pca = PCA(svd_solver='full')
    pca.fit(X)
    assert_almost_equal(pca.explained_variance_ratio_.sum(), 1.0, 3)


def test_pca_arpack_solver():
    # PCA on dense arrays
    X = iris.data
    d = X.shape[1]

    # Loop excluding the extremes, invalid inputs for arpack
    for n_comp in np.arange(1, d):
        pca = PCA(n_components=n_comp, svd_solver='arpack', random_state=0)

        X_r = pca.fit(X).transform(X)
        np.testing.assert_equal(X_r.shape[1], n_comp)

        X_r2 = pca.fit_transform(X)
        assert_array_almost_equal(X_r, X_r2)

        X_r = pca.transform(X)
        assert_array_almost_equal(X_r, X_r2)

        # Test get_covariance and get_precision
        cov = pca.get_covariance()
        precision = pca.get_precision()
        assert_array_almost_equal(np.dot(cov, precision),
                                  np.eye(d), 12)

    pca = PCA(n_components=0, svd_solver='arpack', random_state=0)
    assert_raises(ValueError, pca.fit, X)
    # Check internal state
    assert_equal(pca.n_components,
                 PCA(n_components=0,
                     svd_solver='arpack', random_state=0).n_components)
    assert_equal(pca.svd_solver,
                 PCA(n_components=0,
                     svd_solver='arpack', random_state=0).svd_solver)

    pca = PCA(n_components=d, svd_solver='arpack', random_state=0)
    assert_raises(ValueError, pca.fit, X)
    assert_equal(pca.n_components,
                 PCA(n_components=d,
                     svd_solver='arpack', random_state=0).n_components)
    assert_equal(pca.svd_solver,
                 PCA(n_components=0,
                     svd_solver='arpack', random_state=0).svd_solver)


def test_pca_randomized_solver():
    # PCA on dense arrays
    X = iris.data

    # Loop excluding the 0, invalid for randomized
    for n_comp in np.arange(1, X.shape[1]):
        pca = PCA(n_components=n_comp, svd_solver='randomized', random_state=0)

        X_r = pca.fit(X).transform(X)
        np.testing.assert_equal(X_r.shape[1], n_comp)

        X_r2 = pca.fit_transform(X)
        assert_array_almost_equal(X_r, X_r2)

        X_r = pca.transform(X)
        assert_array_almost_equal(X_r, X_r2)

        # Test get_covariance and get_precision
        cov = pca.get_covariance()
        precision = pca.get_precision()
        assert_array_almost_equal(np.dot(cov, precision),
                                  np.eye(X.shape[1]), 12)

    pca = PCA(n_components=0, svd_solver='randomized', random_state=0)
    assert_raises(ValueError, pca.fit, X)

    pca = PCA(n_components=0, svd_solver='randomized', random_state=0)
    assert_raises(ValueError, pca.fit, X)
    # Check internal state
    assert_equal(pca.n_components,
                 PCA(n_components=0,
                     svd_solver='randomized', random_state=0).n_components)
    assert_equal(pca.svd_solver,
                 PCA(n_components=0,
                     svd_solver='randomized', random_state=0).svd_solver)


def test_no_empty_slice_warning():
    # test if we avoid numpy warnings for computing over empty arrays
    n_components = 10
    n_features = n_components + 2  # anything > n_comps triggered it in 0.16
    X = np.random.uniform(-1, 1, size=(n_components, n_features))
    pca = PCA(n_components=n_components)
    assert_no_warnings(pca.fit, X)


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
    # mean component-wise variance of the remaining 30 features
    X[:, :50] *= 3

    assert_equal(X.shape, (n_samples, n_features))

    # the component-wise variance is thus highly varying:
    assert_greater(X.std(axis=0).std(), 43.8)

    for solver, copy in product(solver_list, (True, False)):
        # whiten the data while projecting to the lower dim subspace
        X_ = X.copy()  # make sure we keep an original across iterations.
        pca = PCA(n_components=n_components, whiten=True, copy=copy,
                  svd_solver=solver, random_state=0, iterated_power=7)
        # test fit_transform
        X_whitened = pca.fit_transform(X_.copy())
        assert_equal(X_whitened.shape, (n_samples, n_components))
        X_whitened2 = pca.transform(X_)
        assert_array_almost_equal(X_whitened, X_whitened2)

        assert_almost_equal(X_whitened.std(axis=0), np.ones(n_components),
                            decimal=6)
        assert_almost_equal(X_whitened.mean(axis=0), np.zeros(n_components))

        X_ = X.copy()
        pca = PCA(n_components=n_components, whiten=False, copy=copy,
                  svd_solver=solver).fit(X_)
        X_unwhitened = pca.transform(X_)
        assert_equal(X_unwhitened.shape, (n_samples, n_components))

        # in that case the output components still have varying variances
        assert_almost_equal(X_unwhitened.std(axis=0).std(), 74.1, 1)
        # we always center, so no test for non-centering.


# Ignore warnings from switching to more power iterations in randomized_svd
@ignore_warnings
def test_explained_variance():
    # Check that PCA output has unit-variance
    rng = np.random.RandomState(0)
    n_samples = 100
    n_features = 80

    X = rng.randn(n_samples, n_features)

    pca = PCA(n_components=2, svd_solver='full').fit(X)
    apca = PCA(n_components=2, svd_solver='arpack', random_state=0).fit(X)
    assert_array_almost_equal(pca.explained_variance_,
                              apca.explained_variance_, 1)
    assert_array_almost_equal(pca.explained_variance_ratio_,
                              apca.explained_variance_ratio_, 3)

    rpca = PCA(n_components=2, svd_solver='randomized', random_state=42).fit(X)
    assert_array_almost_equal(pca.explained_variance_,
                              rpca.explained_variance_, 1)
    assert_array_almost_equal(pca.explained_variance_ratio_,
                              rpca.explained_variance_ratio_, 1)

    # compare to empirical variances
    X_pca = pca.transform(X)
    assert_array_almost_equal(pca.explained_variance_,
                              np.var(X_pca, axis=0))

    X_pca = apca.transform(X)
    assert_array_almost_equal(apca.explained_variance_,
                              np.var(X_pca, axis=0))

    X_rpca = rpca.transform(X)
    assert_array_almost_equal(rpca.explained_variance_, np.var(X_rpca, axis=0),
                              decimal=1)

    # Same with correlated data
    X = datasets.make_classification(n_samples, n_features,
                                     n_informative=n_features-2,
                                     random_state=rng)[0]

    pca = PCA(n_components=2).fit(X)
    rpca = PCA(n_components=2, svd_solver='randomized',
               random_state=rng).fit(X)
    assert_array_almost_equal(pca.explained_variance_ratio_,
                              rpca.explained_variance_ratio_, 5)


def test_pca_check_projection():
    # Test that the projection of data is correct
    rng = np.random.RandomState(0)
    n, p = 100, 3
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    Xt = 0.1 * rng.randn(1, p) + np.array([3, 4, 5])

    for solver in solver_list:
        Yt = PCA(n_components=2, svd_solver=solver).fit(X).transform(Xt)
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
    pca = PCA(n_components=2, svd_solver='full').fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    assert_almost_equal(X, Y_inverse, decimal=3)

    # same as above with whitening (approximate reconstruction)
    for solver in solver_list:
        pca = PCA(n_components=2, whiten=True, svd_solver=solver)
        pca.fit(X)
        Y = pca.transform(X)
        Y_inverse = pca.inverse_transform(Y)
        assert_almost_equal(X, Y_inverse, decimal=3)


def test_pca_validation():
    X = [[0, 1], [1, 0]]
    for solver in solver_list:
        for n_components in [-1, 3]:
            assert_raises(ValueError,
                          PCA(n_components, svd_solver=solver).fit, X)


def test_randomized_pca_check_projection():
    # Test that the projection by randomized PCA on dense data is correct
    rng = np.random.RandomState(0)
    n, p = 100, 3
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5])
    Xt = 0.1 * rng.randn(1, p) + np.array([3, 4, 5])

    Yt = PCA(n_components=2, svd_solver='randomized',
             random_state=0).fit(X).transform(Xt)
    Yt /= np.sqrt((Yt ** 2).sum())

    assert_almost_equal(np.abs(Yt[0][0]), 1., 1)


def test_randomized_pca_check_list():
    # Test that the projection by randomized PCA on list data is correct
    X = [[1.0, 0.0], [0.0, 1.0]]
    X_transformed = PCA(n_components=1, svd_solver='randomized',
                        random_state=0).fit(X).transform(X)
    assert_equal(X_transformed.shape, (2, 1))
    assert_almost_equal(X_transformed.mean(), 0.00, 2)
    assert_almost_equal(X_transformed.std(), 0.71, 2)


def test_randomized_pca_inverse():
    # Test that randomized PCA is inversible on dense data
    rng = np.random.RandomState(0)
    n, p = 50, 3
    X = rng.randn(n, p)  # spherical data
    X[:, 1] *= .00001  # make middle component relatively small
    X += [5, 4, 3]  # make a large mean

    # same check that we can find the original data from the transformed signal
    # (since the data is almost of rank n_components)
    pca = PCA(n_components=2, svd_solver='randomized', random_state=0).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    assert_almost_equal(X, Y_inverse, decimal=2)

    # same as above with whitening (approximate reconstruction)
    pca = PCA(n_components=2, whiten=True, svd_solver='randomized',
              random_state=0).fit(X)
    Y = pca.transform(X)
    Y_inverse = pca.inverse_transform(Y)
    relative_max_delta = (np.abs(X - Y_inverse) / np.abs(X).mean()).max()
    assert_less(relative_max_delta, 1e-5)


def test_pca_dim():
    # Check automated dimensionality setting
    rng = np.random.RandomState(0)
    n, p = 100, 5
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5, 1, 2])
    pca = PCA(n_components='mle', svd_solver='full').fit(X)
    assert_equal(pca.n_components, 'mle')
    assert_equal(pca.n_components_, 1)


def test_infer_dim_1():
    # TODO: explain what this is testing
    # Or at least use explicit variable names...
    n, p = 1000, 5
    rng = np.random.RandomState(0)
    X = (rng.randn(n, p) * .1 + rng.randn(n, 1) * np.array([3, 4, 5, 1, 2]) +
         np.array([1, 0, 7, 4, 6]))
    pca = PCA(n_components=p, svd_solver='full')
    pca.fit(X)
    spect = pca.explained_variance_
    ll = []
    for k in range(p):
        ll.append(_assess_dimension_(spect, k, n, p))
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
    pca = PCA(n_components=p, svd_solver='full')
    pca.fit(X)
    spect = pca.explained_variance_
    assert_greater(_infer_dimension_(spect, n, p), 1)


def test_infer_dim_3():
    n, p = 100, 5
    rng = np.random.RandomState(0)
    X = rng.randn(n, p) * .1
    X[:10] += np.array([3, 4, 5, 1, 2])
    X[10:20] += np.array([6, 0, 7, 2, -1])
    X[30:40] += 2 * np.array([-1, 1, -1, 1, -1])
    pca = PCA(n_components=p, svd_solver='full')
    pca.fit(X)
    spect = pca.explained_variance_
    assert_greater(_infer_dimension_(spect, n, p), 2)


def test_infer_dim_by_explained_variance():
    X = iris.data
    pca = PCA(n_components=0.95, svd_solver='full')
    pca.fit(X)
    assert_equal(pca.n_components, 0.95)
    assert_equal(pca.n_components_, 2)

    pca = PCA(n_components=0.01, svd_solver='full')
    pca.fit(X)
    assert_equal(pca.n_components, 0.01)
    assert_equal(pca.n_components_, 1)

    rng = np.random.RandomState(0)
    # more features than samples
    X = rng.rand(5, 20)
    pca = PCA(n_components=.5, svd_solver='full').fit(X)
    assert_equal(pca.n_components, 0.5)
    assert_equal(pca.n_components_, 2)


def test_pca_score():
    # Test that probabilistic PCA scoring yields a reasonable score
    n, p = 1000, 3
    rng = np.random.RandomState(0)
    X = rng.randn(n, p) * .1 + np.array([3, 4, 5])
    for solver in solver_list:
        pca = PCA(n_components=2, svd_solver=solver)
        pca.fit(X)
        ll1 = pca.score(X)
        h = -0.5 * np.log(2 * np.pi * np.exp(1) * 0.1 ** 2) * p
        np.testing.assert_almost_equal(ll1 / h, 1, 0)


def test_pca_score2():
    # Test that probabilistic PCA correctly separated different datasets
    n, p = 100, 3
    rng = np.random.RandomState(0)
    X = rng.randn(n, p) * .1 + np.array([3, 4, 5])
    for solver in solver_list:
        pca = PCA(n_components=2, svd_solver=solver)
        pca.fit(X)
        ll1 = pca.score(X)
        ll2 = pca.score(rng.randn(n, p) * .2 + np.array([3, 4, 5]))
        assert_greater(ll1, ll2)

        # Test that it gives different scores if whiten=True
        pca = PCA(n_components=2, whiten=True, svd_solver=solver)
        pca.fit(X)
        ll2 = pca.score(X)
        assert_true(ll1 > ll2)


def test_pca_score3():
    # Check that probabilistic PCA selects the right model
    n, p = 200, 3
    rng = np.random.RandomState(0)
    Xl = (rng.randn(n, p) + rng.randn(n, 1) * np.array([3, 4, 5]) +
          np.array([1, 0, 7]))
    Xt = (rng.randn(n, p) + rng.randn(n, 1) * np.array([3, 4, 5]) +
          np.array([1, 0, 7]))
    ll = np.zeros(p)
    for k in range(p):
        pca = PCA(n_components=k, svd_solver='full')
        pca.fit(Xl)
        ll[k] = pca.score(Xt)

    assert_true(ll.argmax() == 1)


def test_svd_solver_auto():
    rng = np.random.RandomState(0)
    X = rng.uniform(size=(1000, 50))

    # case: n_components in (0,1) => 'full'
    pca = PCA(n_components=.5)
    pca.fit(X)
    pca_test = PCA(n_components=.5, svd_solver='full')
    pca_test.fit(X)
    assert_array_almost_equal(pca.components_, pca_test.components_)

    # case: max(X.shape) <= 500 => 'full'
    pca = PCA(n_components=5, random_state=0)
    Y = X[:10, :]
    pca.fit(Y)
    pca_test = PCA(n_components=5, svd_solver='full', random_state=0)
    pca_test.fit(Y)
    assert_array_almost_equal(pca.components_, pca_test.components_)

    # case: n_components >= .8 * min(X.shape) => 'full'
    pca = PCA(n_components=50)
    pca.fit(X)
    pca_test = PCA(n_components=50, svd_solver='full')
    pca_test.fit(X)
    assert_array_almost_equal(pca.components_, pca_test.components_)

    # n_components >= 1 and n_components < .8 * min(X.shape) => 'randomized'
    pca = PCA(n_components=10, random_state=0)
    pca.fit(X)
    pca_test = PCA(n_components=10, svd_solver='randomized', random_state=0)
    pca_test.fit(X)
    assert_array_almost_equal(pca.components_, pca_test.components_)


def test_deprecation_randomized_pca():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))

    depr_message = ("Class RandomizedPCA is deprecated; RandomizedPCA was "
                    "deprecated in 0.18 and will be "
                    "removed in 0.20. Use PCA(svd_solver='randomized') "
                    "instead. The new implementation DOES NOT store "
                    "whiten ``components_``. Apply transform to get them.")

    def fit_deprecated(X):
        global Y
        rpca = RandomizedPCA(random_state=0)
        Y = rpca.fit_transform(X)

    assert_warns_message(DeprecationWarning, depr_message, fit_deprecated, X)
    Y_pca = PCA(svd_solver='randomized', random_state=0).fit_transform(X)
    assert_array_almost_equal(Y, Y_pca)


def test_pca_spase_input():

    X = np.random.RandomState(0).rand(5, 4)
    X = sp.sparse.csr_matrix(X)
    assert(sp.sparse.issparse(X))

    for svd_solver in solver_list:
        pca = PCA(n_components=3, svd_solver=svd_solver)

        assert_raises(TypeError, pca.fit, X)
