# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Virgile Fritsch <virgile.fritsch@inria.fr>
#
# License: BSD Style.

from numpy.testing import assert_almost_equal, assert_array_almost_equal

import numpy as np

from sklearn import datasets
from sklearn.covariance import empirical_covariance, EmpiricalCovariance, \
    ShrunkCovariance, shrunk_covariance, LedoitWolf, ledoit_wolf, OAS, oas, \
    MinCovDet

X = datasets.load_iris().data
X_1d = X[:, 0]
n_samples, n_features = X.shape


def test_covariance():
    """Tests Covariance module on a simple dataset.

    """
    # test covariance fit from data
    cov = EmpiricalCovariance()
    cov.fit(X)
    assert_array_almost_equal(empirical_covariance(X), cov.covariance_, 4)
    assert_almost_equal(cov.error_norm(empirical_covariance(X)), 0)
    assert_almost_equal(
        cov.error_norm(empirical_covariance(X), norm='spectral'), 0)
    assert_almost_equal(
        cov.error_norm(empirical_covariance(X), norm='frobenius'), 0)
    assert_almost_equal(
        cov.error_norm(empirical_covariance(X), scaling=False), 0)
    assert_almost_equal(
        cov.error_norm(empirical_covariance(X), squared=False), 0)
    # Mahalanobis distances computation test
    mahal_dist = cov.mahalanobis(X)
    assert(np.amax(mahal_dist) < 250)
    assert(np.amin(mahal_dist) > 50)

    # test with n_features = 1
    X_1d = X[:, 0].reshape((-1, 1))
    cov = EmpiricalCovariance()
    cov.fit(X_1d)
    assert_array_almost_equal(empirical_covariance(X_1d), cov.covariance_, 4)
    assert_almost_equal(cov.error_norm(empirical_covariance(X_1d)), 0)
    assert_almost_equal(
        cov.error_norm(empirical_covariance(X_1d), norm='spectral'), 0)

    # test integer type
    X_integer = np.asarray([[0, 1], [1, 0]])
    result = np.asarray([[0.25, -0.25], [-0.25, 0.25]])
    assert_array_almost_equal(empirical_covariance(X_integer), result)


def test_shrunk_covariance():
    """Tests ShrunkCovariance module on a simple dataset.

    """
    # compare shrunk covariance obtained from data and from MLE estimate
    cov = ShrunkCovariance(shrinkage=0.5)
    cov.fit(X)
    assert_array_almost_equal(
        shrunk_covariance(empirical_covariance(X), shrinkage=0.5),
        cov.covariance_, 4)

    # same test with shrinkage not provided
    cov = ShrunkCovariance()
    cov.fit(X)
    assert_array_almost_equal(
        shrunk_covariance(empirical_covariance(X)), cov.covariance_, 4)

    # same test with shrinkage = 0 (<==> empirical_covariance)
    cov = ShrunkCovariance(shrinkage=0.)
    cov.fit(X)
    assert_array_almost_equal(empirical_covariance(X), cov.covariance_, 4)

    # test with n_features = 1
    X_1d = X[:, 0].reshape((-1, 1))
    cov = ShrunkCovariance(shrinkage=0.3)
    cov.fit(X_1d)
    assert_array_almost_equal(empirical_covariance(X_1d), cov.covariance_, 4)

    # test shrinkage coeff on a simple data set (without saving precision)
    cov = ShrunkCovariance(shrinkage=0.5, store_precision=False)
    cov.fit(X)
    assert(cov.precision_ is None)


def test_ledoit_wolf():
    """Tests LedoitWolf module on a simple dataset.

    """
    # test shrinkage coeff on a simple data set
    lw = LedoitWolf()
    lw.fit(X, assume_centered=True)
    assert_almost_equal(lw.shrinkage_, 0.00192, 4)
    assert_almost_equal(lw.score(X, assume_centered=True), -2.89795, 4)
    # compare shrunk covariance obtained from data and from MLE estimate
    lw_cov_from_mle, lw_shinkrage_from_mle = ledoit_wolf(X,
                                                        assume_centered=True)
    assert_array_almost_equal(lw_cov_from_mle, lw.covariance_, 4)
    assert_almost_equal(lw_shinkrage_from_mle, lw.shrinkage_)
    # compare estimates given by LW and ShrunkCovariance
    scov = ShrunkCovariance(shrinkage=lw.shrinkage_)
    scov.fit(X, assume_centered=True)
    assert_array_almost_equal(scov.covariance_, lw.covariance_, 4)

    # test with n_features = 1
    X_1d = X[:, 0].reshape((-1, 1))
    lw = LedoitWolf()
    lw.fit(X_1d, assume_centered=True)
    lw_cov_from_mle, lw_shinkrage_from_mle = ledoit_wolf(X_1d,
                                                         assume_centered=True)
    assert_array_almost_equal(lw_cov_from_mle, lw.covariance_, 4)
    assert_almost_equal(lw_shinkrage_from_mle, lw.shrinkage_)
    assert_array_almost_equal((X_1d ** 2).sum() / n_samples, lw.covariance_, 4)

    # test shrinkage coeff on a simple data set (without saving precision)
    lw = LedoitWolf(store_precision=False)
    lw.fit(X, assume_centered=True)
    assert_almost_equal(lw.score(X, assume_centered=True), -2.89795, 4)
    assert(lw.precision_ is None)

    # Same tests without assuming centered data
    # test shrinkage coeff on a simple data set
    lw = LedoitWolf()
    lw.fit(X)
    assert_almost_equal(lw.shrinkage_, 0.007582, 4)
    assert_almost_equal(lw.score(X), 2.243483, 4)
    # compare shrunk covariance obtained from data and from MLE estimate
    lw_cov_from_mle, lw_shinkrage_from_mle = ledoit_wolf(X)
    assert_array_almost_equal(lw_cov_from_mle, lw.covariance_, 4)
    assert_almost_equal(lw_shinkrage_from_mle, lw.shrinkage_)
    # compare estimates given by LW and ShrunkCovariance
    scov = ShrunkCovariance(shrinkage=lw.shrinkage_)
    scov.fit(X)
    assert_array_almost_equal(scov.covariance_, lw.covariance_, 4)

    # test with n_features = 1
    X_1d = X[:, 0].reshape((-1, 1))
    lw = LedoitWolf()
    lw.fit(X_1d)
    lw_cov_from_mle, lw_shinkrage_from_mle = ledoit_wolf(X_1d)
    assert_array_almost_equal(lw_cov_from_mle, lw.covariance_, 4)
    assert_almost_equal(lw_shinkrage_from_mle, lw.shrinkage_)
    assert_array_almost_equal(empirical_covariance(X_1d), lw.covariance_, 4)

    # test shrinkage coeff on a simple data set (without saving precision)
    lw = LedoitWolf(store_precision=False)
    lw.fit(X)
    assert_almost_equal(lw.score(X), 2.2434839, 4)
    assert(lw.precision_ is None)


def test_oas():
    """Tests OAS module on a simple dataset.

    """
    # test shrinkage coeff on a simple data set
    oa = OAS()
    oa.fit(X, assume_centered=True)
    assert_almost_equal(oa.shrinkage_, 0.018740, 4)
    assert_almost_equal(oa.score(X, assume_centered=True), -5.03605, 4)
    # compare shrunk covariance obtained from data and from MLE estimate
    oa_cov_from_mle, oa_shinkrage_from_mle = oas(X, assume_centered=True)
    assert_array_almost_equal(oa_cov_from_mle, oa.covariance_, 4)
    assert_almost_equal(oa_shinkrage_from_mle, oa.shrinkage_)
    # compare estimates given by OAS and ShrunkCovariance
    scov = ShrunkCovariance(shrinkage=oa.shrinkage_)
    scov.fit(X, assume_centered=True)
    assert_array_almost_equal(scov.covariance_, oa.covariance_, 4)

    # test with n_features = 1
    X_1d = X[:, 0].reshape((-1, 1))
    oa = OAS()
    oa.fit(X_1d, assume_centered=True)
    oa_cov_from_mle, oa_shinkrage_from_mle = oas(X_1d, assume_centered=True)
    assert_array_almost_equal(oa_cov_from_mle, oa.covariance_, 4)
    assert_almost_equal(oa_shinkrage_from_mle, oa.shrinkage_)
    assert_array_almost_equal((X_1d ** 2).sum() / n_samples, oa.covariance_, 4)

    # test shrinkage coeff on a simple data set (without saving precision)
    oa = OAS(store_precision=False)
    oa.fit(X, assume_centered=True)
    assert_almost_equal(oa.score(X, assume_centered=True), -5.03605, 4)
    assert(oa.precision_ is None)

    ### Same tests without assuming centered data
    # test shrinkage coeff on a simple data set
    oa = OAS()
    oa.fit(X)
    assert_almost_equal(oa.shrinkage_, 0.020236, 4)
    assert_almost_equal(oa.score(X), 2.079025, 4)
    # compare shrunk covariance obtained from data and from MLE estimate
    oa_cov_from_mle, oa_shinkrage_from_mle = oas(X)
    assert_array_almost_equal(oa_cov_from_mle, oa.covariance_, 4)
    assert_almost_equal(oa_shinkrage_from_mle, oa.shrinkage_)
    # compare estimates given by OAS and ShrunkCovariance
    scov = ShrunkCovariance(shrinkage=oa.shrinkage_)
    scov.fit(X)
    assert_array_almost_equal(scov.covariance_, oa.covariance_, 4)

    # test with n_features = 1
    X_1d = X[:, 0].reshape((-1, 1))
    oa = OAS()
    oa.fit(X_1d)
    oa_cov_from_mle, oa_shinkrage_from_mle = oas(X_1d)
    assert_array_almost_equal(oa_cov_from_mle, oa.covariance_, 4)
    assert_almost_equal(oa_shinkrage_from_mle, oa.shrinkage_)
    assert_array_almost_equal(empirical_covariance(X_1d), oa.covariance_, 4)

    # test shrinkage coeff on a simple data set (without saving precision)
    oa = OAS(store_precision=False)
    oa.fit(X)
    assert_almost_equal(oa.score(X), 2.079025, 4)
    assert(oa.precision_ is None)


def test_mcd():
    """Tests the FastMCD algorithm implementation

    """
    ### Small data set
    # test without outliers (random independant normal data)
    launch_mcd_on_dataset(100, 5, 0, 0.01, 0.1, 80)
    # test with a contaminated data set (medium contamination)
    launch_mcd_on_dataset(100, 5, 20, 0.01, 0.01, 70)
    # test with a contaminated data set (strong contamination)
    launch_mcd_on_dataset(100, 5, 40, 0.1, 0.1, 50)

    ### Medium data set
    launch_mcd_on_dataset(1000, 5, 450, 1e-3, 1e-3, 540)

    ### Large data set
    launch_mcd_on_dataset(1700, 5, 800, 1e-3, 1e-3, 870)

    ### 1D data set
    launch_mcd_on_dataset(500, 1, 100, 0.001, 0.001, 350)


def launch_mcd_on_dataset(
    n_samples, n_features, n_outliers, tol_loc, tol_cov, tol_support):
    """

    """
    rand_gen = np.random.RandomState(0)
    data = rand_gen.randn(n_samples, n_features)
    # add some outliers
    outliers_index = rand_gen.permutation(n_samples)[:n_outliers]
    outliers_offset = 10. * \
        (rand_gen.randint(2, size=(n_outliers, n_features)) - 0.5)
    data[outliers_index] += outliers_offset
    inliers_mask = np.ones(n_samples).astype(bool)
    inliers_mask[outliers_index] = False

    pure_data = data[inliers_mask]
    # compute MCD by fitting an object
    mcd_fit = MinCovDet().fit(data)
    T = mcd_fit.location_
    S = mcd_fit.covariance_
    H = mcd_fit.support_
    # compare with the estimates learnt from the inliers
    error_location = np.mean((pure_data.mean(0) - T) ** 2)
    assert(error_location < tol_loc)
    error_cov = np.mean((empirical_covariance(pure_data) - S) ** 2)
    assert(error_cov < tol_cov)
    assert(np.sum(H) >= tol_support)
