# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Virgile Fritsch <virgile.fritsch@inria.fr>
#
# License: BSD Style.

from numpy.testing import assert_almost_equal, assert_array_almost_equal, \
    assert_raises

import numpy as np

from sklearn import datasets
from sklearn.covariance import empirical_covariance, MinCovDet, \
    EllipticEnvelope

X = datasets.load_iris().data
X_1d = X[:, 0]
n_samples, n_features = X.shape


def test_mcd():
    """Tests the FastMCD algorithm implementation

    """
    ### Small data set
    # test without outliers (random independent normal data)
    launch_mcd_on_dataset(100, 5, 0, 0.01, 0.1, 80)
    # test with a contaminated data set (medium contamination)
    launch_mcd_on_dataset(100, 5, 20, 0.01, 0.01, 70)
    # test with a contaminated data set (strong contamination)
    launch_mcd_on_dataset(100, 5, 40, 0.1, 0.1, 50)

    ### Medium data set
    launch_mcd_on_dataset(1000, 5, 450, 0.1, 0.1, 540)

    ### Large data set
    launch_mcd_on_dataset(1700, 5, 800, 0.1, 0.1, 870)

    ### 1D data set
    launch_mcd_on_dataset(500, 1, 100, 0.001, 0.001, 350)


def launch_mcd_on_dataset(
    n_samples, n_features, n_outliers, tol_loc, tol_cov, tol_support):

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
    mcd_fit = MinCovDet(random_state=rand_gen).fit(data)
    T = mcd_fit.location_
    S = mcd_fit.covariance_
    H = mcd_fit.support_
    # compare with the estimates learnt from the inliers
    error_location = np.mean((pure_data.mean(0) - T) ** 2)
    assert(error_location < tol_loc)
    error_cov = np.mean((empirical_covariance(pure_data) - S) ** 2)
    assert(error_cov < tol_cov)
    assert(np.sum(H) >= tol_support)
    assert_array_almost_equal(mcd_fit.mahalanobis(data), mcd_fit.dist_)


def test_mcd_issue1127():
    # Check that the code does not break with X.shape = (3, 1)
    # (i.e. n_support = n_samples)
    rnd = np.random.RandomState(0)
    X = rnd.normal(size=(3, 1))
    mcd = MinCovDet()
    mcd.fit(X)


def test_outlier_detection():
    rnd = np.random.RandomState(0)
    X = rnd.randn(100, 10)
    clf = EllipticEnvelope(contamination=0.1)
    print clf.threshold
    assert_raises(Exception, clf.predict, X)
    assert_raises(Exception, clf.decision_function, X)
    clf.fit(X)
    y_pred = clf.predict(X)
    decision = clf.decision_function(X, raw_values=True)
    decision_transformed = clf.decision_function(X, raw_values=False)

    assert_array_almost_equal(
        decision, clf.mahalanobis(X))
    assert_array_almost_equal(clf.mahalanobis(X), clf.dist_)
    assert_almost_equal(clf.score(X, np.ones(100)),
                        (100 - y_pred[y_pred == -1].size) / 100.)
    assert(sum(y_pred == -1) == sum(decision_transformed < 0))
