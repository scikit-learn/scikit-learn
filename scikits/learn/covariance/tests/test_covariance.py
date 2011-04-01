# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Virgile Fritsch <virgile.fritsch@inria.fr>
#
# License: BSD Style.

from numpy.testing import assert_almost_equal, assert_array_almost_equal

from .. import Covariance
from .. import ShrunkCovariance, shrunk_covariance
from .. import LedoitWolf, ledoit_wolf

import numpy as np
from scikits.learn import datasets

X = datasets.load_iris().data
n_samples = X.shape[0]

def test_covariance():
    """Tests Covariance module on a simple dataset.
    
    """
    # test covariance fit from data
    cov = Covariance()
    cov.fit(X)
    assert_array_almost_equal(np.dot(X.T, X) / n_samples, cov.covariance_, 4)

def test_shrunk_covariance():
    """Tests ShrunkCovariance module on a simple dataset.
    
    """
    # compare shrunk covariance obtained from data and from MLE estimate
    cov = ShrunkCovariance(shrinkage=0.5)
    cov.fit(X)
    assert_array_almost_equal(
        shrunk_covariance(np.dot(X.T, X) / n_samples, shrinkage=0.5),
        cov.covariance_, 4
        )
    
    # same test with shrinkage not provided
    cov = ShrunkCovariance()
    cov.fit(X)
    assert_array_almost_equal(
        shrunk_covariance(np.dot(X.T, X) / n_samples),
        cov.covariance_, 4
        )

    # test with n_features = 1
    X_1d = X[:,0]
    cov = ShrunkCovariance()
    cov.fit(X_1d)
    assert_array_almost_equal(
        shrunk_covariance(np.dot(X_1d.T, X_1d) / n_samples),
        cov.covariance_, 4
        )

def test_ledoitWolf():
    """Tests LedoitWolf module on a simple dataset.

    """
    # test shrinkage coeff on a simple data set
    lw = LedoitWolf()
    lw.fit(X)
    assert_almost_equal(lw.shrinkage_, 0.00192, 4)
    assert_almost_equal(lw.score(X), -2.89795, 4)

    # compare shrunk covariance obtained from data and from MLE estimate
    lw_cov_from_mle, lw_shinkrage_from_mle = ledoit_wolf(X)
    assert_array_almost_equal(lw_cov_from_mle, lw.covariance_, 4)
    assert_almost_equal(lw_shinkrage_from_mle, lw.shrinkage_)
    
    # compare estimates given by LW and ShrunkCovariance
    scov = ShrunkCovariance(shrinkage=lw.shrinkage_)
    scov.fit(X)
    assert_array_almost_equal(scov.covariance_, lw.covariance_, 4)

    # test with n_features = 1
    X_1d = X[:,0]
    lw = LedoitWolf()
    lw.fit(X_1d)
    lw_cov_from_mle, lw_shinkrage_from_mle = ledoit_wolf(X_1d)
    assert_array_almost_equal(lw_cov_from_mle, lw.covariance_, 4)
    assert_almost_equal(lw_shinkrage_from_mle, lw.shrinkage_)
