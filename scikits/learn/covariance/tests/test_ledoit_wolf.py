# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#
# License: BSD Style.

from numpy.testing import assert_almost_equal, assert_array_almost_equal

from .. import LedoitWolf, Covariance, shrunk_covariance

import numpy as np
from scikits.learn import datasets

X = datasets.load_iris().data
n_samples = X.shape[0]

def test_Covariance():
    """
    Test Covariance on a simple dataset.
    """
    cov = Covariance()
    cov.fit(X)
    assert_array_almost_equal(np.dot(X.T, X) / n_samples, cov.covariance_, 4)

def test_LedoitWolf():
    """
    Test LedoitWolf on a simple dataset.
    """
    lw = LedoitWolf()
    lw.fit(X)
    assert_almost_equal(lw.shrinkage_, 0.00192, 4)
    assert_almost_equal(lw.score(X), -2.89795, 4)
    assert_array_almost_equal(shrunk_covariance(X, shrinkage=lw.shrinkage_),
                              lw.covariance_, 4)
    emp_cov = np.dot(X.T, X) / n_samples
    assert_array_almost_equal(shrunk_covariance(emp_cov,
                    shrinkage=lw.shrinkage_, data_is_cov=True),
                              lw.covariance_, 3)
