# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD Style.

import numpy as np

from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal

from ..bayes import BayesianRidge, ARDRegression


def test_toy_bayesian_ridge_object():
    """
    Test BayesianRegression ridge classifier
    """
    X = np.array([[1], [2], [6], [8], [10]])
    Y = np.array([1, 2, 6, 8, 10])
    clf = BayesianRidge(compute_score = True)
    clf.fit(X, Y)
    Test = [[1], [3], [4]]
    assert(np.abs(clf.predict(Test)-[1, 3, 4]).sum() < 1.e-2) # identity


def test_toy_ard_object():
    """
    Test BayesianRegression ARD classifier
    """
    X = np.array([[1], [2], [3]])
    Y = np.array([1, 2, 3])
    clf = ARDRegression(compute_score = True)
    clf.fit(X, Y)
    Test = [[1], [3], [4]]
    assert(np.abs(clf.predict(Test)-[1, 3, 4]).sum() < 1.e-3) # identity
