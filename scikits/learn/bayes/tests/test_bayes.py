import numpy as np
from scikits.learn.bayes.bayes import bayesian_ridge, BayesianRegression
from numpy.testing import assert_array_almost_equal
from  scikits.learn.datasets.samples_generator import linear

def test_toy():
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    w ,log_likelihood = bayesian_ridge(X, Y)
    assert_array_almost_equal(w, [1])

def test_toy_object():
    """
    Test BayesianRegression classifier
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = BayesianRegression()
    clf.fit(X, Y)
    Test = [[1], [2], [3], [4]]
    assert_array_almost_equal(clf.predict(Test), [1, 2, 3, 4]) # identity

def test_simu_object():
    """
    Test BayesianRegression classifier with simulated linear data
    """
    X,Y  = linear.sparse_uncorrelated(nb_samples=100,nb_features=10)
    clf = BayesianRegression()
    clf.fit(X, Y)
    Xtest,Ytest  = linear.sparse_uncorrelated(nb_samples=100,nb_features=10)
    mse = np.mean((clf.predict(Xtest)-Ytest)**2)
    assert(mse<2.)

