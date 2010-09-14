import numpy as np
from numpy.testing import (assert_array_almost_equal,
                           assert_almost_equal)

from nose.tools import assert_equal, assert_true

from ..lars import lars_path, LeastAngleRegression, LassoLARS, LARS
from ..coordinate_descent import Lasso

from scikits.learn import datasets

n, m = 10, 10
np.random.seed (0)
diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target


def test_1():
    """
    Principle of LARS is to keep covariances tied and decreasing
    """
    max_pred = 10
    alphas_, active, coef_path_ = lars_path(diabetes.data, diabetes.target, max_pred, method="lar")
    for (i, coef_) in enumerate(coef_path_.T):
        res =  y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[ C - eps < abs(cov)])
        if i < max_pred:
            assert ocur == i+1
        else:
            # no more than max_pred variables can go into the active set
            assert ocur == max_pred


def test_lars_lstsq():
    """
    Test that LARS gives least square solution at the end
    of the path
    """
    # test that it arrives to a least squares solution
    alphas_, active, coef_path_ = lars_path(diabetes.data, diabetes.target, method="lar")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_path_.T[-1], coef_lstsq)


def test_lasso_gives_lstsq_solution():
    """
    Test that LARS Lasso gives least square solution at the end
    of the path
    """

    alphas_, active, coef_path_ = lars_path(X, y, max_iter=12, method="lasso")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_lstsq , coef_path_[:,-1])

def test_lasso_lars_vs_lasso_cd(verbose=False):
    """
    Test that LassoLars and Lasso using coordinate descent give the
    same results
    """
    lasso_lars = LassoLARS(alpha=0.1, normalize=False)
    lasso = Lasso(alpha=0.1, fit_intercept=False)
    for alpha in [0.1, 0.01, 0.004]:
        lasso_lars.alpha = alpha
        lasso_lars.fit(X, y, max_iter=12)
        lasso.alpha = alpha
        lasso.fit(X, y, maxit=5000, tol=1e-13)

        # make sure results are the same than with Lasso Coordinate descent
        error = np.linalg.norm(lasso_lars.coef_ - lasso.coef_)
        if verbose:
            print 'Error : ', error
        assert error < 1e-5



################################################################################
# Test LeastAngleRegression
################################################################################

def test_toy():
    """Very simple test on a small dataset"""
    X = [[1, 0, -1.],
         [0, 0, 0],
         [0, 1, .9]]
    Y = [1, 0, -1]

    clf = LeastAngleRegression().fit(X, Y, max_features=1)
    assert_array_almost_equal(clf.coef_, [0, 0, -1.2624], decimal=4)
    assert_array_almost_equal(clf.alphas_, [1.4135, 0.1510], decimal=4)
    assert_array_almost_equal(clf.alphas_.shape, clf.coef_path_.shape[1])
    assert np.linalg.norm(clf.predict(X) - Y) < 0.3
    

    clf = LeastAngleRegression().fit(X, Y) # implicitly max_features=2
    assert_array_almost_equal(clf.coef_, [0, -.0816, -1.34412], decimal=4)
    assert_array_almost_equal(clf.alphas_, \
                 [ 1.41356,   .1510,   3.919e-16], decimal=4)
    assert_array_almost_equal(clf.alphas_.shape, clf.coef_path_.shape[1])
    assert_array_almost_equal(clf.predict(X), np.array(Y))

    # TODO: check that Lasso with coordinate descent finds the same
    # coefficients


def test_toy2():
    """
    Check that LAR responds well to features (columns) that are zero
    arrays.
    """
    X = [[0, 0],
         [0, 1]]
    Y = [0, 1]
    clf = LeastAngleRegression().fit(X, Y)
    assert_array_almost_equal(clf.coef_, [0., .70], decimal=2)


def test_feature_selection():
    n_samples, n_features = 442, 100

    # deterministic test
    np.random.seed(0)

    # generate random input set
    X = np.random.randn(n_samples, n_features)

    # generate a ground truth model with only the first 10 features being non
    # zeros (the other features are not correlated to Y and should be ignored by
    # the L1 regularizer)
    coef_ = np.random.randn(n_features)
    coef_[10:] = 0.0

    # generate the grand truth Y from the model and Y
    Y = np.dot(X, coef_)
    Y += np.random.normal(Y.shape) # make some (label) noise!

    # fit the model assuming that will allready know that only 10 out of
    # n_features are contributing to Y
    clf = LeastAngleRegression().fit(X, Y, max_features=10)

    # ensure that only the first 10 coefs are non zeros and the remaining set to
    # null, as in the ground thrutg model
    assert_equal((clf.coef_[:10] == 0.0).sum(), 0)
    assert_equal((clf.coef_[10:] != 0.0).sum(), 0)

    # train again, but this time without knowing in advance how many features
    # are useful:
    clf = LeastAngleRegression().fit(X, Y, max_features=None)
    assert_equal((clf.coef_[:10] == 0.0).sum(), 0)
    assert_equal((clf.coef_[10:] != 0.0).sum(), 89)

    # explicitly set to zero parameters really close to zero (manual
    # thresholding)
    sparse_coef = np.where(abs(clf.coef_) < 1e-13,
                           np.zeros(clf.coef_.shape),
                           clf.coef_)

    # we find again that only the first 10 features are useful
    assert_equal((sparse_coef[:10] == 0.0).sum(), 0)
    assert_equal((sparse_coef[10:] != 0.0).sum(), 0)


def test_predict(verbose=False):
    """
    Just see if predicted values are close from known response.
    """
    n, m = 10, 20
    np.random.seed(0)
    X = np.random.randn(n, m)
    Y = np.random.randn(n)
    Y = Y - Y.mean() # center response

    Y_ = LeastAngleRegression().fit(X, Y, fit_intercept=False).predict(X)
    if verbose:
        print np.linalg.norm(Y - Y_)
    assert np.linalg.norm(Y-Y_) < 1e-10


    # the same but with an intercept
    Y = Y + 10.
    Y_ = LeastAngleRegression().fit(X, Y).predict(X)
    assert np.linalg.norm(Y-Y_) < 1e-10


def test_sparse_coding():
    """Use LARS as a sparse encoder w.r.t to a given fixed dictionary"""
    n_samples, n_features = 42, 100
    n_dictionary = 30
    max_features = 5

    # generate random input set
    X = np.random.randn(n_samples, n_features)

    # take the first vectors of the dataset as a dictionnary
    D = X[:n_dictionary].T

    assert_equal(D.shape, (n_features, n_dictionary))

    def sparse_encode(vector):
        return LeastAngleRegression().fit(
            D, vector, max_features=max_features, normalize=False,
            fit_intercept=False
        ).coef_

    def sparse_decode(vector):
        return np.dot(D, vector)

    # sparse encode each vector and check the quality of the encoded
    # representation
    for i, x_i in enumerate(X):
        # compute a sparse representation of x_i using the dictionary D
        c_i = sparse_encode(x_i)

        assert_equal(c_i.shape, (n_dictionary,))

        if i < n_dictionary:
            # x_i is one of the vector of D hence there is a trivial sparse code
            expected_code = np.zeros(n_dictionary)
            expected_code[i] = 1.0

            assert_almost_equal(c_i, expected_code, decimal=2)
        else:
            # x_i does not exactly belong to the dictionary, hence up to 5
            # components of the dictionary are used to represent it
            assert_true((c_i != 0.0).sum() <= max_features)



