""" Test the graph_lasso module.
"""
import sys

import numpy as np
from scipy import linalg

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_less

from sklearn.covariance import (graph_lasso, GraphLasso, GraphLassoCV,
                                empirical_covariance)
from sklearn.datasets.samples_generator import make_sparse_spd_matrix
from sklearn.externals.six.moves import StringIO
from sklearn.utils import check_random_state
from sklearn import datasets


def test_graph_lasso(random_state=0):
    # Sample data from a sparse multivariate normal
    dim = 20
    n_samples = 100
    random_state = check_random_state(random_state)
    prec = make_sparse_spd_matrix(dim, alpha=.95,
                                  random_state=random_state)
    cov = linalg.inv(prec)
    X = random_state.multivariate_normal(np.zeros(dim), cov, size=n_samples)
    emp_cov = empirical_covariance(X)

    for alpha in (0., .1, .25):
        covs = dict()
        icovs = dict()
        for method in ('cd', 'lars'):
            cov_, icov_, costs = graph_lasso(emp_cov, alpha=alpha, mode=method,
                                             return_costs=True)
            covs[method] = cov_
            icovs[method] = icov_
            costs, dual_gap = np.array(costs).T
            # Check that the costs always decrease (doesn't hold if alpha == 0)
            if not alpha == 0:
                assert_array_less(np.diff(costs), 0)
        # Check that the 2 approaches give similar results
        assert_array_almost_equal(covs['cd'], covs['lars'], decimal=4)
        assert_array_almost_equal(icovs['cd'], icovs['lars'], decimal=4)

    # Smoke test the estimator
    model = GraphLasso(alpha=.25).fit(X)
    model.score(X)
    assert_array_almost_equal(model.covariance_, covs['cd'], decimal=4)
    assert_array_almost_equal(model.covariance_, covs['lars'], decimal=4)

    # For a centered matrix, assume_centered could be chosen True or False
    # Check that this returns indeed the same result for centered data
    Z = X - X.mean(0)
    precs = list()
    for assume_centered in (False, True):
        prec_ = GraphLasso(assume_centered=assume_centered).fit(Z).precision_
        precs.append(prec_)
    assert_array_almost_equal(precs[0], precs[1])


def test_graph_lasso_iris():
    # Hard-coded solution from R glasso package for alpha=1.0
    # The iris datasets in R and scikit-learn do not match in a few places,
    # these values are for the scikit-learn version.
    cov_R = np.array([
        [0.68112222, 0.0, 0.2651911, 0.02467558],
        [0.00, 0.1867507, 0.0, 0.00],
        [0.26519111, 0.0, 3.0924249, 0.28774489],
        [0.02467558, 0.0, 0.2877449, 0.57853156]
        ])
    icov_R = np.array([
        [1.5188780, 0.0, -0.1302515, 0.0],
        [0.0, 5.354733, 0.0, 0.0],
        [-0.1302515, 0.0, 0.3502322, -0.1686399],
        [0.0, 0.0, -0.1686399, 1.8123908]
        ])
    X = datasets.load_iris().data
    emp_cov = empirical_covariance(X)
    for method in ('cd', 'lars'):
        cov, icov = graph_lasso(emp_cov, alpha=1.0, return_costs=False,
                                mode=method)
        assert_array_almost_equal(cov, cov_R)
        assert_array_almost_equal(icov, icov_R)


def test_graph_lasso_iris_singular():
    # Small subset of rows to test the rank-deficient case
    # Need to choose samples such that none of the variances are zero
    indices = np.arange(10, 13)

    # Hard-coded solution from R glasso package for alpha=0.01
    cov_R = np.array([
        [0.08, 0.056666662595, 0.00229729713223, 0.00153153142149],
        [0.056666662595, 0.082222222222, 0.00333333333333, 0.00222222222222],
        [0.002297297132, 0.003333333333, 0.00666666666667, 0.00009009009009],
        [0.001531531421, 0.002222222222, 0.00009009009009, 0.00222222222222]
    ])
    icov_R = np.array([
        [24.42244057, -16.831679593, 0.0, 0.0],
        [-16.83168201, 24.351841681, -6.206896552, -12.5],
        [0.0, -6.206896171, 153.103448276, 0.0],
        [0.0, -12.499999143, 0.0, 462.5]
    ])
    X = datasets.load_iris().data[indices, :]
    emp_cov = empirical_covariance(X)
    for method in ('cd', 'lars'):
        cov, icov = graph_lasso(emp_cov, alpha=0.01, return_costs=False,
                                mode=method)
        assert_array_almost_equal(cov, cov_R, decimal=5)
        assert_array_almost_equal(icov, icov_R, decimal=5)


def test_graph_lasso_cv(random_state=1):
    # Sample data from a sparse multivariate normal
    dim = 5
    n_samples = 6
    random_state = check_random_state(random_state)
    prec = make_sparse_spd_matrix(dim, alpha=.96,
                                  random_state=random_state)
    cov = linalg.inv(prec)
    X = random_state.multivariate_normal(np.zeros(dim), cov, size=n_samples)
    # Capture stdout, to smoke test the verbose mode
    orig_stdout = sys.stdout
    try:
        sys.stdout = StringIO()
        # We need verbose very high so that Parallel prints on stdout
        GraphLassoCV(verbose=100, alphas=5, tol=1e-1).fit(X)
    finally:
        sys.stdout = orig_stdout

    # Smoke test with specified alphas
    GraphLassoCV(alphas=[0.8, 0.5], tol=1e-1, n_jobs=1).fit(X)
