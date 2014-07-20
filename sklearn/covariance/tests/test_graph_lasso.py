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

    for alpha in (0., .01, .1):
        covs = dict()
        for method in ('cd', 'lars'):
            cov_, _, costs = graph_lasso(emp_cov, alpha=alpha,
                                         return_costs=True)
            covs[method] = cov_
            costs, dual_gap = np.array(costs).T
            # Check that the costs always decrease (doesn't hold if alpha == 0)
            if not alpha == 0:
                assert_array_less(np.diff(costs), 0)
        # Check that the 2 approaches give similar results
        assert_array_almost_equal(covs['cd'], covs['lars'])

    # Smoke test the estimator
    model = GraphLasso(alpha=.1).fit(X)
    model.score(X)
    assert_array_almost_equal(model.covariance_, covs['cd'])
    assert_array_almost_equal(model.covariance_, covs['lars'])

    # For a centered matrix, assume_centered could be chosen True or False
    # Check that this returns indeed the same result for centered data
    Z = X - X.mean(0)
    precs = list()
    for assume_centered in (False, True):
        prec_ = GraphLasso(assume_centered=assume_centered).fit(Z).precision_
        precs.append(prec_)
    assert_array_almost_equal(precs[0], precs[1])


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
        GraphLassoCV(verbose=100, alphas=3, tol=1e-1).fit(X)
    finally:
        sys.stdout = orig_stdout

    # Smoke test with specified alphas
    GraphLassoCV(alphas=[0.8, 0.5], tol=1e-1, n_jobs=1).fit(X)
