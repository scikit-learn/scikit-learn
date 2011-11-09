""" Test the g_lasso module.
"""
import numpy as np
from scipy import linalg

from sklearn.covariance import graph_lasso, GraphLasso, GraphLassoCV
from sklearn.datasets.samples_generator import make_sparse_spd_matrix
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

    for alpha in (.1, .01):
        covs = dict()
        for method in ('cd', 'lars'):
            cov_, _, costs = graph_lasso(X, alpha=.1, return_costs=True)
            covs[method] = cov_
            costs, dual_gap = np.array(costs).T
            # Check that the costs always decrease
            np.testing.assert_array_less(np.diff(costs), 0)
        # Check that the 2 approaches give similar results
        np.testing.assert_allclose(covs['cd'], covs['lars'])

    # Smoke test the estimator
    model = GraphLasso(alpha=.1).fit(X)
    np.testing.assert_allclose(model.covariance_, covs['cd'])


def test_graph_lasso_cv(random_state=1):
    # Sample data from a sparse multivariate normal
    dim = 5
    n_samples = 6
    random_state = check_random_state(random_state)
    prec = make_sparse_spd_matrix(dim, alpha=.96,
                                  random_state=random_state)
    cov = linalg.inv(prec)
    X = random_state.multivariate_normal(np.zeros(dim), cov, size=n_samples)
    GraphLassoCV(verbose=10, alphas=3).fit(X)

