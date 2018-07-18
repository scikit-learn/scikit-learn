""" Test the graphical_lasso module.
"""
import sys

import numpy as np
from scipy import linalg
import pytest

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_less
from sklearn.utils.testing import assert_warns_message

from sklearn.covariance import (graphical_lasso, GraphicalLasso,
                                GraphicalLassoCV, empirical_covariance)
from sklearn.datasets.samples_generator import make_sparse_spd_matrix
from sklearn.externals.six.moves import StringIO
from sklearn.utils import check_random_state
from sklearn import datasets
from sklearn.utils.fixes import PY3_OR_LATER

from numpy.testing import assert_equal


def test_graphical_lasso(random_state=0):
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
            cov_, icov_, costs = graphical_lasso(emp_cov, return_costs=True,
                                                 alpha=alpha, mode=method)
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
    model = GraphicalLasso(alpha=.25).fit(X)
    model.score(X)
    assert_array_almost_equal(model.covariance_, covs['cd'], decimal=4)
    assert_array_almost_equal(model.covariance_, covs['lars'], decimal=4)

    # For a centered matrix, assume_centered could be chosen True or False
    # Check that this returns indeed the same result for centered data
    Z = X - X.mean(0)
    precs = list()
    for assume_centered in (False, True):
        prec_ = GraphicalLasso(
            assume_centered=assume_centered).fit(Z).precision_
        precs.append(prec_)
    assert_array_almost_equal(precs[0], precs[1])


def test_graphical_lasso_iris():
    # Hard-coded solution from R glasso package for alpha=1.0
    # (need to set penalize.diagonal to FALSE)
    cov_R = np.array([
        [0.68112222, 0.0000000, 0.265820, 0.02464314],
        [0.00000000, 0.1887129, 0.000000, 0.00000000],
        [0.26582000, 0.0000000, 3.095503, 0.28697200],
        [0.02464314, 0.0000000, 0.286972, 0.57713289]
        ])
    icov_R = np.array([
        [1.5190747, 0.000000, -0.1304475, 0.0000000],
        [0.0000000, 5.299055, 0.0000000, 0.0000000],
        [-0.1304475, 0.000000, 0.3498624, -0.1683946],
        [0.0000000, 0.000000, -0.1683946, 1.8164353]
        ])
    X = datasets.load_iris().data
    emp_cov = empirical_covariance(X)
    for method in ('cd', 'lars'):
        cov, icov = graphical_lasso(emp_cov, alpha=1.0, return_costs=False,
                                    mode=method)
        assert_array_almost_equal(cov, cov_R)
        assert_array_almost_equal(icov, icov_R)


def test_graphical_lasso_iris_singular():
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
        cov, icov = graphical_lasso(emp_cov, alpha=0.01, return_costs=False,
                                    mode=method)
        assert_array_almost_equal(cov, cov_R, decimal=5)
        assert_array_almost_equal(icov, icov_R, decimal=5)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
def test_graphical_lasso_cv(random_state=1):
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
        GraphicalLassoCV(verbose=100, alphas=5, tol=1e-1).fit(X)
    finally:
        sys.stdout = orig_stdout

    # Smoke test with specified alphas
    GraphicalLassoCV(alphas=[0.8, 0.5], tol=1e-1, n_jobs=1).fit(X)


@pytest.mark.filterwarnings('ignore: You should specify a value')  # 0.22
@pytest.mark.skipif(not PY3_OR_LATER,
                    reason='On Python 2 DeprecationWarning is not issued for some unkown reason.')
def test_deprecated_grid_scores(random_state=1):
    dim = 5
    n_samples = 6
    random_state = check_random_state(random_state)
    prec = make_sparse_spd_matrix(dim, alpha=.96,
                                  random_state=random_state)
    cov = linalg.inv(prec)
    X = random_state.multivariate_normal(np.zeros(dim), cov, size=n_samples)
    graphical_lasso = GraphicalLassoCV(alphas=[0.8, 0.5], tol=1e-1, n_jobs=1)
    graphical_lasso.fit(X)

    depr_message = ("Attribute grid_scores was deprecated in version "
                    "0.19 and will be removed in 0.21. Use "
                    "``grid_scores_`` instead")

    with pytest.warns(DeprecationWarning, match=depr_message):
        assert_equal(graphical_lasso.grid_scores, graphical_lasso.grid_scores_)
