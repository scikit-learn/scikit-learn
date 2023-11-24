# import statements
import numpy as np
import pytest

# qualified import statements
from math import prod
from scipy.stats import uniform

# sklearn imports
from sklearn.polynomial_chaos._adaptive import BasisIncrementStrategy, GerstnerGriebel

from sklearn.linear_model import LassoCV
from sklearn.polynomial_chaos import PolynomialChaosRegressor
from sklearn.utils import check_random_state

# Check from string
def test_from_string():

    # passes
    strategy = BasisIncrementStrategy.from_string("gerstner_griebel")
    assert isinstance(strategy, GerstnerGriebel)

    # unknown basis growth strategy
    with pytest.raises(ValueError, match="unknown strategy"):
        BasisIncrementStrategy.from_string("wololo")

# Check adaptive construction
def test_adaptive_construction():
    random_state = check_random_state(17)
    a = np.array([1, 2, 5, 10, 20, 50, 100, 500])
    dimension = len(a)
    solver = LassoCV(fit_intercept=False, tol=1e-4)
    pce = PolynomialChaosRegressor(uniform(), degree=0, solver=solver)
    N = 66
    X = uniform().rvs((N, dimension), random_state=random_state)
    y = prod((abs(4*X_j - 2) + a_j) / (1 + a_j) for a_j, X_j in zip(a, X.T))
    pce.fit(X, y, max_iter=30)
    exact = np.array([0.6342, 0.2945, 0.0756, 0.0227, 0.0062, 0.0011, 0.0003, 
                      0.0000])
    rel_err = np.linalg.norm(pce.total_sens() - exact) / np.linalg.norm(exact)
    assert rel_err < 0.1

# Check old/active split
def test_old_active_split():
    random_state = check_random_state(17)
    distribution = uniform(0, 1)
    pce = PolynomialChaosRegressor(distribution, degree=1)
    X = distribution.rvs((66, 1), random_state=random_state)
    y = np.prod((3*X**2 + 1)/2, axis=1)
    pce.fit(X, y, max_iter=2)
    assert not np.any(pce.strategy_.old[0]) # (0, 0, ..., 0) is in old set