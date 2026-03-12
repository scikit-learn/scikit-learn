import numpy as np
import pytest
from sklearn.linear_model import POPSRegression, BayesianRidge
from sklearn.utils.estimator_checks import parametrize_with_checks

@parametrize_with_checks([POPSRegression()])
def test_pops_check_estimator(estimator, check):
    check(estimator)

def test_pops_zero_noise_limit():
    # In low-noise limit, BayesianRidge uncertainty vanishes,
    # but POPS should maintain misspecification uncertainty.
    rng = np.random.RandomState(42)
    X = rng.randn(100, 5)
    w = rng.randn(5)
    # Slightly non-linear to induce misspecification
    y = X @ w + 0.1 * X[:, 0]**2 + 1e-6 * rng.randn(100)
    
    pops = POPSRegression().fit(X, y)
    br = BayesianRidge().fit(X, y)
    
    _, y_std_pops = pops.predict(X, return_std=True)
    _, y_std_br = br.predict(X, return_std=True)
    
    # POPS uncertainty should be larger than BayesianRidge in this regime
    assert np.mean(y_std_pops) > np.mean(y_std_br)

def test_pops_return_bounds():
    rng = np.random.RandomState(42)
    X = rng.randn(20, 2)
    y = X @ [1, 2] + 0.1 * rng.randn(20)
    
    pops = POPSRegression().fit(X, y)
    y_mean, y_max, y_min = pops.predict(X, return_bounds=True)
    
    assert y_max.shape == (20,)
    assert y_min.shape == (20,)
    assert np.all(y_max >= y_mean)
    assert np.all(y_min <= y_mean)

def test_pops_posterior_modes():
    rng = np.random.RandomState(42)
    X = rng.randn(20, 2)
    y = X @ [1, 2] + 0.1 * rng.randn(20)
    
    for posterior in ["hypercube", "ensemble"]:
        pops = POPSRegression(posterior=posterior).fit(X, y)
        assert hasattr(pops, "misspecification_sigma_")
