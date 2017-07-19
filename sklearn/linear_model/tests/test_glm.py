import numpy as np

from sklearn.linear_model.glm import (
    # Link, IdentityLink,
    LogLink,
    TweedieDistribution,
    NormalDistribution, PoissonDistribution,
    GammaDistribution, InverseGaussianDistribution,
    # GeneralizedHyperbolicSecand,
    GeneralizedLinearRegressor)

from sklearn.utils.testing import (
    # assert_equal,
    assert_array_equal, assert_array_almost_equal)


def test_family_bounds():
    """Test the valid range of distributions
    """
    family = NormalDistribution()
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [True, True, True])

    family = PoissonDistribution()
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, True, True])

    family = TweedieDistribution(power=1.5)
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, True, True])

    family = GammaDistribution()
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, False, True])

    family = InverseGaussianDistribution()
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, False, True])

    family = TweedieDistribution(power=4.5)
    result = family.in_y_range([-1, 0, 1])
    assert_array_equal(result, [False, False, True])


def test_glm_identiy_regression():
    """Test linear regression on a simple dataset
    """
    coef = [1, 2]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.dot(X, coef)
    families = (
        NormalDistribution(), PoissonDistribution(),
        GammaDistribution(), InverseGaussianDistribution(),
        TweedieDistribution(power=1.5), TweedieDistribution(power=4.5))
    for solver in ['irls', 'lbfgs', 'newton-cg']:
        for family in families:
            glm = GeneralizedLinearRegressor(
                family=family, fit_intercept=False, solver=solver)
            res = glm.fit(X, y)
            assert_array_almost_equal(res.coef_, coef)


def test_glm_log_regression():
    """Test linear regression on a simple dataset
    """
    coef = [1, 2]
    X = np.array([[1, 1, 1, 1, 1], [0, 1, 2, 3, 4]]).T
    y = np.exp(np.dot(X, coef))
    families = (
        NormalDistribution(), PoissonDistribution(),
        GammaDistribution(), InverseGaussianDistribution(),
        TweedieDistribution(power=1.5), TweedieDistribution(power=4.5))
    for solver in ['irls', 'lbfgs', 'newton-cg']:
        for family in families:
            glm = GeneralizedLinearRegressor(
                family=family, link=LogLink(), fit_intercept=False,
                solver=solver, start_params='ols')
            res = glm.fit(X, y)
            assert_array_almost_equal(res.coef_, coef)


# TODO: Test compatibility with R's glm, glmnet
