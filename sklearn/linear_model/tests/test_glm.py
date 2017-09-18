import numpy as np

from sklearn.linear_model.glm import (
    Link,
    IdentityLink,
    LogLink,
    TweedieDistribution,
    NormalDistribution, PoissonDistribution,
    GammaDistribution, InverseGaussianDistribution,
    GeneralizedHyperbolicSecand,
    GeneralizedLinearRegressor)
from sklearn.linear_model.ridge import Ridge

from sklearn.utils.testing import (
    assert_equal, assert_almost_equal,
    assert_array_equal, assert_array_almost_equal)


def test_link_properties():
    """Test link inverse and derivative
    """
    rng = np.random.RandomState(0)
    x = rng.rand(100)*100
    from sklearn.linear_model.glm import Link
    for link in vars()['Link'].__subclasses__():
        link = link()
        assert_almost_equal(link.link(link.inverse(x)), x, decimal=10)
        assert_almost_equal(link.inverse_derivative(link.link(x)),
                            1/link.derivative(x), decimal=10)


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


def test_deviance_zero():
    """Test deviance(y,y) = 0 for different families
    """
    for family in [NormalDistribution(), PoissonDistribution(),
                   GammaDistribution(), InverseGaussianDistribution(),
                   TweedieDistribution(power=-2.5),
                   TweedieDistribution(power=-1),
                   TweedieDistribution(power=1.5),
                   TweedieDistribution(power=2.5),
                   TweedieDistribution(power=4),
                   GeneralizedHyperbolicSecand()]:
        assert_almost_equal(family.deviance(0.1, 0.1), 0, decimal=10)
        assert_almost_equal(family.deviance(1.5, 1.5), 0, decimal=10)


def test_glm_identiy_regression():
    """Test GLM regression with identity link on a simple dataset
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
                alpha=0, family=family, fit_intercept=False, solver=solver)
            res = glm.fit(X, y)
            assert_array_almost_equal(res.coef_, coef)


def test_glm_log_regression():
    """Test GLM regression with log link on a simple dataset
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
                alpha=0, family=family, link=LogLink(), fit_intercept=False,
                solver=solver, start_params='least_squares')
            res = glm.fit(X, y)
            assert_array_almost_equal(res.coef_, coef)


def test_normal_ridge():
    """Test ridge regression for Normal distributions

    Compare to test_ridge in test_ridge.py.
    """
    rng = np.random.RandomState(0)
    alpha = 1.0

    # With more samples than features
    n_samples, n_features, n_predict = 6, 5, 10
    y = rng.randn(n_samples)
    X = rng.randn(n_samples, n_features)
    T = rng.randn(n_predict, n_features)

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=True)
    ridge.fit(X, y)
    for solver in ['irls', 'lbfgs', 'newton-cg']:
        glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0,
                                         family='normal', link='identity',
                                         fit_intercept=True, solver=solver)
        glm.fit(X, y)
        assert_equal(glm.coef_.shape, (X.shape[1], ))
        assert_array_almost_equal(glm.coef_, ridge.coef_)
        assert_almost_equal(glm.intercept_, ridge.intercept_)
        assert_array_almost_equal(glm.predict(T), ridge.predict(T))

    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=False, normalize=False)
    ridge.fit(X, y)
    glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0,
                                     family='normal', link='identity',
                                     fit_intercept=False, solver='irls')
    glm.fit(X, y)
    assert_equal(glm.coef_.shape, (X.shape[1], ))
    assert_array_almost_equal(glm.coef_, ridge.coef_)
    assert_almost_equal(glm.intercept_, ridge.intercept_)
    assert_array_almost_equal(glm.predict(T), ridge.predict(T))

    # With more features than samples
    n_samples, n_features, n_predict = 5, 10, 10
    y = rng.randn(n_samples)
    X = rng.randn(n_samples, n_features)
    T = rng.randn(n_predict, n_features)

    # GLM has 1/(2*n) * Loss + 1/2*L2, Ridge has Loss + L2
    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=True)
    ridge.fit(X, y)
    for solver in ['irls', 'lbfgs', 'newton-cg']:
        glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0,
                                         family='normal', link='identity',
                                         fit_intercept=True, solver=solver)
        glm.fit(X, y)
        assert_equal(glm.coef_.shape, (X.shape[1], ))
        assert_array_almost_equal(glm.coef_, ridge.coef_)
        assert_almost_equal(glm.intercept_, ridge.intercept_)
        assert_array_almost_equal(glm.predict(T), ridge.predict(T))

    ridge = Ridge(alpha=alpha*n_samples, fit_intercept=False, normalize=False)
    ridge.fit(X, y)
    glm = GeneralizedLinearRegressor(alpha=1.0, l1_ratio=0,
                                     family='normal', link='identity',
                                     fit_intercept=False, solver='irls')
    glm.fit(X, y)
    assert_equal(glm.coef_.shape, (X.shape[1], ))
    assert_array_almost_equal(glm.coef_, ridge.coef_)
    assert_almost_equal(glm.intercept_, ridge.intercept_)
    assert_array_almost_equal(glm.predict(T), ridge.predict(T))


# TODO: Test compatibility with R's glm, glmnet
