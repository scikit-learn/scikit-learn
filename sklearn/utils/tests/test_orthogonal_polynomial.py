# import statements
import numpy as np
import pytest

# import distributions
from scipy.stats import beta, expon, maxwell, norm, uniform

# sklearn imports
from sklearn.utils import check_random_state
from sklearn.utils._orthogonal_polynomial import (
    Hermite,
    Jacobi,
    Laguerre,
    Legendre,
    Polynomial,
)


# check constructors
@pytest.mark.parametrize(
    "polynomial, nargs, error_type",
    [
        (Hermite, 0, TypeError),
        (Jacobi, 2, ValueError),
        (Laguerre, 0, TypeError),
        (Legendre, 0, TypeError),
    ],
)
def test_constructors(polynomial, nargs, error_type):
    (
        # pass
        polynomial(*[1] * nargs) if nargs > 0 else polynomial()
    )

    # constructor with arguments throws error
    with pytest.raises(error_type):
        polynomial(1)


# check vandermonde matrices
@pytest.mark.parametrize(
    "polynomial, points, degree, expected",
    [
        (
            Hermite(),
            [-2, 0, 2],
            2,
            np.array([[1.0, -2.0, 3.0], [1.0, 0.0, -1.0], [1.0, 2.0, 3.0]]),
        ),
        (
            Jacobi(2, 2),
            [-1, 0, 1],
            2,
            np.array([[1, -3, 6], [1, 0, -1], [1, 3, 6]]),
        ),
        (
            Laguerre(),
            [0, 1, 2],
            2,
            np.array([[1, 1, 1], [1, 0, -0.5], [1, -1, -1]]),
        ),
        (
            Legendre(),
            [-1, 0, 1],
            2,
            np.array([[1, -1, 1], [1, 0, -0.5], [1, 1, 1]]),
        ),
    ],
)
def test_vandermonde(polynomial, points, degree, expected):
    # generate vandermonde matrices
    V = polynomial.vandermonde(points, degree)

    # check
    assert np.all(V == expected)


# check vandermonde points
def test_vandermonde_points():
    # non-numeric input throws error
    with pytest.raises(ValueError):
        Legendre().vandermonde("wololo", 2)

    # 2d points throws error
    with pytest.raises(ValueError):
        points = np.array([[-1, 0, 1], [-0.5, 0, 0.5]])
        Legendre().vandermonde(points, 2)

    # passes
    Legendre().vandermonde((-1 / 3, 2 / 5, 7 / 8), 2)


# check vandermonde degree
def test_vandermonde_degree():
    # non-integer degree throws error
    with pytest.raises(ValueError, match="wololo"):
        Legendre().vandermonde([-1, 0, 1], "wololo")

    # degree < 0 throws error
    with pytest.raises(ValueError, match="-1"):
        Legendre().vandermonde([-1, 0, 1], -1)

    # degree = 0 passes
    Legendre().vandermonde([-1, 0, 1], 0)


# test lookup from distribution
@pytest.mark.parametrize(
    "distribution, polynomial_type",
    [
        (norm(), Hermite),
        (beta(2, 2), Jacobi),
        (expon(), Laguerre),
        (uniform(), Legendre),
    ],
)
def test_from_distribution(distribution, polynomial_type):
    # unfrozen distribution type throws error
    with pytest.raises(ValueError, match="dist"):
        Polynomial.from_distribution(uniform)

    # unknown distribution type throws error
    with pytest.raises(ValueError, match="type"):
        Polynomial.from_distribution(maxwell())

    # passes
    assert isinstance(Polynomial.from_distribution(distribution), polynomial_type)


# test scale features from distribution
@pytest.mark.parametrize(
    "distribution, mean, std",
    [
        (norm(loc=1, scale=0.1), 0, 1),
        (beta(2, 2, loc=3, scale=3), 0, np.sqrt(20) / 10),
        (expon(loc=-3, scale=10), 1, 1),
        (uniform(loc=7.5, scale=1), 0, np.sqrt(3) / 3),
    ],
)
def test_scale_features_from_distribution(distribution, mean, std):
    # passes
    n = 100_000
    random_state = check_random_state(17)
    X_1d = distribution.rvs((n, 1), random_state=random_state)
    X_2d = distribution.rvs((n, 2), random_state=random_state)
    polynomial = Polynomial.from_distribution(distribution)
    X_1d_trans = polynomial.scale_features_from_distribution(X_1d, distribution)
    X_2d_trans = polynomial.scale_features_from_distribution(X_2d, distribution)
    assert abs(X_1d_trans.mean() - mean) < 0.01
    assert abs(X_1d_trans.std() - std) < 0.01
    assert abs(X_2d_trans.mean() - mean) < 0.01
    assert abs(X_2d_trans.std() - std) < 0.01

    # unfrozen distribution type throws error
    with pytest.raises(ValueError, match="dist"):
        polynomial.scale_features_from_distribution(X_1d, uniform)

    # 3d array throws error
    with pytest.raises(ValueError, match="2d array"):
        X_3d = distribution.rvs((n, 2, 2), random_state=random_state)
        polynomial.scale_features_from_distribution(X_3d, distribution)

    # unmatched distribution type
    with pytest.raises(ValueError, match="maxwell"):
        polynomial.scale_features_from_distribution(maxwell.rvs((n, 1)), maxwell())


# test norm
@pytest.mark.parametrize(
    "polynomial, norm",
    [
        (Hermite(), 1.4142135623730951),
        (Jacobi(2, 2), 1.1547005383792515),
        (Laguerre(), 1),
        (Legendre(), 0.4472135954999579),
    ],
)
def test_norm(polynomial, norm):
    # non-integer degree throws error
    with pytest.raises(ValueError, match="wololo"):
        polynomial.norm("wololo")

    # degree < 0 throws error
    with pytest.raises(ValueError, match="degree"):
        polynomial.norm(-1)

    # test degree 2
    assert abs(polynomial.norm(2) - norm) < 1e-15


# compare analytical norm to numerically computed values
@pytest.mark.parametrize(
    "distribution, polynomial, tolerance",
    [
        (uniform(-1, 2), Legendre(), 1e-3),
        (norm(), Hermite(), 2e-2),
        (expon(), Laguerre(), 1e-1),
        (beta(2, 5, loc=-1, scale=2), Jacobi(shape_params=(2, 5)), 1e-3),
        (beta(2, 2, loc=-1, scale=2), Jacobi(shape_params=(2, 2)), 1e-3),
        (beta(2, 8, loc=-1, scale=2), Jacobi(shape_params=(2, 8)), 1e-3),
        (beta(5, 5, loc=-1, scale=2), Jacobi(shape_params=(5, 5)), 1e-3),
    ],
)
def test_norm_numerically(distribution, polynomial, tolerance):
    N = 2**22
    degree = 3
    random_state = check_random_state(17)
    X = distribution.rvs((N, 1), random_state=random_state)
    V = polynomial.vandermonde(X, degree)
    for j in range(degree + 1):
        normsq = polynomial._norm_squared(j)
        assert abs(normsq - np.mean(V[:, j] ** 2)) / normsq < tolerance


# special test for Jacobi
def test_jacobi():
    # negative alpha raises error
    with pytest.raises(ValueError, match="alpha"):
        Jacobi(0, 1)

    # negative beta raises error
    with pytest.raises(ValueError, match="beta"):
        Jacobi(1, 0)

    # negative alpha raises error
    with pytest.raises(ValueError, match="alpha"):
        Jacobi(shape_params=(2, 1))

    # negative beta raises error
    with pytest.raises(ValueError, match="beta"):
        Jacobi(shape_params=(1, 2))

    # wrong number of shape parameters throws error
    with pytest.raises(ValueError, match="2 parameters"):
        Jacobi(shape_params=(2,))

    # wrong argument type for alpha throws error
    with pytest.raises(ValueError, match="float or int"):
        Jacobi("a", 1)

    # wrong argument type for beta throws error
    with pytest.raises(ValueError, match="float or int"):
        Jacobi(1, "a")

    # specifying alpha, beta and shape parameters throws error
    with pytest.raises(ValueError, match="both"):
        Jacobi(0.5, 0.5, shape_params=(2, 1))


# test print method
def test_print():
    # default
    polynomial = Legendre()
    assert "Legendre" in str(polynomial)

    # Jacobi
    polynomial = Jacobi(2, 5)
    assert "Jacobi" in str(polynomial)
    assert "2" in str(polynomial)
    assert "5" in str(polynomial)
