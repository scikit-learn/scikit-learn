# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest
from scipy.stats import beta, expon, maxwell, norm, uniform

from sklearn.utils import check_random_state
from sklearn.utils._orthogonal_polynomial import (
    Hermite,
    Jacobi,
    Laguerre,
    Legendre,
    Polynomial,
)


# Check constructors
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
        # Pass
        polynomial(*[1] * nargs) if nargs > 0 else polynomial()
    )

    # Constructor with arguments throws error
    with pytest.raises(error_type):
        polynomial(1)


# Check vandermonde matrices
@pytest.mark.parametrize(
    "polynomial, points, degree, expected",
    [
        # Orthogonal polynomials
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
        # Orthonormal polynomials
        (
            Hermite(normalize=True),
            [-2, 0, 2],
            2,
            np.array(
                [
                    [1.0, -2.0, 3.0 / np.sqrt(2)],
                    [1.0, 0.0, -1.0 / np.sqrt(2)],
                    [1.0, 2.0, 3.0 / np.sqrt(2)],
                ]
            ),
        ),
        (
            Jacobi(2, 2, normalize=True),
            [-1, 0, 1],
            2,
            np.array(
                [
                    [1, -np.sqrt(7), np.sqrt(27)],
                    [1, 0, -np.sqrt(3 / 4)],
                    [1, np.sqrt(7), np.sqrt(27)],
                ]
            ),
        ),
        (
            Laguerre(normalize=True),
            [0, 1, 2],
            2,
            np.array([[1, 1, 1], [1, 0, -0.5], [1, -1, -1]]),
        ),
        (
            Legendre(normalize=True),
            [-1, 0, 1],
            2,
            np.array(
                [
                    [1, -np.sqrt(3), np.sqrt(5)],
                    [1, 0, -np.sqrt(5) / 2],
                    [1, np.sqrt(3), np.sqrt(5)],
                ]
            ),
        ),
    ],
)
def test_vandermonde(polynomial, points, degree, expected):
    # Generate vandermonde matrices
    V = polynomial.vandermonde(points, degree)

    # Check
    assert np.all(np.abs(V - expected) < 1e-15)


# Check vandermonde points
def test_vandermonde_points():
    # Non-numeric input throws error
    with pytest.raises(ValueError):
        Legendre().vandermonde("wololo", 2)

    # 2d points throws error
    with pytest.raises(ValueError):
        points = np.array([[-1, 0, 1], [-0.5, 0, 0.5]])
        Legendre().vandermonde(points, 2)

    # Passes
    Legendre().vandermonde((-1 / 3, 2 / 5, 7 / 8), 2)


# Check vandermonde degree
def test_vandermonde_degree():
    # Non-integer degree throws error
    with pytest.raises(ValueError, match="wololo"):
        Legendre().vandermonde([-1, 0, 1], "wololo")

    # Degree < 0 throws error
    with pytest.raises(ValueError, match="-1"):
        Legendre().vandermonde([-1, 0, 1], -1)

    # Degree = 0 passes
    Legendre().vandermonde([-1, 0, 1], 0)


# Test lookup from distribution
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
    # Unfrozen distribution type throws error
    with pytest.raises(ValueError, match="dist"):
        Polynomial.from_distribution(uniform)

    # Unknown distribution type throws error
    with pytest.raises(ValueError, match="type"):
        Polynomial.from_distribution(maxwell())

    # Passes
    assert isinstance(Polynomial.from_distribution(distribution), polynomial_type)


# Test scale features from distribution
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
    # Passes
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

    # Unfrozen distribution type throws error
    with pytest.raises(ValueError, match="dist"):
        polynomial.scale_features_from_distribution(X_1d, uniform)

    # 3d array throws error
    with pytest.raises(ValueError, match="2d array"):
        X_3d = distribution.rvs((n, 2, 2), random_state=random_state)
        polynomial.scale_features_from_distribution(X_3d, distribution)

    # Unmatched distribution type
    with pytest.raises(ValueError, match="maxwell"):
        polynomial.scale_features_from_distribution(maxwell.rvs((n, 1)), maxwell())


# Test norm
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
    # Non-integer degree throws error
    with pytest.raises(ValueError, match="wololo"):
        polynomial.norm("wololo")

    # Degree < 0 throws error
    with pytest.raises(ValueError, match="degree"):
        polynomial.norm(-1)

    # Test degree 2
    assert abs(polynomial.norm(2) - norm) < 1e-15


# Compare analytical norm to numerically computed values
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


# Special test for Jacobi
def test_jacobi():
    # Negative alpha raises error
    with pytest.raises(ValueError, match="alpha"):
        Jacobi(0, 1)

    # Negative beta raises error
    with pytest.raises(ValueError, match="beta"):
        Jacobi(1, 0)

    # Negative alpha raises error
    with pytest.raises(ValueError, match="alpha"):
        Jacobi(shape_params=(2, 1))

    # Negative beta raises error
    with pytest.raises(ValueError, match="beta"):
        Jacobi(shape_params=(1, 2))

    # Wrong number of shape parameters throws error
    with pytest.raises(ValueError, match="2 parameters"):
        Jacobi(shape_params=(2,))

    # Wrong argument type for alpha throws error
    with pytest.raises(ValueError, match="float or int"):
        Jacobi("a", 1)

    # Wrong argument type for beta throws error
    with pytest.raises(ValueError, match="float or int"):
        Jacobi(1, "a")

    # Specifying alpha, beta and shape parameters throws error
    with pytest.raises(ValueError, match="both"):
        Jacobi(0.5, 0.5, shape_params=(2, 1))


# Test print method
def test_print():
    # Default
    polynomial = Legendre()
    assert "Legendre" in str(polynomial)

    # Jacobi
    polynomial = Jacobi(2, 5)
    assert "Jacobi" in str(polynomial)
    assert "2" in str(polynomial)
    assert "5" in str(polynomial)
