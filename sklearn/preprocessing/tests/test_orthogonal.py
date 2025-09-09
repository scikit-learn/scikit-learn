# import statements
import numpy as np
import pytest

# sklearn imports
from sklearn.preprocessing import OrthogonalPolynomialFeatures
from sklearn.utils._multiindexset import TotalDegree


# check fit_transform
def test_fit_transform():
    X = np.linspace(0, 1, num=6).reshape(3, 2)
    poly = OrthogonalPolynomialFeatures()
    X_trans = poly.fit_transform(X)
    X_exact = np.array(
        [
            [1.0, 0.0, -0.5, 0.2, 0.0, -0.44],
            [1.0, 0.4, -0.26, 0.6, 0.24, 0.04],
            [1.0, 0.8, 0.46, 1.0, 0.8, 1.0],
        ]
    )
    assert np.linalg.norm(X_trans - X_exact) < 1e-15


# check fit_transform with different polynomial types
def test_fit_transform_hybrid():
    X = np.linspace(0, 1, num=6).reshape(3, 2)
    poly = OrthogonalPolynomialFeatures(polynomial=("Legendre", "Hermite"))
    X_trans = poly.fit_transform(X)
    X_exact = np.array(
        [
            [1.0, 0.0, -0.5, 0.2, 0.0, -0.96],
            [1.0, 0.4, -0.26, 0.6, 0.24, -0.64],
            [1.0, 0.8, 0.46, 1.0, 0.8, 0.0],
        ]
    )
    assert np.linalg.norm(X_trans - X_exact) < 1e-15


# check polynomial argument
def test_fit_polynomial():
    # unmatched number of input features and number of polynomials throws error
    with pytest.raises(ValueError, match="polynomial"):
        X = np.linspace(0, 1, num=9).reshape(3, 3)
        poly = OrthogonalPolynomialFeatures(polynomial=("Legendre", "Hermite"))
        poly.fit_transform(X)

    # unknown polynomial type throws error
    with pytest.raises(ValueError, match="polynomial"):
        X = np.linspace(0, 1, num=6).reshape(3, 2)
        poly = OrthogonalPolynomialFeatures(polynomial="ab")
        poly.fit_transform(X)

    # unknown polynomial type throws error
    with pytest.raises(ValueError, match="polynomial"):
        X = np.linspace(0, 1, num=6).reshape(3, 2)
        poly = OrthogonalPolynomialFeatures(polynomial=3)
        poly.fit_transform(X)

    # unknown polynomial type throws error
    with pytest.raises(ValueError, match="polynomial"):
        X = np.linspace(0, 1, num=6).reshape(3, 2)
        poly = OrthogonalPolynomialFeatures(polynomial=("Legendre", 3))
        poly.fit_transform(X)


# check degree / truncation specification
@pytest.mark.parametrize(
    "truncation, n_terms",
    [
        ("full_tensor", 512),
        ("total_degree", 120),
        ("hyperbolic_cross", 38),
        ("Zaremba_cross", 98),
    ],
)
def test_fit_predefined_multiindex_set_shape(truncation, n_terms):
    X = np.linspace(0, 1, num=6).reshape(2, 3)
    poly = OrthogonalPolynomialFeatures(7, truncation=truncation)
    X_trans = poly.fit_transform(X)
    assert X_trans.shape[0] == 2
    assert X_trans.shape[1] == n_terms


# check degree / truncation / weights specification
@pytest.mark.parametrize(
    "truncation, n_terms",
    [
        ("full_tensor", 28),
        ("total_degree", 17),
        ("hyperbolic_cross", 11),
        ("Zaremba_cross", 15),
    ],
)
def test_fit_predefined_multiindex_set_shape_weighted(truncation, n_terms):
    X = np.linspace(0, 1, num=6).reshape(3, 2)
    poly = OrthogonalPolynomialFeatures(6, truncation=truncation, weights=(1, 3 / 5))
    X_trans = poly.fit_transform(X)
    assert X_trans.shape[0] == 3
    assert X_trans.shape[1] == n_terms


# check custom multiindices
@pytest.mark.parametrize(
    "multiindices",
    [
        TotalDegree(2, 4).indices(),
        list(TotalDegree(2, 4).indices()),
        np.vstack(list(TotalDegree(2, 4).indices())),
    ],
)
def test_fit_custom_multiindices(multiindices):
    X = np.linspace(0, 1, num=6).reshape(3, 2)
    poly = OrthogonalPolynomialFeatures(6, multiindices=multiindices)
    X_trans = poly.fit_transform(X)
    assert X_trans.shape[0] == 3
    assert X_trans.shape[1] == 15


# check transform
def test_transform():
    with pytest.raises(ValueError, match="fit"):
        X = np.linspace(0, 1, num=6).reshape(3, 2)
        poly = OrthogonalPolynomialFeatures()
        poly.transform(X)


# check feature names
def test_feature_names():
    X = np.linspace(0, 1, num=6).reshape(3, 2)
    poly = OrthogonalPolynomialFeatures()

    # unfitted throws error
    with pytest.raises(ValueError, match="fit"):
        poly.get_feature_names_out(X)

    # passes
    poly.fit_transform(X)

    # check output features
    names = poly.get_feature_names_out()
    assert len(names) == 6
    for name in names:
        assert "Legendre" in name
        assert "x0" in name

    # check output features with custom input feature names
    names = poly.get_feature_names_out(input_features=("a", "b"))
    assert len(names) == 6
    for name in names:
        assert "Legendre" in name
        assert "a" in name
        assert "b" in name

    # check invalid number of input features
    with pytest.raises(ValueError, match="features"):
        poly.get_feature_names_out(input_features=("a", "b", "c"))
