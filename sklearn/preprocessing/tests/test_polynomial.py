import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest

from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import KBinsDiscretizer, SplineTransformer


# TODO: add PolynomialFeatures if it moves to _polynomial.py
@pytest.mark.parametrize("est", (SplineTransformer,))
def test_polynomial_and_spline_array_order(est):
    """Test that output array has the given order."""
    X = np.arange(10).reshape(5, 2)

    def is_c_contiguous(a):
        return np.isfortran(a.T)

    assert is_c_contiguous(est().fit_transform(X))
    assert is_c_contiguous(est(order="C").fit_transform(X))
    assert np.isfortran(est(order="F").fit_transform(X))


@pytest.mark.parametrize(
    "params, err_msg",
    [
        ({"degree": -1}, "degree must be a non-negative integer."),
        ({"degree": 2.5}, "degree must be a non-negative integer."),
        ({"degree": "string"}, "degree must be a non-negative integer."),
        ({"n_knots": 1}, "n_knots must be a positive integer >= 2."),
        ({"n_knots": 1}, "n_knots must be a positive integer >= 2."),
        ({"n_knots": 2.5}, "n_knots must be a positive integer >= 2."),
        ({"n_knots": "string"}, "n_knots must be a positive integer >= 2."),
        ({"knots": "string"}, "Expected 2D array, got scalar array instead:"),
        ({"knots": [1, 2]}, "Expected 2D array, got 1D array instead:"),
        (
            {"knots": [[1]]},
            r"Number of knots, knots.shape\[0\], must be >= 2.",
        ),
        (
            {"knots": [[1, 5], [2, 6]]},
            r"knots.shape\[1\] == n_features is violated.",
        ),
        (
            {"knots": [[1], [1], [2]]},
            "knots must be sorted without duplicates.",
        ),
        ({"knots": [[2], [1]]}, "knots must be sorted without duplicates."),
        (
            {"extrapolation": None},
            "extrapolation must be one of 'error', 'constant', 'linear' or "
            "'continue'.",
        ),
        (
            {"extrapolation": 1},
            "extrapolation must be one of 'error', 'constant', 'linear' or "
            "'continue'.",
        ),
        (
            {"extrapolation": "string"},
            "extrapolation must be one of 'error', 'constant', 'linear' or "
            "'continue'.",
        ),
        ({"include_bias": None}, "include_bias must be bool."),
        ({"include_bias": 1}, "include_bias must be bool."),
        ({"include_bias": "string"}, "include_bias must be bool."),
    ],
)
def test_spline_transformer_input_validation(params, err_msg):
    """Test that we raise errors for invalid input in SplineTransformer."""
    X = [[1], [2]]

    with pytest.raises(ValueError, match=err_msg):
        SplineTransformer(**params).fit(X)


def test_spline_transformer_manual_knot_input():
    """Test that array-like knot positions in SplineTransformer are accepted.
    """
    X = np.arange(20).reshape(10, 2)
    knots = [[0.5, 1], [1.5, 2], [5, 10]]
    st1 = SplineTransformer(degree=3, knots=knots).fit(X)
    knots = np.asarray(knots)
    st2 = SplineTransformer(degree=3, knots=knots).fit(X)
    for i in range(X.shape[1]):
        assert_allclose(st1.bsplines_[i].t, st2.bsplines_[i].t)


def test_spline_transformer_feature_names():
    """Test that SplineTransformer generates correct features name."""
    X = np.arange(20).reshape(10, 2)
    splt = SplineTransformer(n_knots=3, degree=3, include_bias=True).fit(X)
    feature_names = splt.get_feature_names()
    assert_array_equal(
        feature_names,
        [
            "x0_sp_0",
            "x0_sp_1",
            "x0_sp_2",
            "x0_sp_3",
            "x0_sp_4",
            "x1_sp_0",
            "x1_sp_1",
            "x1_sp_2",
            "x1_sp_3",
            "x1_sp_4",
        ],
    )

    splt = SplineTransformer(n_knots=3, degree=3, include_bias=False).fit(X)
    feature_names = splt.get_feature_names(["a", "b"])
    assert_array_equal(
        feature_names,
        [
            "a_sp_0",
            "a_sp_1",
            "a_sp_2",
            "a_sp_3",
            "b_sp_0",
            "b_sp_1",
            "b_sp_2",
            "b_sp_3",
        ],
    )


@pytest.mark.parametrize("degree", range(1, 5))
@pytest.mark.parametrize("n_knots", range(3, 5))
@pytest.mark.parametrize("knots", ["uniform", "quantile"])
def test_spline_transformer_unity_decomposition(degree, n_knots, knots):
    """Test that B-splines are indeed a decomposition of unity.

    Splines basis functions must sum up to 1 per row, if we stay in between
    boundaries.
    """
    X = np.linspace(0, 1, 100)[:, None]
    # make the boundaries 0 and 1 part of X_train, for sure.
    X_train = np.r_[[[0]], X[::2, :], [[1]]]
    X_test = X[1::2, :]
    splt = SplineTransformer(
        n_knots=n_knots, degree=degree, knots=knots, include_bias=True
    )
    splt.fit(X_train)
    for X in [X_train, X_test]:
        assert_allclose(np.sum(splt.transform(X), axis=1), 1)


@pytest.mark.parametrize(["bias", "intercept"], [(True, False), (False, True)])
def test_spline_transformer_linear_regression(bias, intercept):
    """Test that B-splines fit a sinusodial curve pretty well."""
    X = np.linspace(0, 10, 100)[:, None]
    y = np.sin(X[:, 0]) + 2  # +2 to avoid the value 0 in assert_allclose
    pipe = Pipeline(
        steps=[
            (
                "spline",
                SplineTransformer(
                    n_knots=15,
                    degree=3,
                    include_bias=bias,
                    extrapolation="constant",
                ),
            ),
            ("ols", LinearRegression(fit_intercept=intercept)),
        ]
    )
    pipe.fit(X, y)
    assert_allclose(pipe.predict(X), y, rtol=1e-3)


@pytest.mark.parametrize(["bias", "intercept"], [(True, False), (False, True)])
@pytest.mark.parametrize("degree", [1, 2, 3, 4, 5])
def test_spline_transformer_extrapolation(bias, intercept, degree):
    """Test that B-spline extrapolation works correctly."""
    # we use a straight line for that
    X = np.linspace(-1, 1, 100)[:, None]
    y = X.squeeze()

    # 'constant'
    pipe = Pipeline(
        [
            [
                "spline",
                SplineTransformer(
                    n_knots=4,
                    degree=degree,
                    include_bias=bias,
                    extrapolation="constant",
                ),
            ],
            ["ols", LinearRegression(fit_intercept=intercept)],
        ]
    )
    pipe.fit(X, y)
    assert_allclose(pipe.predict([[-10], [5]]), [-1, 1])

    # 'linear'
    pipe = Pipeline(
        [
            [
                "spline",
                SplineTransformer(
                    n_knots=4,
                    degree=degree,
                    include_bias=bias,
                    extrapolation="linear",
                ),
            ],
            ["ols", LinearRegression(fit_intercept=intercept)],
        ]
    )
    pipe.fit(X, y)
    assert_allclose(pipe.predict([[-10], [5]]), [-10, 5])

    # 'error'
    splt = SplineTransformer(
        n_knots=4, degree=degree, include_bias=bias, extrapolation="error"
    )
    splt.fit(X)
    with pytest.raises(ValueError):
        splt.transform([[-10]])
    with pytest.raises(ValueError):
        splt.transform([[5]])


def test_spline_transformer_kbindiscretizer():
    """Test that a B-spline of degree=0 is equivalent to KBinsDiscretizer."""
    rng = np.random.RandomState(97531)
    X = rng.randn(200).reshape(200, 1)
    n_bins = 5
    n_knots = n_bins + 1

    splt = SplineTransformer(
        n_knots=n_knots, degree=0, knots="quantile", include_bias=True
    )
    splines = splt.fit_transform(X)

    kbd = KBinsDiscretizer(
        n_bins=n_bins, encode="onehot-dense", strategy="quantile"
    )
    kbins = kbd.fit_transform(X)

    # Though they should be exactly equal, we test approximately with high
    # accuracy.
    assert_allclose(splines, kbins, rtol=1e-13)


@pytest.mark.parametrize("n_knots", [5, 10])
@pytest.mark.parametrize("include_bias", [True, False])
@pytest.mark.parametrize("degree", [3, 5])
def test_spline_transformer_n_features_out(n_knots, include_bias, degree):
    """Test that transform results in n_features_out_ features."""
    splt = SplineTransformer(
        n_knots=n_knots,
        degree=degree,
        include_bias=include_bias
    )
    X = np.linspace(0, 1, 10)[:, None]
    splt.fit(X)

    assert splt.transform(X).shape[1] == splt.n_features_out_
