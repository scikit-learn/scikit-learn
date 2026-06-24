import warnings

import numpy as np
import pytest

from sklearn import config_context, datasets
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.compose import TransformedTargetRegressor
from sklearn.datasets import make_regression
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LinearRegression, OrthogonalMatchingPursuit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler
from sklearn.utils._testing import assert_allclose

friedman = datasets.make_friedman1(random_state=0)


def test_transform_target_regressor_error():
    X, y = friedman
    # provide a transformer and functions at the same time
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(),
        transformer=StandardScaler(),
        func=np.exp,
        inverse_func=np.log,
    )
    with pytest.raises(
        ValueError,
        match="'transformer' and functions 'func'/'inverse_func' cannot both be set.",
    ):
        regr.fit(X, y)
    # fit with sample_weight with a regressor which does not support it
    sample_weight = np.ones((y.shape[0],))
    regr = TransformedTargetRegressor(
        regressor=OrthogonalMatchingPursuit(), transformer=StandardScaler()
    )
    with pytest.raises(
        TypeError,
        match=r"fit\(\) got an unexpected keyword argument 'sample_weight'",
    ):
        regr.fit(X, y, sample_weight=sample_weight)

    # one of (func, inverse_func) is given but the other one is not
    regr = TransformedTargetRegressor(func=np.exp)
    with pytest.raises(
        ValueError,
        match="When 'func' is provided, 'inverse_func' must also be provided",
    ):
        regr.fit(X, y)

    regr = TransformedTargetRegressor(inverse_func=np.log)
    with pytest.raises(
        ValueError,
        match="When 'inverse_func' is provided, 'func' must also be provided",
    ):
        regr.fit(X, y)


def test_transform_target_regressor_invertible():
    X, y = friedman
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.sqrt,
        inverse_func=np.log,
        check_inverse=True,
    )
    with pytest.warns(
        UserWarning,
        match=(r"The provided functions.* are not strictly inverse of each other"),
    ):
        regr.fit(X, y)
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(), func=np.sqrt, inverse_func=np.log
    )
    regr.set_params(check_inverse=False)

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        regr.fit(X, y)


def _check_standard_scaled(y, y_pred):
    y_mean = np.mean(y, axis=0)
    y_std = np.std(y, axis=0)
    assert_allclose((y - y_mean) / y_std, y_pred)


def _check_shifted_by_one(y, y_pred):
    assert_allclose(y + 1, y_pred)


def test_transform_target_regressor_functions():
    X, y = friedman
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(), func=np.log, inverse_func=np.exp
    )
    y_pred = regr.fit(X, y).predict(X)
    # check the transformer output
    y_tran = regr.transformer_.transform(y.reshape(-1, 1)).squeeze()
    assert_allclose(np.log(y), y_tran)
    assert_allclose(
        y, regr.transformer_.inverse_transform(y_tran.reshape(-1, 1)).squeeze()
    )
    assert y.shape == y_pred.shape
    assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
    # check the regressor output
    lr = LinearRegression().fit(X, regr.func(y))
    assert_allclose(regr.regressor_.coef_.ravel(), lr.coef_.ravel())


def test_transform_target_regressor_functions_multioutput():
    X = friedman[0]
    y = np.vstack((friedman[1], friedman[1] ** 2 + 1)).T
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(), func=np.log, inverse_func=np.exp
    )
    y_pred = regr.fit(X, y).predict(X)
    # check the transformer output
    y_tran = regr.transformer_.transform(y)
    assert_allclose(np.log(y), y_tran)
    assert_allclose(y, regr.transformer_.inverse_transform(y_tran))
    assert y.shape == y_pred.shape
    assert_allclose(y_pred, regr.inverse_func(regr.regressor_.predict(X)))
    # check the regressor output
    lr = LinearRegression().fit(X, regr.func(y))
    assert_allclose(regr.regressor_.coef_.ravel(), lr.coef_.ravel())


@pytest.mark.parametrize(
    "X,y", [friedman, (friedman[0], np.vstack((friedman[1], friedman[1] ** 2 + 1)).T)]
)
def test_transform_target_regressor_1d_transformer(X, y):
    # All transformer in scikit-learn expect 2D data. FunctionTransformer with
    # validate=False lift this constraint without checking that the input is a
    # 2D vector. We check the consistency of the data shape using a 1D and 2D y
    # array.
    transformer = FunctionTransformer(
        func=lambda x: x + 1, inverse_func=lambda x: x - 1
    )
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(), transformer=transformer
    )
    y_pred = regr.fit(X, y).predict(X)
    assert y.shape == y_pred.shape
    # consistency forward transform
    y_tran = regr.transformer_.transform(y)
    _check_shifted_by_one(y, y_tran)
    assert y.shape == y_pred.shape
    # consistency inverse transform
    assert_allclose(y, regr.transformer_.inverse_transform(y_tran).squeeze())
    # consistency of the regressor
    lr = LinearRegression()
    transformer2 = clone(transformer)
    lr.fit(X, transformer2.fit_transform(y))
    y_lr_pred = lr.predict(X)
    assert_allclose(y_pred, transformer2.inverse_transform(y_lr_pred))
    assert_allclose(regr.regressor_.coef_, lr.coef_)


@pytest.mark.parametrize(
    "X,y", [friedman, (friedman[0], np.vstack((friedman[1], friedman[1] ** 2 + 1)).T)]
)
def test_transform_target_regressor_2d_transformer(X, y):
    # Check consistency with transformer accepting only 2D array and a 1D/2D y
    # array.
    transformer = StandardScaler()
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(), transformer=transformer
    )
    y_pred = regr.fit(X, y).predict(X)
    assert y.shape == y_pred.shape
    # consistency forward transform
    if y.ndim == 1:  # create a 2D array and squeeze results
        y_tran = regr.transformer_.transform(y.reshape(-1, 1))
    else:
        y_tran = regr.transformer_.transform(y)
    _check_standard_scaled(y, y_tran.squeeze())
    assert y.shape == y_pred.shape
    # consistency inverse transform
    assert_allclose(y, regr.transformer_.inverse_transform(y_tran).squeeze())
    # consistency of the regressor
    lr = LinearRegression()
    transformer2 = clone(transformer)
    if y.ndim == 1:  # create a 2D array and squeeze results
        lr.fit(X, transformer2.fit_transform(y.reshape(-1, 1)).squeeze())
        y_lr_pred = lr.predict(X).reshape(-1, 1)
        y_pred2 = transformer2.inverse_transform(y_lr_pred).squeeze()
    else:
        lr.fit(X, transformer2.fit_transform(y))
        y_lr_pred = lr.predict(X)
        y_pred2 = transformer2.inverse_transform(y_lr_pred)

    assert_allclose(y_pred, y_pred2)
    assert_allclose(regr.regressor_.coef_, lr.coef_)


def test_transform_target_regressor_2d_transformer_multioutput():
    # Check consistency with transformer accepting only 2D array and a 2D y
    # array.
    X = friedman[0]
    y = np.vstack((friedman[1], friedman[1] ** 2 + 1)).T
    transformer = StandardScaler()
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(), transformer=transformer
    )
    y_pred = regr.fit(X, y).predict(X)
    assert y.shape == y_pred.shape
    # consistency forward transform
    y_tran = regr.transformer_.transform(y)
    _check_standard_scaled(y, y_tran)
    assert y.shape == y_pred.shape
    # consistency inverse transform
    assert_allclose(y, regr.transformer_.inverse_transform(y_tran).squeeze())
    # consistency of the regressor
    lr = LinearRegression()
    transformer2 = clone(transformer)
    lr.fit(X, transformer2.fit_transform(y))
    y_lr_pred = lr.predict(X)
    assert_allclose(y_pred, transformer2.inverse_transform(y_lr_pred))
    assert_allclose(regr.regressor_.coef_, lr.coef_)


def test_transform_target_regressor_3d_target():
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/18866
    # Check with a 3D target with a transformer that reshapes the target
    X = friedman[0]
    y = np.tile(friedman[1].reshape(-1, 1, 1), [1, 3, 2])

    def flatten_data(data):
        return data.reshape(data.shape[0], -1)

    def unflatten_data(data):
        return data.reshape(data.shape[0], -1, 2)

    transformer = FunctionTransformer(func=flatten_data, inverse_func=unflatten_data)
    regr = TransformedTargetRegressor(
        regressor=LinearRegression(), transformer=transformer
    )
    y_pred = regr.fit(X, y).predict(X)
    assert y.shape == y_pred.shape


def test_transform_target_regressor_multi_to_single():
    X = friedman[0]
    y = np.transpose([friedman[1], (friedman[1] ** 2 + 1)])

    def func(y):
        out = np.sqrt(y[:, 0] ** 2 + y[:, 1] ** 2)
        return out[:, np.newaxis]

    def inverse_func(y):
        return y

    tt = TransformedTargetRegressor(
        func=func, inverse_func=inverse_func, check_inverse=False
    )
    tt.fit(X, y)
    y_pred_2d_func = tt.predict(X)
    assert y_pred_2d_func.shape == (100, 1)

    # force that the function only return a 1D array
    def func(y):
        return np.sqrt(y[:, 0] ** 2 + y[:, 1] ** 2)

    tt = TransformedTargetRegressor(
        func=func, inverse_func=inverse_func, check_inverse=False
    )
    tt.fit(X, y)
    y_pred_1d_func = tt.predict(X)
    assert y_pred_1d_func.shape == (100, 1)

    assert_allclose(y_pred_1d_func, y_pred_2d_func)


class DummyCheckerArrayTransformer(TransformerMixin, BaseEstimator):
    def fit(self, X, y=None):
        assert isinstance(X, np.ndarray)
        return self

    def transform(self, X):
        assert isinstance(X, np.ndarray)
        return X

    def inverse_transform(self, X):
        assert isinstance(X, np.ndarray)
        return X


class DummyCheckerListRegressor(DummyRegressor):
    def fit(self, X, y, sample_weight=None):
        assert isinstance(X, list)
        return super().fit(X, y, sample_weight)

    def predict(self, X):
        assert isinstance(X, list)
        return super().predict(X)


def test_transform_target_regressor_ensure_y_array():
    # check that the target ``y`` passed to the transformer will always be a
    # numpy array. Similarly, if ``X`` is passed as a list, we check that the
    # predictor receive as it is.
    X, y = friedman
    tt = TransformedTargetRegressor(
        transformer=DummyCheckerArrayTransformer(),
        regressor=DummyCheckerListRegressor(),
        check_inverse=False,
    )
    tt.fit(X.tolist(), y.tolist())
    tt.predict(X.tolist())
    with pytest.raises(AssertionError):
        tt.fit(X, y.tolist())
    with pytest.raises(AssertionError):
        tt.predict(X)


class DummyTransformer(TransformerMixin, BaseEstimator):
    """Dummy transformer which count how many time fit was called."""

    def __init__(self, fit_counter=0):
        self.fit_counter = fit_counter

    def fit(self, X, y=None):
        self.fit_counter += 1
        return self

    def transform(self, X):
        return X

    def inverse_transform(self, X):
        return X


@pytest.mark.parametrize("check_inverse", [False, True])
def test_transform_target_regressor_count_fit(check_inverse):
    # regression test for gh-issue #11618
    # check that we only call a single time fit for the transformer
    X, y = friedman
    ttr = TransformedTargetRegressor(
        transformer=DummyTransformer(), check_inverse=check_inverse
    )
    ttr.fit(X, y)
    assert ttr.transformer_.fit_counter == 1


class DummyRegressorWithExtraFitParams(DummyRegressor):
    def fit(self, X, y, sample_weight=None, check_input=True):
        # on the test below we force this to false, we make sure this is
        # actually passed to the regressor
        assert not check_input
        return super().fit(X, y, sample_weight)


def test_transform_target_regressor_pass_fit_parameters():
    X, y = friedman
    regr = TransformedTargetRegressor(
        regressor=DummyRegressorWithExtraFitParams(), transformer=DummyTransformer()
    )

    regr.fit(X, y, check_input=False)
    assert regr.transformer_.fit_counter == 1


def test_transform_target_regressor_route_pipeline():
    X, y = friedman

    regr = TransformedTargetRegressor(
        regressor=DummyRegressorWithExtraFitParams(), transformer=DummyTransformer()
    )
    estimators = [("normalize", StandardScaler()), ("est", regr)]

    pip = Pipeline(estimators)
    pip.fit(X, y, **{"est__check_input": False})

    assert regr.transformer_.fit_counter == 1


class DummyRegressorWithExtraPredictParams(DummyRegressor):
    def predict(self, X, check_input=True):
        # In the test below we make sure that the check input parameter is
        # passed as false
        self.predict_called = True
        assert not check_input
        return super().predict(X)


def test_transform_target_regressor_pass_extra_predict_parameters():
    # Checks that predict kwargs are passed to regressor.
    X, y = friedman
    regr = TransformedTargetRegressor(
        regressor=DummyRegressorWithExtraPredictParams(), transformer=DummyTransformer()
    )

    regr.fit(X, y)
    regr.predict(X, check_input=False)
    assert regr.regressor_.predict_called


@pytest.mark.parametrize("output_format", ["pandas", "polars"])
def test_transform_target_regressor_not_warns_with_global_output_set(output_format):
    """Test that TransformedTargetRegressor will not raise warnings if
    set_config(transform_output="pandas"/"polars") is set globally; regression test for
    issue #29361."""
    X, y = datasets.make_regression()
    y = np.abs(y) + 1
    with config_context(transform_output=output_format):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            TransformedTargetRegressor(
                regressor=LinearRegression(), func=np.log, inverse_func=np.exp
            ).fit(X, y)


class ValidateDimensionRegressor(BaseEstimator):
    """A regressor that expects the target to have a specific number of dimensions."""

    def __init__(self, ndim):
        self.ndim = ndim

    def fit(self, X, y):
        assert y.ndim == self.ndim

    def predict(self, X):
        pass  # pragma: no cover


@pytest.mark.parametrize("ndim", [1, 2])
def test_transform_target_regressor_preserves_input_shape(ndim):
    """Check that TransformedTargetRegressor internally preserves the shape of the input

    non-regression test for issue #26530.
    """
    X, y = datasets.make_regression(n_samples=10, n_features=5, random_state=42)
    if ndim == 2:
        y = y.reshape(-1, 1)

    regr = TransformedTargetRegressor(regressor=ValidateDimensionRegressor(ndim))
    regr.fit(X, y)


@config_context(enable_metadata_routing=True)
def test_transform_target_regressor_metadata_routing_default_estimator():
    """Test that metadata request is set on the default regressor"""
    X, y = make_regression()
    ttr = TransformedTargetRegressor()
    ttr.fit(X, y, sample_weight=np.ones(shape=(X.shape[0],)))


# ---------------------------------------------------------------------------
# Bias correction tests
# ---------------------------------------------------------------------------


def _make_log_regression_data(n_samples=500, noise_scale=0.5, random_state=42):
    """Generate data where the true relationship is Y = exp(linear + noise).

    A log transform on Y linearises the problem, but inverse-transforming
    the predicted mean in log-space with exp() introduces a systematic
    downward bias (Jensen's inequality).
    """
    rng = np.random.RandomState(random_state)
    X = rng.randn(n_samples, 3)
    z_true = 2 * X[:, 0] + X[:, 1] - 0.5 * X[:, 2]
    z_noisy = z_true + noise_scale * rng.randn(n_samples)
    y = np.exp(z_noisy)
    return X, y


def test_bias_correction_multiplicative_reduces_bias():
    X, y = _make_log_regression_data()

    tt_none = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction=None,
    )
    tt_mult = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction="multiplicative",
    )

    tt_none.fit(X, y)
    tt_mult.fit(X, y)

    pred_none = tt_none.predict(X)
    pred_mult = tt_mult.predict(X)

    bias_none = np.abs(np.mean(pred_none) - np.mean(y))
    bias_mult = np.abs(np.mean(pred_mult) - np.mean(y))

    assert bias_mult < bias_none, (
        f"Multiplicative correction did not reduce bias: "
        f"uncorrected={bias_none:.4f}, corrected={bias_mult:.4f}"
    )


def test_bias_correction_additive_reduces_bias():
    X, y = _make_log_regression_data()

    tt_none = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction=None,
    )
    tt_add = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction="additive",
    )

    tt_none.fit(X, y)
    tt_add.fit(X, y)

    pred_none = tt_none.predict(X)
    pred_add = tt_add.predict(X)

    bias_none = np.abs(np.mean(pred_none) - np.mean(y))
    bias_add = np.abs(np.mean(pred_add) - np.mean(y))

    assert bias_add < bias_none, (
        f"Additive correction did not reduce bias: "
        f"uncorrected={bias_none:.4f}, corrected={bias_add:.4f}"
    )


def test_bias_correction_taylor_reduces_bias():
    X, y = _make_log_regression_data()

    tt_none = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction=None,
    )
    tt_taylor = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction="taylor",
    )

    tt_none.fit(X, y)
    tt_taylor.fit(X, y)

    pred_none = tt_none.predict(X)
    pred_taylor = tt_taylor.predict(X)

    bias_none = np.abs(np.mean(pred_none) - np.mean(y))
    bias_taylor = np.abs(np.mean(pred_taylor) - np.mean(y))

    assert bias_taylor < bias_none, (
        f"Taylor correction did not reduce bias: "
        f"uncorrected={bias_none:.4f}, corrected={bias_taylor:.4f}"
    )


def test_bias_correction_taylor_per_sample():
    """Taylor correction should produce different adjustments for predictions
    at different magnitudes, unlike the global modes."""
    X, y = _make_log_regression_data()

    tt = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction="taylor",
    )
    tt.fit(X, y)

    tt_none = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction=None,
    )
    tt_none.fit(X, y)

    pred_taylor = tt.predict(X)
    pred_none = tt_none.predict(X)

    correction = pred_taylor - pred_none
    assert np.std(correction) > 1e-10, (
        "Taylor correction should vary across samples"
    )


def test_bias_correction_none_default():
    """bias_correction=None should behave identically to the original code."""
    X, y = _make_log_regression_data()

    tt_default = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
    )
    tt_explicit = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction=None,
    )

    tt_default.fit(X, y)
    tt_explicit.fit(X, y)

    assert_allclose(tt_default.predict(X), tt_explicit.predict(X))
    assert tt_default.bias_correction_factor_ is None
    assert tt_default.residual_variance_ is None


def test_bias_correction_identity_transform():
    """With an identity transform there is no Jensen bias, so correction
    should be effectively a no-op (factor ~1 or ~0, Hessian ~0)."""
    X, y = _make_log_regression_data()

    for mode in ("additive", "multiplicative", "taylor"):
        tt = TransformedTargetRegressor(
            regressor=LinearRegression(),
            bias_correction=mode,
        )
        tt_none = TransformedTargetRegressor(
            regressor=LinearRegression(),
            bias_correction=None,
        )

        tt.fit(X, y)
        tt_none.fit(X, y)

        assert_allclose(
            tt.predict(X),
            tt_none.predict(X),
            rtol=1e-5,
            err_msg=f"Identity transform with mode={mode} should be a no-op",
        )


def test_bias_correction_multioutput():
    """Per-column factors and variance for 2D targets."""
    rng = np.random.RandomState(0)
    X = rng.randn(200, 3)
    z = X @ np.array([[1, 2], [0.5, -1], [0.3, 0.7]])
    noise = 0.3 * rng.randn(200, 2)
    y = np.exp(z + noise)

    for mode in ("additive", "multiplicative", "taylor"):
        tt = TransformedTargetRegressor(
            regressor=LinearRegression(),
            func=np.log,
            inverse_func=np.exp,
            bias_correction=mode,
        )
        tt.fit(X, y)
        pred = tt.predict(X)
        assert pred.shape == y.shape, f"Shape mismatch for mode={mode}"

        if mode in ("additive", "multiplicative"):
            factor = tt.bias_correction_factor_
            assert hasattr(factor, "__len__") and len(factor) == 2, (
                f"Expected per-column factor for mode={mode}, got {factor}"
            )
        else:
            var = tt.residual_variance_
            assert hasattr(var, "__len__") and len(var) == 2, (
                f"Expected per-column variance for mode={mode}, got {var}"
            )


def test_bias_correction_multiplicative_near_zero_warns():
    """When inverse-transformed training predictions have near-zero mean,
    multiplicative mode should warn and default the factor to 1.0."""
    rng = np.random.RandomState(42)
    X = rng.randn(100, 2)
    y = rng.randn(100)

    def to_centered(x):
        return x

    def from_centered(x):
        return x * 1e-12

    tt = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=to_centered,
        inverse_func=from_centered,
        check_inverse=False,
        bias_correction="multiplicative",
    )
    with pytest.warns(
        UserWarning,
        match="near zero",
    ):
        tt.fit(X, y)
    assert_allclose(tt.bias_correction_factor_, 1.0, atol=1e-5)


def test_bias_correction_taylor_nonsmooth_warns():
    """When the Hessian is non-finite (e.g., from a piecewise transform),
    Taylor mode should warn and skip correction for affected samples."""
    rng = np.random.RandomState(42)
    X = rng.randn(100, 2)
    y = np.abs(rng.randn(100)) + 1

    def piecewise_func(x):
        return np.where(x > 1.5, np.log(x), x - 1.5 + np.log(1.5))

    def piecewise_inv(x):
        threshold = np.log(1.5)
        return np.where(x > threshold, np.exp(x), x + 1.5 - np.log(1.5))

    tt = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=piecewise_func,
        inverse_func=piecewise_inv,
        check_inverse=False,
        bias_correction="taylor",
    )
    tt.fit(X, y)
    pred = tt.predict(X)
    assert pred.shape == (100,)
    assert np.all(np.isfinite(pred))


def test_bias_correction_with_power_transformer():
    """All modes should work with a transformer object, not just func/inverse_func."""
    from sklearn.preprocessing import PowerTransformer

    rng = np.random.RandomState(42)
    X = rng.randn(200, 3)
    y = np.abs(X @ np.array([1.0, 0.5, 0.3]) + 0.2 * rng.randn(200)) + 1.0

    for mode in ("additive", "multiplicative", "taylor"):
        tt = TransformedTargetRegressor(
            regressor=LinearRegression(),
            transformer=PowerTransformer(method="yeo-johnson"),
            bias_correction=mode,
        )
        tt.fit(X, y)
        pred = tt.predict(X)
        assert pred.shape == y.shape, f"Shape mismatch for mode={mode}"
        assert np.all(np.isfinite(pred)), f"Non-finite preds for mode={mode}"


def test_bias_correction_clone_preserves_param():
    """clone() should keep the bias_correction param but drop fitted attrs."""
    tt = TransformedTargetRegressor(
        regressor=LinearRegression(),
        func=np.log,
        inverse_func=np.exp,
        bias_correction="multiplicative",
    )
    X, y = _make_log_regression_data()
    tt.fit(X, y)

    assert tt.bias_correction_factor_ is not None

    tt_cloned = clone(tt)
    assert tt_cloned.bias_correction == "multiplicative"
    assert not hasattr(tt_cloned, "bias_correction_factor_")
    assert not hasattr(tt_cloned, "residual_variance_")


def test_bias_correction_factor_attribute_shape():
    """bias_correction_factor_ should be a scalar for 1D y and an array
    matching n_outputs for 2D y."""
    X, y_1d = _make_log_regression_data(n_samples=200)

    for mode in ("additive", "multiplicative"):
        tt_1d = TransformedTargetRegressor(
            regressor=LinearRegression(),
            func=np.log,
            inverse_func=np.exp,
            bias_correction=mode,
        )
        tt_1d.fit(X, y_1d)
        assert np.isscalar(tt_1d.bias_correction_factor_), (
            f"Expected scalar factor for 1D y with mode={mode}"
        )

    rng = np.random.RandomState(0)
    y_2d = np.column_stack([y_1d, np.abs(rng.randn(200)) + 1])

    for mode in ("additive", "multiplicative"):
        tt_2d = TransformedTargetRegressor(
            regressor=LinearRegression(),
            func=np.log,
            inverse_func=np.exp,
            bias_correction=mode,
        )
        tt_2d.fit(X, y_2d)
        factor = tt_2d.bias_correction_factor_
        assert hasattr(factor, "__len__") and len(factor) == 2, (
            f"Expected 2-element factor for 2D y with mode={mode}, got {factor}"
        )
