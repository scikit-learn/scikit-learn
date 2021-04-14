import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest
from pytest import approx
from scipy.optimize import (
    minimize,
    minimize_scalar,
    newton,
)
from scipy.special import logit

from sklearn._loss.link import _inclusive_low_high
from sklearn._loss.loss import (
    _LOSSES,
    AbsoluteError,
    BinaryCrossEntropy,
    CategoricalCrossEntropy,
    HalfGammaLoss,
    HalfPoissonLoss,
    HalfSquaredError,
    HalfTweedieLoss,
    PinballLoss,
)
from sklearn.utils import assert_all_finite
from sklearn.utils._testing import skip_if_32bit
from sklearn.utils.fixes import sp_version, parse_version


ALL_LOSSES = list(_LOSSES.values())

LOSS_INSTANCES = [loss() for loss in ALL_LOSSES]
# HalfTweedieLoss(power=1.5) is already there as default
LOSS_INSTANCES += [
    PinballLoss(quantile=0.25),
    HalfTweedieLoss(power=-1.5),
    HalfTweedieLoss(power=0),
    HalfTweedieLoss(power=1),
    HalfTweedieLoss(power=2),
    HalfTweedieLoss(power=3.0),
]


def loss_instance_name(loss):
    name = loss.__class__.__name__
    if hasattr(loss, "quantile"):
        name += f"(quantile={loss.quantile})"
    elif hasattr(loss, "power"):
        name += f"(power={loss.power})"
    return name


def random_y_true_raw_prediction(
    loss, n_samples, y_bound=(-100, 100), raw_bound=(-5, 5), seed=42
):
    """Random generate y_true and raw_prediction in valid range."""
    rng = np.random.RandomState(seed)
    if loss.n_classes <= 2:
        raw_prediction = rng.uniform(
            low=raw_bound[0], high=raw_bound[0], size=n_samples
        )
        # generate a y_true in valid range
        low, high = _inclusive_low_high(loss.interval_y_true)
        low = max(low, y_bound[0])
        high = min(high, y_bound[1])
        y_true = rng.uniform(low, high, size=n_samples)
        # set some values at special boundaries
        if (
            loss.interval_y_true.low == 0
            and loss.interval_y_true.low_inclusive
        ):
            y_true[:: (n_samples // 3)] = 0
        if (
            loss.interval_y_true.high == 1
            and loss.interval_y_true.high_inclusive
        ):
            y_true[1:: (n_samples // 3)] = 1
    else:
        raw_prediction = np.empty((n_samples, loss.n_classes))
        raw_prediction.flat[:] = rng.uniform(
            low=raw_bound[0],
            high=raw_bound[1],
            size=n_samples * loss.n_classes,
        )
        y_true = np.arange(n_samples).astype(float) % loss.n_classes

    return y_true, raw_prediction


def numerical_derivative(func, x, eps):
    """Helper function for numerical (first) derivatives."""
    # For numerical derivatives, see
    # https://en.wikipedia.org/wiki/Numerical_differentiation
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient
    # We use central finite differences of accuracy 4.
    h = np.full_like(x, fill_value=eps)
    f_minus_2h = func(x - 2 * h)
    f_minus_1h = func(x - h)
    f_plus_1h = func(x + h)
    f_plus_2h = func(x + 2 * h)
    return (-f_plus_2h + 8 * f_plus_1h - 8 * f_minus_1h + f_minus_2h) / (
        12.0 * eps
    )


@pytest.mark.parametrize("loss", LOSS_INSTANCES, ids=loss_instance_name)
def test_loss_boundary(loss):
    # make sure low and high are always within the interval, used for linspace
    if loss.n_classes is None or loss.n_classes <= 2:
        low, high = _inclusive_low_high(loss.interval_y_true)
        y_true = np.linspace(low, high, num=10)
    else:
        y_true = np.linspace(0, 9, num=10)

    # add boundaries if they are included
    if loss.interval_y_true.low_inclusive:
        y_true = np.r_[y_true, loss.interval_y_true.low]
    if loss.interval_y_true.high_inclusive:
        y_true = np.r_[y_true, loss.interval_y_true.high]

    assert loss.in_y_true_range(y_true)

    low, high = _inclusive_low_high(loss.interval_y_pred)
    if loss.n_classes is None or loss.n_classes <= 2:
        y_pred = np.linspace(low, high, num=10)
    else:
        y_pred = np.empty((10, 3))
        y_pred[:, 0] = np.linspace(low, high, num=10)
        y_pred[:, 1] = 0.5 * (1 - y_pred[:, 0])
        y_pred[:, 2] = 0.5 * (1 - y_pred[:, 0])

    assert loss.in_y_pred_range(y_pred)

    # calculating losses should not fail
    raw_prediction = loss.link(y_pred)
    loss.loss(y_true=y_true, raw_prediction=raw_prediction)


@pytest.mark.parametrize(
    "loss, y_true_success, y_true_fail",
    [
        (HalfSquaredError(), [-100, 0, 0.1, 100], [-np.inf, np.inf]),
        (AbsoluteError(), [-100, 0, 0.1, 100], [-np.inf, np.inf]),
        (PinballLoss(), [-100, 0, 0.1, 100], [-np.inf, np.inf]),
        (HalfPoissonLoss(), [0, 0.1, 100], [-np.inf, -3, -0.1, np.inf]),
        (HalfGammaLoss(), [0.1, 100], [-np.inf, -3, -0.1, 0, np.inf]),
        (HalfTweedieLoss(power=-3), [-100, 0, 0.1, 100], [-np.inf, np.inf]),
        (HalfTweedieLoss(power=0), [-100, 0, 0.1, 100], [-np.inf, np.inf]),
        (
            HalfTweedieLoss(power=1.5),
            [0, 0.1, 100],
            [-np.inf, -3, -0.1, np.inf],
        ),
        (HalfTweedieLoss(power=2), [0.1, 100], [-np.inf, -3, -0.1, 0, np.inf]),
        (HalfTweedieLoss(power=3), [0.1, 100], [-np.inf, -3, -0.1, 0, np.inf]),
        (BinaryCrossEntropy(), [0, 0.5, 1], [-np.inf, -1, 2, np.inf]),
        (CategoricalCrossEntropy(), [0.0, 1.0, 2], [-np.inf, -1, 1.1, np.inf]),
    ],
)
def test_loss_boundary_y_true(loss, y_true_success, y_true_fail):
    # Test boundaries of y_true for loss functions.
    for y in y_true_success:
        assert loss.in_y_true_range(np.array([y]))
    for y in y_true_fail:
        assert not loss.in_y_true_range(np.array([y]))


@pytest.mark.parametrize(
    "loss, y_pred_success, y_pred_fail",
    [
        (HalfSquaredError(), [-100, 0, 0.1, 100], [-np.inf, np.inf]),
        (AbsoluteError(), [-100, 0, 0.1, 100], [-np.inf, np.inf]),
        (PinballLoss(), [-100, 0, 0.1, 100], [-np.inf, np.inf]),
        (HalfPoissonLoss(), [0.1, 100], [-np.inf, -3, -0.1, 0, np.inf]),
        (HalfGammaLoss(), [0.1, 100], [-np.inf, -3, -0.1, 0, np.inf]),
        (
            HalfTweedieLoss(power=-3),
            [0.1, 100],
            [-np.inf, -3, -0.1, 0, np.inf],
        ),
        (HalfTweedieLoss(power=0), [0.1, 100], [-np.inf, -3, -0.1, 0, np.inf]),
        (
            HalfTweedieLoss(power=1.5),
            [0.1, 100],
            [-np.inf, -3, -0.1, 0, np.inf],
        ),
        (HalfTweedieLoss(power=2), [0.1, 100], [-np.inf, -3, -0.1, 0, np.inf]),
        (HalfTweedieLoss(power=3), [0.1, 100], [-np.inf, -3, -0.1, 0, np.inf]),
        (BinaryCrossEntropy(), [0.1, 0.5], [-np.inf, 0, 1, np.inf]),
        (CategoricalCrossEntropy(), [0.1, 0.5], [-np.inf, 0, 1, np.inf]),
    ],
)
def test_loss_boundary_y_pred(loss, y_pred_success, y_pred_fail):
    # Test boundaries of y_pred for loss functions.
    for y in y_pred_success:
        assert loss.in_y_pred_range(np.array([y]))
    for y in y_pred_fail:
        assert not loss.in_y_pred_range(np.array([y]))


@pytest.mark.parametrize("loss", ALL_LOSSES)
@pytest.mark.parametrize("dtype_in", [np.float32, np.float64])
@pytest.mark.parametrize("dtype_out", [np.float32, np.float64])
@pytest.mark.parametrize("sample_weight", [None, 1])
@pytest.mark.parametrize("out1", [None, 1])
@pytest.mark.parametrize("out2", [None, 1])
@pytest.mark.parametrize("n_threads", [1, 2])
def test_loss_dtype(
    loss, dtype_in, dtype_out, sample_weight, out1, out2, n_threads
):
    # Test that loss accepts if all input arrays are either all float32 or all
    # float64, and all output arrays are either all float32 or all float64.
    loss = loss()
    if loss.n_classes <= 2:
        # generate a y_true in valid range
        low, high = _inclusive_low_high(loss.interval_y_true, dtype=dtype_in)
        y_true = np.array([0.5 * (high - low)], dtype=dtype_in)
        raw_prediction = np.array([0.0], dtype=dtype_in)
    else:
        y_true = np.array([0], dtype=dtype_in)
        raw_prediction = np.full(
            shape=(1, loss.n_classes), fill_value=0.0, dtype=dtype_in
        )

    if sample_weight is not None:
        sample_weight = np.array([2.0], dtype=dtype_in)
    if out1 is not None:
        out1 = np.empty_like(y_true, dtype=dtype_out)
    if out2 is not None:
        out2 = np.empty_like(raw_prediction, dtype=dtype_out)

    loss.loss(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        loss=out1,
        n_threads=n_threads,
    )
    loss.gradient(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        gradient=out2,
        n_threads=n_threads,
    )
    loss.loss_gradient(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        loss=out1,
        gradient=out2,
        n_threads=n_threads,
    )
    if out1 is not None and loss.n_classes >= 3:
        out1 = np.empty_like(raw_prediction, dtype=dtype_out)
    loss.gradient_hessian(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        gradient=out1,
        hessian=out2,
        n_threads=n_threads,
    )


@pytest.mark.parametrize("loss", LOSS_INSTANCES, ids=loss_instance_name)
@pytest.mark.parametrize("sample_weight", [None, "range"])
def test_loss_same_as_C_functions(loss, sample_weight):
    y_true, raw_prediction = random_y_true_raw_prediction(
        loss=loss,
        n_samples=20,
        y_bound=(-100, 100),
        raw_bound=(-10, 10),
        seed=42,
    )
    if sample_weight == "range":
        sample_weight = np.linspace(1, y_true.shape[0], num=y_true.shape[0])

    out_l1 = np.empty_like(y_true)
    out_l2 = np.empty_like(y_true)
    out_g1 = np.empty_like(raw_prediction)
    out_g2 = np.empty_like(raw_prediction)
    out_h1 = np.empty_like(raw_prediction)
    out_h2 = np.empty_like(raw_prediction)
    assert_allclose(
        loss.loss(
            y_true=y_true,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            loss=out_l1,
        ),
        loss._loss(
            y_true=y_true,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            loss=out_l2,
        ),
    )
    assert_allclose(
        loss.gradient(
            y_true=y_true,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            gradient=out_g1,
        ),
        loss._gradient(
            y_true=y_true,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            gradient=out_g2,
        ),
    )
    loss.loss_gradient(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        loss=out_l1,
        gradient=out_g1,
    )
    loss._loss_gradient(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        loss=out_l2,
        gradient=out_g2,
    )
    assert_allclose(out_l1, out_l2)
    assert_allclose(out_g1, out_g2)
    loss.gradient_hessian(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        gradient=out_g1,
        hessian=out_h1,
    )
    loss._gradient_hessian(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        gradient=out_g2,
        hessian=out_h2,
    )
    assert_allclose(out_g1, out_g2)
    assert_allclose(out_h1, out_h2)


@pytest.mark.parametrize("loss", LOSS_INSTANCES, ids=loss_instance_name)
@pytest.mark.parametrize("sample_weight", [None, "range"])
def test_loss_gradients_are_the_same(loss, sample_weight):
    # Test that loss and gradient are the same across different functions.
    # Also test that output arguments contain correct result.
    y_true, raw_prediction = random_y_true_raw_prediction(
        loss=loss,
        n_samples=20,
        y_bound=(-100, 100),
        raw_bound=(-10, 10),
        seed=42,
    )
    if sample_weight == "range":
        sample_weight = np.linspace(1, y_true.shape[0], num=y_true.shape[0])

    out_l1 = np.empty_like(y_true)
    out_l2 = np.empty_like(y_true)
    out_g1 = np.empty_like(raw_prediction)
    out_g2 = np.empty_like(raw_prediction)
    out_g3 = np.empty_like(raw_prediction)
    out_h3 = np.empty_like(raw_prediction)

    l1 = loss.loss(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        loss=out_l1,
    )
    g1 = loss.gradient(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        gradient=out_g1,
    )
    l2, g2 = loss.loss_gradient(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        loss=out_l2,
        gradient=out_g2,
    )
    g3, h3 = loss.gradient_hessian(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
        gradient=out_g3,
        hessian=out_h3,
    )
    assert_allclose(l1, l2)
    assert_array_equal(l1, out_l1)
    assert np.shares_memory(l1, out_l1)
    assert_array_equal(l2, out_l2)
    assert np.shares_memory(l2, out_l2)
    assert_allclose(g1, g2)
    assert_allclose(g1, g3)
    assert_array_equal(g1, out_g1)
    assert np.shares_memory(g1, out_g1)
    assert_array_equal(g2, out_g2)
    assert np.shares_memory(g2, out_g2)
    assert_array_equal(g3, out_g3)
    assert np.shares_memory(g3, out_g3)

    if hasattr(loss, "gradient_proba"):
        assert loss.n_classes >= 3  # only for CategoricalCrossEntropy
        out_g4 = np.empty_like(raw_prediction)
        out_proba = np.empty_like(raw_prediction)
        g4, proba = loss.gradient_proba(
            y_true=y_true,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
            gradient=out_g4,
            proba=out_proba,
        )
        assert_allclose(g1, out_g4)
        assert_allclose(g1, g4)
        assert_allclose(proba, out_proba)
        assert_allclose(np.sum(proba, axis=1), 1)


@pytest.mark.parametrize("loss", LOSS_INSTANCES, ids=loss_instance_name)
@pytest.mark.parametrize("sample_weight", ["ones", "random"])
def test_sample_weight_multiplies_gradients(loss, sample_weight):
    # Make sure that passing sample weights to the gradient and hessians
    # computation methods is equivalent to multiplying by the weights.

    n_samples = 100
    y_true, raw_prediction = random_y_true_raw_prediction(
        loss=loss,
        n_samples=n_samples,
        y_bound=(-100, 100),
        raw_bound=(-5, 5),
        seed=42,
    )

    if sample_weight == "ones":
        sample_weight = np.ones(shape=n_samples, dtype=np.float64)
    else:
        rng = np.random.RandomState(42)
        sample_weight = rng.normal(size=n_samples).astype(np.float64)

    baseline_prediction = loss.fit_intercept_only(
        y_true=y_true, sample_weight=None
    )

    if loss.n_classes <= 2:
        raw_prediction = np.zeros(
            shape=(n_samples,), dtype=baseline_prediction.dtype
        )
    else:
        raw_prediction = np.zeros(
            shape=(n_samples, loss.n_classes), dtype=baseline_prediction.dtype
        )
    raw_prediction += baseline_prediction

    gradient, hessian = loss.gradient_hessian(
        y_true=y_true, raw_prediction=raw_prediction, sample_weight=None
    )

    gradient_sw, hessian_sw = loss.gradient_hessian(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
    )

    if loss.n_classes <= 2:
        assert_allclose(gradient * sample_weight, gradient_sw)
        assert_allclose(hessian * sample_weight, hessian_sw)
    else:
        assert_allclose(gradient * sample_weight[:, None], gradient_sw)
        assert_allclose(hessian * sample_weight[:, None], hessian_sw)


@pytest.mark.parametrize("loss", LOSS_INSTANCES, ids=loss_instance_name)
@pytest.mark.parametrize("sample_weight", [None, "range"])
def test_loss_of_perfect_prediction(loss, sample_weight):
    # Test that loss of y_true = y_pred plus constant_to_optimal_zero sums up
    # to zero.
    if loss.n_classes <= 2:
        # Use small values such that exp(value) is not nan.
        raw_prediction = np.array([-10, -0.1, 0, 0.1, 3, 10])
        y_true = loss.inverse(raw_prediction)
    else:
        # CategoricalCrossEntropy
        y_true = np.arange(loss.n_classes).astype(float)
        # raw_prediction with entries -exp(10), but +exp(10) on the diagonal
        # this is close enough to np.inf which would produce nan
        raw_prediction = np.full(
            shape=(loss.n_classes, loss.n_classes),
            fill_value=-np.exp(10),
            dtype=float,
        )
        raw_prediction.flat[:: loss.n_classes + 1] = np.exp(10)

    if sample_weight == "range":
        sample_weight = np.linspace(1, y_true.shape[0], num=y_true.shape[0])

    loss_value = loss.loss(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
    )
    constant_term = loss.constant_to_optimal_zero(
        y_true=y_true, sample_weight=sample_weight
    )
    # Comparing loss_value + constant_term to zero would result in large
    # round-off errors.
    assert_allclose(loss_value, -constant_term, atol=1e-14, rtol=1e-15)


@pytest.mark.parametrize("loss", LOSS_INSTANCES, ids=loss_instance_name)
@pytest.mark.parametrize("sample_weight", [None, "range"])
def test_gradients_hessians_numerically(loss, sample_weight):
    # Test that gradients are computed correctly by comparing to numerical
    # derivatives of loss functions.
    # Test that hessians are correct by numerical derivative of gradients.
    n_samples = 20
    y_true, raw_prediction = random_y_true_raw_prediction(
        loss=loss,
        n_samples=n_samples,
        y_bound=(-100, 100),
        raw_bound=(-5, 5),
        seed=42,
    )

    if sample_weight == "range":
        sample_weight = np.linspace(1, y_true.shape[0], num=y_true.shape[0])

    g, h = loss.gradient_hessian(
        y_true=y_true,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
    )

    assert g.shape == raw_prediction.shape
    assert h.shape == raw_prediction.shape

    if loss.n_classes <= 2:

        def loss_func(x):
            return loss.loss(
                y_true=y_true, raw_prediction=x, sample_weight=sample_weight,
            )

        g_numeric = numerical_derivative(loss_func, raw_prediction, eps=1e-6)
        assert_allclose(g, g_numeric, rtol=5e-6, atol=1e-10)

        def grad_func(x):
            return loss.gradient(
                y_true=y_true, raw_prediction=x, sample_weight=sample_weight,
            )

        h_numeric = numerical_derivative(grad_func, raw_prediction, eps=1e-6)
        if loss.approx_hessian:
            assert np.all(h >= h_numeric)
        else:
            assert_allclose(h, h_numeric, rtol=5e-6, atol=1e-10)
    else:
        # For multiclass loss, we should only change the predictions of the
        # class for which the derivative is taken for, e.g. offset[:, k] = eps
        # for class k.
        # As a softmax is computed, offsetting the whole array by a constant
        # would have no effect on the probabilities, and thus on the loss.
        for k in range(loss.n_classes):

            def loss_func(x):
                raw = raw_prediction.copy()
                raw[:, k] = x
                return loss.loss(
                    y_true=y_true,
                    raw_prediction=raw,
                    sample_weight=sample_weight,
                )

            g_numeric = numerical_derivative(
                loss_func, raw_prediction[:, k], eps=1e-5
            )
            assert_allclose(g[:, k], g_numeric, rtol=5e-6, atol=1e-10)

            def grad_func(x):
                raw = raw_prediction.copy()
                raw[:, k] = x
                return loss.gradient(
                    y_true=y_true,
                    raw_prediction=raw,
                    sample_weight=sample_weight,
                )[:, k]

            h_numeric = numerical_derivative(
                grad_func, raw_prediction[:, k], eps=1e-6
            )
            if loss.approx_hessian:
                assert np.all(h >= h_numeric)
            else:
                assert_allclose(h[:, k], h_numeric, rtol=5e-6, atol=1e-10)


@pytest.mark.parametrize(
    "loss, x0, y_true",
    [
        ("squared_error", -2.0, 42),
        ("squared_error", 117.0, 1.05),
        ("squared_error", 0.0, 0.0),
        # The argmin of binary_crossentropy for y_true=0 and y_true=1 is resp.
        # -inf and +inf due to logit, cf. "complete separation". Therefore, we
        # use 0 < y_true < 1.
        ("binary_crossentropy", 0.3, 0.1),
        ("binary_crossentropy", -12, 0.2),
        ("binary_crossentropy", 30, 0.9),
        ("poisson_loss", 12.0, 1.0),
        ("poisson_loss", 0.0, 2.0),
        ("poisson_loss", -22.0, 10.0),
    ],
)
@pytest.mark.skipif(
    sp_version == parse_version("1.2.0"),
    reason="bug in scipy 1.2.0, see scipy issue #9608",
)
@skip_if_32bit
def test_derivatives(loss, x0, y_true):
    # Check that gradients are zero when the loss is minimized on a single
    # value/sample using Halley's method with the first and second order
    # derivatives computed by the Loss instance.
    # Note that methods of Loss instances operate on arrays while the newton
    # root finder expects a scalar or a one-element array for this purpose.

    loss = _LOSSES[loss](sample_weight=None)
    y_true = np.array([y_true], dtype=np.float64)
    x0 = np.array([x0], dtype=np.float64)

    def func(x: np.ndarray) -> np.ndarray:
        # Add constant term such that loss has its minimum at zero, which is
        # required by the newton method.
        return loss.loss(
            y_true=y_true, raw_prediction=x
        ) + loss.constant_to_optimal_zero(y_true=y_true)

    def fprime(x: np.ndarray) -> np.ndarray:
        return loss.gradient(y_true=y_true, raw_prediction=x)

    def fprime2(x: np.ndarray) -> np.ndarray:
        return loss.gradient_hessian(y_true=y_true, raw_prediction=x)[1]

    optimum = newton(
        func,
        x0=x0,
        fprime=fprime,
        fprime2=fprime2,
        maxiter=100,
        tol=5e-8,
    )

    # Need to ravel arrays because assert_allclose requires matching dimensions
    y_true = y_true.ravel()
    optimum = optimum.ravel()
    assert_allclose(loss.inverse(optimum), y_true)
    assert_allclose(func(optimum), 0, atol=1e-14)
    assert_allclose(
        loss.gradient(y_true=y_true, raw_prediction=optimum), 0, atol=5e-7
    )


@pytest.mark.parametrize("loss", LOSS_INSTANCES, ids=loss_instance_name)
@pytest.mark.parametrize("sample_weight", [None, "range"])
def test_loss_intercept_only(loss, sample_weight):
    # Test that fit_intercept_only returns the argmin of the loss and that the
    # gradient is zero.
    n_samples = 50
    if loss.n_classes <= 2:
        y_true = loss.inverse(np.linspace(-4, 4, num=n_samples))
    else:
        y_true = np.arange(n_samples).astype(float) % loss.n_classes
        y_true[::5] = 0  # exceedance of class 0

    if sample_weight == "range":
        sample_weight = np.linspace(0.1, 2, num=n_samples)

    a = loss.fit_intercept_only(y_true=y_true, sample_weight=sample_weight)

    # find minimum by optimization
    def fun(x):
        if loss.n_classes <= 2:
            raw_prediction = np.full(shape=(n_samples), fill_value=x)
        else:
            raw_prediction = np.ascontiguousarray(
                np.broadcast_to(x, shape=(n_samples, loss.n_classes))
            )
        return loss(
            y_true=y_true,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
        )

    if loss.n_classes <= 2:
        opt = minimize_scalar(fun, tol=1e-7, options={"maxiter": 100})
        grad = loss.gradient(
            y_true=y_true,
            raw_prediction=np.full_like(y_true, a),
            sample_weight=sample_weight,
        )
        assert a.shape == tuple()  # scalar
        assert a.dtype == y_true.dtype
        assert_all_finite(a)
        a == approx(opt.x, rel=1e-7)
        grad.sum() == approx(0, abs=1e-12)
    else:
        # constraint corresponds to sum(raw_prediction) = 0
        # without the constraint, we would need to apply
        # loss.symmetrize_raw_prediction to opt.x before comparing
        # TODO: With scipy 1.1.0, one could use
        # LinearConstraint(np.ones((1, loss.n_classes)), 0, 0)
        opt = minimize(
            fun,
            np.empty((loss.n_classes)),
            tol=1e-13,
            options={"maxiter": 100},
            method="SLSQP",
            constraints={
                "type": "eq",
                "fun": lambda x: np.ones((1, loss.n_classes)) @ x
            },
        )
        grad = loss.gradient(
            y_true=y_true,
            raw_prediction=np.tile(a, (n_samples, 1)),
            sample_weight=sample_weight,
        )
        assert a.dtype == y_true.dtype
        assert_all_finite(a)
        assert_allclose(a, opt.x, rtol=5e-6, atol=1e-12)
        assert_allclose(grad.sum(axis=0), 0, atol=1e-12)


@pytest.mark.parametrize(
    "loss, func, link, low, high, random_dist",
    [
        (HalfSquaredError, np.mean, "identity", None, None, "normal"),
        (AbsoluteError, np.median, "identity", None, None, "normal"),
        (HalfPoissonLoss, np.mean, np.log, 0, None, "poisson"),
        (BinaryCrossEntropy, np.mean, logit, 0, 1, "binomial"),
    ],
)
def test_specific_fit_intercept_only(loss, func, link, low, high, random_dist):
    rng = np.random.RandomState(0)
    loss = loss()
    if random_dist == "binomial":
        y_train = rng.binomial(1, 0.5, size=100)
    else:
        y_train = getattr(rng, random_dist)(size=100)
    baseline_prediction = loss.fit_intercept_only(y_true=y_train)
    # Make sure baseline prediction is the expected one, i.e. func, e.g.
    # mean or median.
    assert_all_finite(baseline_prediction)
    if link == "identity":
        assert baseline_prediction == approx(func(y_train))
        assert_allclose(loss.inverse(baseline_prediction), baseline_prediction)
    else:
        assert baseline_prediction == approx(link(func(y_train)))

    # Test baseline at boundary
    if low is not None:
        y_train.fill(low)
        baseline_prediction = loss.fit_intercept_only(y_true=y_train)
        assert_all_finite(baseline_prediction)
    if high is not None:
        y_train.fill(high)
        baseline_prediction = loss.fit_intercept_only(y_true=y_train)
        assert_all_finite(baseline_prediction)


def test_categorical_crossentropy_fit_intercept_only():
    rng = np.random.RandomState(0)
    n_classes = 4
    loss = CategoricalCrossEntropy(n_classes=n_classes)
    # Same logic as test_single_fit_intercept_only. Here inverse link function
    # = softmax and link function = log - symmetry term
    y_train = rng.randint(0, n_classes + 1, size=100).astype(np.float64)
    baseline_prediction = loss.fit_intercept_only(y_true=y_train)
    assert baseline_prediction.shape == (n_classes,)
    p = np.zeros(n_classes, dtype=y_train.dtype)
    for k in range(n_classes):
        p[k] = (y_train == k).mean()
    assert_allclose(baseline_prediction, np.log(p) - np.mean(np.log(p)))
    assert_allclose(baseline_prediction[None, :], loss.link(p[None, :]))

    for y_train in (np.zeros(shape=10), np.ones(shape=10)):
        y_train = y_train.astype(np.float64)
        baseline_prediction = loss.fit_intercept_only(y_true=y_train)
        assert baseline_prediction.dtype == y_train.dtype
        assert_all_finite(baseline_prediction)


def test_binary_and_categorical_crossentropy():
    # Test that CategoricalCrossEntropy with n_classes = 2 is the same as
    # BinaryCrossEntropy
    rng = np.random.RandomState(0)
    n_samples = 20
    bce = BinaryCrossEntropy()
    cce = CategoricalCrossEntropy(n_classes=2)
    y_train = rng.randint(0, 2, size=n_samples).astype(np.float64)
    raw_prediction = rng.normal(size=n_samples)
    raw_cce = np.empty((n_samples, 2))
    raw_cce[:, 0] = -0.5 * raw_prediction
    raw_cce[:, 1] = 0.5 * raw_prediction
    assert_allclose(
        bce.loss(y_true=y_train, raw_prediction=raw_prediction),
        cce.loss(y_true=y_train, raw_prediction=raw_cce)
    )
