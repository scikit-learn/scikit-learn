import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn import config_context
from sklearn._loss.link import (
    _LINKS,
    HalfLogitLink,
    Interval,
    MultinomialLogit,
    _inclusive_low_high,
)
from sklearn.utils._array_api import (
    _convert_to_numpy,
    yield_namespace_device_dtype_combinations,
)
from sklearn.utils._testing import _array_api_for_tests

LINK_FUNCTIONS = list(_LINKS.values())


def test_interval_raises():
    """Test that interval with low > high raises ValueError."""
    with pytest.raises(
        ValueError, match="One must have low <= high; got low=1, high=0."
    ):
        Interval(1, 0, False, False)


@pytest.mark.parametrize(
    "interval",
    [
        Interval(0, 1, False, False),
        Interval(0, 1, False, True),
        Interval(0, 1, True, False),
        Interval(0, 1, True, True),
        Interval(-float("inf"), float("inf"), False, False),
        Interval(-float("inf"), float("inf"), False, True),
        Interval(-float("inf"), float("inf"), True, False),
        Interval(-float("inf"), float("inf"), True, True),
        Interval(-10, -1, False, False),
        Interval(-10, -1, False, True),
        Interval(-10, -1, True, False),
        Interval(-10, -1, True, True),
    ],
)
def test_is_in_range(interval):
    """Test that low and high are always within the interval used for linspace."""
    low, high = _inclusive_low_high(interval)
    x = np.linspace(low, high, num=10)

    assert interval.includes(x)

    # x contains lower bound
    assert interval.includes(np.r_[x, interval.low]) == interval.low_inclusive

    # x contains upper bound
    assert interval.includes(np.r_[x, interval.high]) == interval.high_inclusive

    # x contains upper and lower bound
    assert interval.includes(np.r_[x, interval.low, interval.high]) == (
        interval.low_inclusive and interval.high_inclusive
    )


@pytest.mark.parametrize("link", LINK_FUNCTIONS)
def test_link_inverse_identity(link, global_random_seed):
    """Test that link of inverse gives identity."""
    rng = np.random.RandomState(global_random_seed)
    link = link()
    n_samples, n_classes = 100, None
    # The values for `raw_prediction` are limited from -20 to 20 because in the
    # class `LogitLink` the term `expit(x)` comes very close to 1 for large
    # positive x and therefore loses precision.
    if link.is_multiclass:
        n_classes = 10
        raw_prediction = rng.uniform(low=-20, high=20, size=(n_samples, n_classes))
        if isinstance(link, MultinomialLogit):
            raw_prediction = link.symmetrize_raw_prediction(raw_prediction)
    elif isinstance(link, HalfLogitLink):
        raw_prediction = rng.uniform(low=-10, high=10, size=(n_samples))
    else:
        raw_prediction = rng.uniform(low=-20, high=20, size=(n_samples))

    assert_allclose(link.link(link.inverse(raw_prediction)), raw_prediction)
    y_pred = link.inverse(raw_prediction)
    assert_allclose(link.inverse(link.link(y_pred)), y_pred)


@pytest.mark.parametrize(
    "namespace, device_name, dtype_name",
    yield_namespace_device_dtype_combinations(),
)
@pytest.mark.parametrize("link", LINK_FUNCTIONS)
def test_link_inverse_array_api(
    namespace, device_name, dtype_name, link, global_random_seed
):
    """Test that link and inverse link give same result for array API inputs."""
    rng = np.random.RandomState(global_random_seed)
    link = link()
    n_samples, n_classes = 100, None
    # The values for `raw_prediction` are limited from -20 to 20 because in the
    # class `LogitLink` the term `expit(x)` comes very close to 1 for large
    # positive x and therefore loses precision.
    if link.is_multiclass:
        n_classes = 10
        raw_prediction = rng.uniform(low=-20, high=20, size=(n_samples, n_classes))
        if isinstance(link, MultinomialLogit):
            raw_prediction = link.symmetrize_raw_prediction(raw_prediction)
    elif isinstance(link, HalfLogitLink):
        raw_prediction = rng.uniform(low=-10, high=10, size=(n_samples))
    else:
        raw_prediction = rng.uniform(low=-20, high=20, size=(n_samples))

    xp, device = _array_api_for_tests(namespace, device_name)
    if dtype_name != "float64":
        raw_prediction *= 0.5  # avoid overflow
        rtol = 1e-3 if n_classes else 1e-4
    else:
        rtol = 1e-8

    with config_context(array_api_dispatch=True):
        raw_prediction_xp = xp.asarray(raw_prediction.astype(dtype_name), device=device)
        assert_allclose(
            _convert_to_numpy(link.inverse(raw_prediction_xp), xp=xp),
            link.inverse(raw_prediction),
            rtol=rtol,
        )

        y_pred = link.inverse(raw_prediction)
        y_pred_xp = xp.asarray(y_pred.astype(dtype_name), device=device)
        assert_allclose(
            _convert_to_numpy(link.link(y_pred_xp), xp=xp),
            link.link(y_pred),
            rtol=rtol,
        )
