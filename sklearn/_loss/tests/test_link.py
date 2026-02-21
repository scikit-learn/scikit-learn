import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn._loss.link import (
    _LINKS,
    HalfLogitLink,
    Interval,
    MultinomialLogit,
    _inclusive_low_high,
)
from sklearn.utils._array_api import (
    _convert_to_numpy,
    _get_namespace_device_dtype_ids,
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
    "namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize(
    "interval",
    [
        Interval(0, 1, False, False),
        Interval(0, 1, False, True),
        Interval(0, 1, True, False),
        Interval(0, 1, True, True),
        Interval(-np.inf, np.inf, False, False),
        Interval(-np.inf, np.inf, False, True),
        Interval(-np.inf, np.inf, True, False),
        Interval(-np.inf, np.inf, True, True),
        Interval(-10, -1, False, False),
        Interval(-10, -1, False, True),
        Interval(-10, -1, True, False),
        Interval(-10, -1, True, True),
    ],
)
def test_is_in_range(namespace, device, dtype_name, interval):
    """Test that low and high are always within the interval used for linspace."""
    xp = _array_api_for_tests(namespace, device)
    dtype = xp.float32 if dtype_name == "float32" else xp.float64

    low, high = _inclusive_low_high(interval, dtype=dtype, xp=xp)

    x = xp.linspace(low, high, num=10, device=device)
    assert interval.includes(x)

    # x contains lower bound
    assert interval.includes(xp.concat((x, [interval.low]))) == interval.low_inclusive

    # x contains upper bound
    assert interval.includes(xp.concat((x, [interval.high]))) == interval.high_inclusive

    # x contains upper and lower bound
    assert interval.includes(xp.concat((x, [interval.low, interval.high]))) == (
        interval.low_inclusive and interval.high_inclusive
    )


@pytest.mark.parametrize(
    "namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize("link", LINK_FUNCTIONS)
def test_link_inverse_identity(namespace, device, dtype_name, link, global_random_seed):
    """Test that link of inverse gives identity."""
    xp = _array_api_for_tests(namespace, device)
    dtype = xp.float32 if dtype_name == "float32" else xp.float64
    rng = np.random.RandomState(global_random_seed)
    link = link()
    n_samples, n_classes = 100, None
    # The values for `raw_prediction` are limited from -20 to 20 because in the
    # class `LogitLink` the term `expit(x)` comes very close to 1 for large
    # positive x and therefore loses precision.
    if link.is_multiclass:
        n_classes = 10
        raw_prediction = rng.uniform(low=-20, high=20, size=(n_samples, n_classes))
    elif isinstance(link, HalfLogitLink):
        raw_prediction = rng.uniform(low=-10, high=10, size=(n_samples))
    else:
        raw_prediction = rng.uniform(low=-20, high=20, size=(n_samples))

    if dtype_name == "float32":
        raw_prediction *= 0.5  # avoid overflow
    raw_prediction = xp.asarray(raw_prediction, dtype=dtype, device=device)

    if isinstance(link, MultinomialLogit):
        raw_prediction = link.symmetrize_raw_prediction(raw_prediction)

    assert_allclose(
        _convert_to_numpy(link.link(link.inverse(raw_prediction)), xp=xp),
        _convert_to_numpy(raw_prediction, xp=xp),
    )
    y_pred = link.inverse(raw_prediction)
    assert_allclose(
        _convert_to_numpy(link.inverse(link.link(y_pred)), xp=xp),
        _convert_to_numpy(y_pred, xp=xp),
    )
