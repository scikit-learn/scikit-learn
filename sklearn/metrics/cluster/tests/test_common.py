from functools import partial
from itertools import chain

import numpy as np
import pytest

from sklearn._config import config_context
from sklearn.metrics.cluster import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    calinski_harabasz_score,
    completeness_score,
    davies_bouldin_score,
    fowlkes_mallows_score,
    homogeneity_score,
    mutual_info_score,
    normalized_mutual_info_score,
    rand_score,
    silhouette_score,
    v_measure_score,
)
from sklearn.utils._array_api import (
    _atol_for_type,
    _convert_to_numpy,
    _get_namespace_device_dtype_ids,
    yield_namespace_device_dtype_combinations,
)
from sklearn.utils._testing import _array_api_for_tests, assert_allclose

# Dictionaries of metrics
# ------------------------
# The goal of having those dictionaries is to have an easy way to call a
# particular metric and associate a name to each function:
#   - SUPERVISED_METRICS: all supervised cluster metrics - (when given a
# ground truth value)
#   - UNSUPERVISED_METRICS: all unsupervised cluster metrics
#
# Those dictionaries will be used to test systematically some invariance
# properties, e.g. invariance toward several input layout.
#

SUPERVISED_METRICS = {
    "adjusted_mutual_info_score": adjusted_mutual_info_score,
    "adjusted_rand_score": adjusted_rand_score,
    "rand_score": rand_score,
    "completeness_score": completeness_score,
    "homogeneity_score": homogeneity_score,
    "mutual_info_score": mutual_info_score,
    "normalized_mutual_info_score": normalized_mutual_info_score,
    "v_measure_score": v_measure_score,
    "fowlkes_mallows_score": fowlkes_mallows_score,
}

UNSUPERVISED_METRICS = {
    "silhouette_score": silhouette_score,
    "silhouette_manhattan": partial(silhouette_score, metric="manhattan"),
    "calinski_harabasz_score": calinski_harabasz_score,
    "davies_bouldin_score": davies_bouldin_score,
}

# Lists of metrics with common properties
# ---------------------------------------
# Lists of metrics with common properties are used to test systematically some
# functionalities and invariance, e.g. SYMMETRIC_METRICS lists all metrics
# that are symmetric with respect to their input argument y_true and y_pred.
#
# --------------------------------------------------------------------
# Symmetric with respect to their input arguments y_true and y_pred.
# Symmetric metrics only apply to supervised clusters.
SYMMETRIC_METRICS = [
    "adjusted_rand_score",
    "rand_score",
    "v_measure_score",
    "mutual_info_score",
    "adjusted_mutual_info_score",
    "normalized_mutual_info_score",
    "fowlkes_mallows_score",
]

NON_SYMMETRIC_METRICS = ["homogeneity_score", "completeness_score"]

# Metrics whose upper bound is 1
NORMALIZED_METRICS = [
    "adjusted_rand_score",
    "rand_score",
    "homogeneity_score",
    "completeness_score",
    "v_measure_score",
    "adjusted_mutual_info_score",
    "fowlkes_mallows_score",
    "normalized_mutual_info_score",
]


rng = np.random.RandomState(0)
y1 = rng.randint(3, size=30)
y2 = rng.randint(3, size=30)


def test_symmetric_non_symmetric_union():
    assert sorted(SYMMETRIC_METRICS + NON_SYMMETRIC_METRICS) == sorted(
        SUPERVISED_METRICS
    )


@pytest.mark.parametrize(
    "metric_name, y1, y2", [(name, y1, y2) for name in SYMMETRIC_METRICS]
)
def test_symmetry(metric_name, y1, y2):
    metric = SUPERVISED_METRICS[metric_name]
    assert metric(y1, y2) == pytest.approx(metric(y2, y1))


@pytest.mark.parametrize(
    "metric_name, y1, y2", [(name, y1, y2) for name in NON_SYMMETRIC_METRICS]
)
def test_non_symmetry(metric_name, y1, y2):
    metric = SUPERVISED_METRICS[metric_name]
    assert metric(y1, y2) != pytest.approx(metric(y2, y1))


@pytest.mark.parametrize("metric_name", NORMALIZED_METRICS)
def test_normalized_output(metric_name):
    upper_bound_1 = [0, 0, 0, 1, 1, 1]
    upper_bound_2 = [0, 0, 0, 1, 1, 1]
    metric = SUPERVISED_METRICS[metric_name]
    assert metric([0, 0, 0, 1, 1], [0, 0, 0, 1, 2]) > 0.0
    assert metric([0, 0, 1, 1, 2], [0, 0, 1, 1, 1]) > 0.0
    assert metric([0, 0, 0, 1, 2], [0, 1, 1, 1, 1]) < 1.0
    assert metric([0, 0, 0, 1, 2], [0, 1, 1, 1, 1]) < 1.0
    assert metric(upper_bound_1, upper_bound_2) == pytest.approx(1.0)

    lower_bound_1 = [0, 0, 0, 0, 0, 0]
    lower_bound_2 = [0, 1, 2, 3, 4, 5]
    score = np.array(
        [metric(lower_bound_1, lower_bound_2), metric(lower_bound_2, lower_bound_1)]
    )
    assert not (score < 0).any()


@pytest.mark.parametrize("metric_name", chain(SUPERVISED_METRICS, UNSUPERVISED_METRICS))
def test_permute_labels(metric_name):
    # All clustering metrics do not change score due to permutations of labels
    # that is when 0 and 1 exchanged.
    y_label = np.array([0, 0, 0, 1, 1, 0, 1])
    y_pred = np.array([1, 0, 1, 0, 1, 1, 0])
    if metric_name in SUPERVISED_METRICS:
        metric = SUPERVISED_METRICS[metric_name]
        score_1 = metric(y_pred, y_label)
        assert_allclose(score_1, metric(1 - y_pred, y_label))
        assert_allclose(score_1, metric(1 - y_pred, 1 - y_label))
        assert_allclose(score_1, metric(y_pred, 1 - y_label))
    else:
        metric = UNSUPERVISED_METRICS[metric_name]
        X = np.random.randint(10, size=(7, 10))
        score_1 = metric(X, y_pred)
        assert_allclose(score_1, metric(X, 1 - y_pred))


@pytest.mark.parametrize("metric_name", chain(SUPERVISED_METRICS, UNSUPERVISED_METRICS))
# For all clustering metrics Input parameters can be both
# in the form of arrays lists, positive, negative or string
def test_format_invariance(metric_name):
    y_true = [0, 0, 0, 0, 1, 1, 1, 1]
    y_pred = [0, 1, 2, 3, 4, 5, 6, 7]

    def generate_formats(y):
        y = np.array(y)
        yield y, "array of ints"
        yield y.tolist(), "list of ints"
        yield [str(x) + "-a" for x in y.tolist()], "list of strs"
        yield (
            np.array([str(x) + "-a" for x in y.tolist()], dtype=object),
            "array of strs",
        )
        yield y - 1, "including negative ints"
        yield y + 1, "strictly positive ints"

    if metric_name in SUPERVISED_METRICS:
        metric = SUPERVISED_METRICS[metric_name]
        score_1 = metric(y_true, y_pred)
        y_true_gen = generate_formats(y_true)
        y_pred_gen = generate_formats(y_pred)
        for (y_true_fmt, fmt_name), (y_pred_fmt, _) in zip(y_true_gen, y_pred_gen):
            assert score_1 == metric(y_true_fmt, y_pred_fmt)
    else:
        metric = UNSUPERVISED_METRICS[metric_name]
        X = np.random.randint(10, size=(8, 10))
        score_1 = metric(X, y_true)
        assert score_1 == metric(X.astype(float), y_true)
        y_true_gen = generate_formats(y_true)
        for y_true_fmt, fmt_name in y_true_gen:
            assert score_1 == metric(X, y_true_fmt)


@pytest.mark.parametrize("metric", SUPERVISED_METRICS.values())
def test_single_sample(metric):
    # only the supervised metrics support single sample
    for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        metric([i], [j])


@pytest.mark.parametrize(
    "metric_name, metric_func", dict(SUPERVISED_METRICS, **UNSUPERVISED_METRICS).items()
)
def test_inf_nan_input(metric_name, metric_func):
    if metric_name in SUPERVISED_METRICS:
        invalids = [
            ([0, 1], [np.inf, np.inf]),
            ([0, 1], [np.nan, np.nan]),
            ([0, 1], [np.nan, np.inf]),
        ]
    else:
        X = np.random.randint(10, size=(2, 10))
        invalids = [(X, [np.inf, np.inf]), (X, [np.nan, np.nan]), (X, [np.nan, np.inf])]
    with pytest.raises(ValueError, match=r"contains (NaN|infinity)"):
        for args in invalids:
            metric_func(*args)


@pytest.mark.parametrize("name", chain(SUPERVISED_METRICS, UNSUPERVISED_METRICS))
def test_returned_value_consistency(name):
    """Ensure that the returned values of all metrics are consistent.

    It can only be a float. It should not be a numpy float64 or float32.
    """

    rng = np.random.RandomState(0)
    X = rng.randint(10, size=(20, 10))
    labels_true = rng.randint(0, 3, size=(20,))
    labels_pred = rng.randint(0, 3, size=(20,))

    if name in SUPERVISED_METRICS:
        metric = SUPERVISED_METRICS[name]
        score = metric(labels_true, labels_pred)
    else:
        metric = UNSUPERVISED_METRICS[name]
        score = metric(X, labels_pred)

    assert isinstance(score, float)
    assert not isinstance(score, (np.float64, np.float32))


def check_array_api_metric(
    metric, array_namespace, device, dtype_name, a_np, b_np, **metric_kwargs
):
    xp = _array_api_for_tests(array_namespace, device)

    a_xp = xp.asarray(a_np, device=device)
    b_xp = xp.asarray(b_np, device=device)

    metric_np = metric(a_np, b_np, **metric_kwargs)

    # When array API dispatch is disabled, and np.asarray works (for example PyTorch
    # with CPU device), calling the metric function with such numpy compatible inputs
    # should work (albeit by implicitly converting to numpy arrays instead of
    # dispatching to the array library).
    try:
        np.asarray(a_xp)
        np.asarray(b_xp)
        numpy_as_array_works = True
    except (TypeError, RuntimeError, ValueError):
        # PyTorch with CUDA device and CuPy raise TypeError consistently.
        # array-api-strict chose to raise RuntimeError instead. NumPy raises
        # a ValueError if the `__array__` dunder does not return an array.
        # Exception type may need to be updated in the future for other libraries.
        numpy_as_array_works = False

    def _check_metric_matches(metric_a, metric_b, convert_a=False):
        if convert_a:
            metric_a = _convert_to_numpy(xp.asarray(metric_a), xp)
        assert_allclose(metric_a, metric_b, atol=_atol_for_type(dtype_name))

    if numpy_as_array_works:
        metric_xp = metric(a_xp, b_xp, **metric_kwargs)

        # Handle cases where multiple return values are not of the same shape,
        # e.g. precision_recall_curve:
        _check_metric_matches(metric_xp, metric_np)

        metric_xp_mixed_1 = metric(a_np, b_xp, **metric_kwargs)
        _check_metric_matches(metric_xp_mixed_1, metric_np)

        metric_xp_mixed_2 = metric(a_xp, b_np, **metric_kwargs)
        _check_metric_matches(metric_xp_mixed_2, metric_np)

    with config_context(array_api_dispatch=True):
        metric_xp = metric(a_xp, b_xp, **metric_kwargs)
        _check_metric_matches(metric_xp, metric_np, convert_a=True)


def check_array_api_unsupervised_metric(metric, array_namespace, device, dtype_name):
    y_pred = np.array([1, 0, 1, 0, 1, 1, 0])
    X = np.random.randint(10, size=(7, 10))

    check_array_api_metric(
        metric,
        array_namespace,
        device,
        dtype_name,
        a_np=X,
        b_np=y_pred,
    )


array_api_metric_checkers = {
    calinski_harabasz_score: [
        check_array_api_unsupervised_metric,
    ]
}


def yield_metric_checker_combinations(metric_checkers=array_api_metric_checkers):
    for metric, checkers in metric_checkers.items():
        for checker in checkers:
            yield metric, checker


@pytest.mark.parametrize(
    "array_namespace, device, dtype_name",
    yield_namespace_device_dtype_combinations(),
    ids=_get_namespace_device_dtype_ids,
)
@pytest.mark.parametrize("metric, check_func", yield_metric_checker_combinations())
def test_array_api_compliance(metric, array_namespace, device, dtype_name, check_func):
    check_func(metric, array_namespace, device, dtype_name)
