# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import Counter
from importlib.util import find_spec
from itertools import product

import numpy as np
import pytest

from sklearn.callback import MetricMonitor
from sklearn.callback.tests._utils import (
    MaxIterEstimator,
    MetaEstimator,
    WhileEstimator,
)
from sklearn.metrics import mean_pinball_loss, mean_squared_error


@pytest.mark.parametrize(
    "metric, metric_params",
    [(mean_squared_error, None), (mean_pinball_loss, {"alpha": 0.6})],
)
@pytest.mark.parametrize(
    "on",
    ["train_set", "validation_set", "both"],
)
def test_metric_monitor_logged_values(metric, metric_params, on):
    """Test that the correct values are logged with a simple estimator.

    The type of the log output depends on the availability of the pandas library.
    """
    max_iter = 3
    n_dim = 5
    n_samples = 3
    estimator = MaxIterEstimator(max_iter=max_iter)
    callback = MetricMonitor(metric, metric_params=metric_params, on=on)
    estimator.set_callbacks(callback)
    rng = np.random.RandomState(0)
    X_train, y_train = rng.uniform(size=(n_dim, n_samples)), rng.uniform(size=n_dim)
    X_val, y_val = rng.uniform(size=(n_dim, n_samples)), rng.uniform(size=n_dim)

    estimator.fit(X=X_train, y=y_train, X_val=X_val, y_val=y_val)
    estimator_name = estimator.__class__.__name__

    metric_params = metric_params or dict()

    expected_log_dict = {
        f"0_{estimator_name}_fit": [],
        f"1_{estimator_name}_iteration": [],
        metric.__name__: [],
        "on": [],
    }
    for i in range(max_iter):
        if on in ("train_set", "both"):
            expected_log_dict[f"0_{estimator_name}_fit"].append(0)
            expected_log_dict[f"1_{estimator_name}_iteration"].append(i)
            expected_log_dict["on"].append("train_set")
            expected_log_dict[metric.__name__].append(
                metric(y_train, X_train.mean(axis=1) * (i + 1), **metric_params)
            )
        if on in ("validation_set", "both"):
            expected_log_dict[f"0_{estimator_name}_fit"].append(0)
            expected_log_dict[f"1_{estimator_name}_iteration"].append(i)
            expected_log_dict["on"].append("validation_set")
            expected_log_dict[metric.__name__].append(
                metric(y_val, X_val.mean(axis=1) * (i + 1), **metric_params)
            )

    try:  # Check dataframes logs if pandas is installed
        import pandas as pd

        run_id, log_df = callback.logs

        expected_log_df = pd.DataFrame(expected_log_dict)
        expected_log_df = expected_log_df.set_index(
            [
                col
                for col in expected_log_df.columns
                if col not in (metric.__name__, "on")
            ]
        )

        assert log_df.equals(expected_log_df)
        assert np.array_equal(log_df.index.names, expected_log_df.index.names)

    except ImportError:  # Check dict of lists logs if pandas is not installed
        logs = callback.logs
        assert set(logs.keys()) == set(expected_log_dict.keys()).union(set(["run"]))
        for key, val in logs.items():
            if key != "run":
                assert val == expected_log_dict[key]


def test_no_predict_error():
    """Test the error raised when the estimator does not have a predict method."""
    estimator = WhileEstimator()
    callback = MetricMonitor(mean_pinball_loss, metric_params={"alpha": 0.6})
    estimator.set_callbacks(callback)

    with pytest.raises(ValueError, match="does not have a predict method"):
        estimator.fit()


def test_wrong_kwarg_error():
    """Test the error raised when giving wrong kwargs for the metric."""
    with pytest.raises(ValueError, match="cannot be used with the function"):
        MetricMonitor(mean_pinball_loss, metric_params={"wrong_name": 0.6})


@pytest.mark.parametrize("prefer", ["processes", "threads"])
@pytest.mark.parametrize(
    "metric, metric_params",
    [(mean_squared_error, None), (mean_pinball_loss, {"alpha": 0.6})],
)
@pytest.mark.parametrize(
    "on",
    ["train_set", "validation_set", "both"],
)
def test_metric_monitor_logged_values_meta_estimator(prefer, metric, metric_params, on):
    """Test that the correct values are logged with a meta-estimator.

    The type of the log output depends on the availability of the pandas library.
    """
    n_outer = 3
    n_inner = 2
    max_iter = 4
    n_dim = 5
    n_samples = 3
    rng = np.random.RandomState(0)
    X_train, y_train = rng.uniform(size=(n_dim, n_samples)), rng.uniform(size=n_dim)
    X_val, y_val = rng.uniform(size=(n_dim, n_samples)), rng.uniform(size=n_dim)
    callback = MetricMonitor(metric, metric_params=metric_params, on=on)
    est = MaxIterEstimator(max_iter=max_iter)
    est.set_callbacks(callback)
    meta_est = MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=2, prefer=prefer
    )
    meta_est_name = meta_est.__class__.__name__
    est_name = est.__class__.__name__

    meta_est.fit(X=X_train, y=y_train, X_val=X_val, y_val=y_val)

    metric_params = metric_params or dict()
    expected_log_dict = {
        metric.__name__: [],
        "on": [],
        f"0_{meta_est_name}_fit": [],
        f"1_{meta_est_name}_outer": [],
        f"2_{meta_est_name}_inner|{est_name}_fit": [],
        f"3_{est_name}_iteration": [],
    }

    for i_outer, i_inner in product(range(n_outer), range(n_inner)):
        for i_estimator_iteration in range(max_iter):
            if on in ("train_set", "both"):
                expected_log_dict[metric.__name__].append(
                    metric(
                        y_train,
                        X_train.mean(axis=1) * (i_estimator_iteration + 1),
                        **metric_params,
                    )
                )
                expected_log_dict["on"].append("train_set")
                expected_log_dict[f"0_{meta_est_name}_fit"].append(0)
                expected_log_dict[f"1_{meta_est_name}_outer"].append(i_outer)
                expected_log_dict[f"2_{meta_est_name}_inner|{est_name}_fit"].append(
                    i_inner
                )
                expected_log_dict[f"3_{est_name}_iteration"].append(
                    i_estimator_iteration
                )

            if on in ("validation_set", "both"):
                expected_log_dict[metric.__name__].append(
                    metric(
                        y_val,
                        X_val.mean(axis=1) * (i_estimator_iteration + 1),
                        **metric_params,
                    )
                )
                expected_log_dict["on"].append("validation_set")
                expected_log_dict[f"0_{meta_est_name}_fit"].append(0)
                expected_log_dict[f"1_{meta_est_name}_outer"].append(i_outer)
                expected_log_dict[f"2_{meta_est_name}_inner|{est_name}_fit"].append(
                    i_inner
                )
                expected_log_dict[f"3_{est_name}_iteration"].append(
                    i_estimator_iteration
                )

    try:  # Check dataframes logs if pandas is installed
        import pandas as pd

        run_id, log_df = callback.logs

        expected_log_df = pd.DataFrame(expected_log_dict)
        expected_log_df = expected_log_df.set_index(
            [
                col
                for col in expected_log_df.columns
                if col not in (metric.__name__, "on")
            ]
        )

        assert np.array_equal(log_df.index.names, expected_log_df.index.names)
        assert log_df.equals(expected_log_df)

    except ImportError:  # Check dict of lists logs if pandas is not installed
        logs = callback.logs
        assert set(logs.keys()) == set(expected_log_dict.keys()).union(set(["run"]))
        for key, val in logs.items():
            if key != "run":
                # Verify list equality up to a permutation because the parallelization
                # of the meta-est can change the logging order.
                assert Counter(val) == Counter(expected_log_dict[key])


def test_get_logs_output_type():
    """Test the type of the get_logs output."""
    estimator = MaxIterEstimator()
    callback = MetricMonitor(mean_squared_error)
    estimator.set_callbacks(callback)
    estimator.fit()

    if find_spec("pandas"):
        assert isinstance(callback.logs, tuple)

    else:
        assert isinstance(callback.logs, dict)

    estimator.fit()
    logs = callback.logs

    assert isinstance(logs, list)
    assert len(logs) == 2
