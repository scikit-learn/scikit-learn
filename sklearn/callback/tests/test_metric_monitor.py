# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from itertools import product

import numpy as np
import pandas as pd
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
def test_metric_monitor(metric, metric_params):
    max_iter = 3
    n_dim = 5
    n_samples = 3
    intercept = 1
    estimator = MaxIterEstimator(intercept=intercept, max_iter=max_iter)
    callback_train = MetricMonitor(
        metric, metric_params=metric_params, on_validation=False
    )
    callback_val = MetricMonitor(
        metric, metric_params=metric_params, on_validation=True
    )
    estimator.set_callbacks([callback_train, callback_val])
    rng = np.random.RandomState(0)
    X_train, y_train = rng.uniform(size=(n_dim, n_samples)), rng.uniform(size=n_dim)
    X_val, y_val = rng.uniform(size=(n_dim, n_samples)), rng.uniform(size=n_dim)

    estimator.fit(X=X_train, y=y_train, X_val=X_val, y_val=y_val)

    metric_params = metric_params or dict()
    log_train = callback_train.get_logs()
    assert len(log_train) == 1
    run_id_train, log_train = next(iter(log_train.items()))
    assert f"{estimator.__class__.__name__}_{id(estimator)}_" in run_id_train

    expected_log_train = pd.DataFrame(
        [
            {
                f"0_{estimator.__class__.__name__}_fit": 0,
                f"1_{estimator.__class__.__name__}_fit_iter": i,
                metric.__name__: metric(
                    y_train, X_train.mean(axis=1) * (i + 1) + intercept, **metric_params
                ),
            }
            for i in range(max_iter)
        ]
    )
    expected_log_train = expected_log_train.set_index(
        [col for col in expected_log_train.columns if col != metric.__name__]
    )
    assert log_train.equals(expected_log_train)
    assert np.array_equal(log_train.index.names, expected_log_train.index.names)

    log_val = callback_val.get_logs()
    assert len(log_val) == 1
    run_id_val, log_val = next(iter(log_val.items()))
    assert f"{estimator.__class__.__name__}_{id(estimator)}_" in run_id_val

    expected_log_val = pd.DataFrame(
        [
            {
                f"0_{estimator.__class__.__name__}_fit": 0,
                f"1_{estimator.__class__.__name__}_fit_iter": i,
                metric.__name__: metric(
                    y_val, X_val.mean(axis=1) * (i + 1) + intercept, **metric_params
                ),
            }
            for i in range(max_iter)
        ]
    )
    expected_log_val = expected_log_val.set_index(
        [col for col in expected_log_val.columns if col != metric.__name__]
    )
    assert log_val.equals(expected_log_val)
    assert np.array_equal(log_val.index.names, expected_log_val.index.names)


def test_no_predict_error():
    estimator = WhileEstimator()
    callback = MetricMonitor(mean_pinball_loss, metric_params={"alpha": 0.6})
    estimator.set_callbacks(callback)

    with pytest.raises(ValueError, match="does not have a predict method"):
        estimator.fit()


def test_wrong_kwarg_error():
    with pytest.raises(ValueError, match="cannot be used with the function"):
        MetricMonitor(mean_pinball_loss, metric_params={"wrong_name": 0.6})


@pytest.mark.parametrize("prefer", ["processes", "threads"])
@pytest.mark.parametrize(
    "metric, metric_params",
    [(mean_squared_error, None), (mean_pinball_loss, {"alpha": 0.6})],
)
def test_within_meta_estimator(prefer, metric, metric_params):
    n_outer = 3
    n_inner = 2
    max_iter = 4
    n_dim = 5
    n_samples = 3
    rng = np.random.RandomState(0)
    X_train, y_train = rng.uniform(size=(n_dim, n_samples)), rng.uniform(size=n_dim)
    X_val, y_val = rng.uniform(size=(n_dim, n_samples)), rng.uniform(size=n_dim)
    callback_train = MetricMonitor(
        metric, metric_params=metric_params, on_validation=False
    )
    callback_val = MetricMonitor(
        metric, metric_params=metric_params, on_validation=True
    )
    est = MaxIterEstimator(max_iter=max_iter)
    est.set_callbacks([callback_train, callback_val])
    meta_est = MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=2, prefer=prefer
    )

    meta_est.fit(X=X_train, y=y_train, X_val=X_val, y_val=y_val)

    metric_params = metric_params or dict()
    expected_log_train = []
    expected_log_val = []
    for i_outer, i_inner in product(range(n_outer), range(n_inner)):
        est = MaxIterEstimator()
        for i_estimator_fit_iter in range(max_iter):
            setattr(est, "n_iter_", i_estimator_fit_iter + 1)
            expected_log_train.append(
                {
                    metric.__name__: metric(
                        y_train, est.predict(X_train), **metric_params
                    ),
                    f"0_{meta_est.__class__.__name__}_fit": 0,
                    f"1_{meta_est.__class__.__name__}_outer": i_outer,
                    f"2_{meta_est.__class__.__name__}_inner|"
                    f"{est.__class__.__name__}_fit": i_inner,
                    f"3_{est.__class__.__name__}_fit_iter": i_estimator_fit_iter,
                }
            )
            expected_log_val.append(
                {
                    metric.__name__: metric(y_val, est.predict(X_val), **metric_params),
                    f"0_{meta_est.__class__.__name__}_fit": 0,
                    f"1_{meta_est.__class__.__name__}_outer": i_outer,
                    f"2_{meta_est.__class__.__name__}_inner|"
                    f"{est.__class__.__name__}_fit": i_inner,
                    f"3_{est.__class__.__name__}_fit_iter": i_estimator_fit_iter,
                }
            )
    expected_log_train = pd.DataFrame(expected_log_train)
    expected_log_train = expected_log_train.set_index(
        [col for col in expected_log_train.columns if col != metric.__name__]
    )
    expected_log_val = pd.DataFrame(expected_log_val)
    expected_log_val = expected_log_val.set_index(
        [col for col in expected_log_val.columns if col != metric.__name__]
    )

    log_train = callback_train.get_logs()
    assert len(log_train) == 1
    run_id_train, log_train = next(iter(log_train.items()))
    log_val = callback_val.get_logs()
    assert len(log_val) == 1
    run_id_val, log_val = next(iter(log_val.items()))

    assert f"{meta_est.__class__.__name__}_{id(meta_est)}_" in run_id_train
    assert f"{meta_est.__class__.__name__}_{id(meta_est)}_" in run_id_val
    assert len(log_train) == len(expected_log_train)
    assert len(log_val) == len(expected_log_val)
    assert np.array_equal(log_train.index.names, expected_log_train.index.names)
    assert np.array_equal(log_val.index.names, expected_log_val.index.names)
    assert log_train.equals(expected_log_train)
    assert log_val.equals(expected_log_val)
