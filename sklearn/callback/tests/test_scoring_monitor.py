# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest

from sklearn.callback import ScoringMonitor
from sklearn.callback.tests._utils import (
    MaxIterEstimator,
    MetaEstimator,
    WhileEstimator,
)
from sklearn.datasets import make_regression
from sklearn.metrics import check_scoring, make_scorer
from sklearn.utils._testing import assert_allclose


def _make_expected_output_MaxIterEstimator(max_iter, scoring, as_frame, X, y, depth=0):
    """Generate the expected output of a ScoringMonitor on a MaxIterEstimator."""
    scorer = check_scoring(None, scoring)
    if isinstance(scoring, str):
        score_names = [scoring]
    elif callable(scoring):
        score_names = ["score"]
    else:  # multi-scorer
        score_names = list(scorer._scorers.keys())

    est_name = MaxIterEstimator.__name__

    expected_log = []
    # fit loop iterations
    for i in range(max_iter):
        fitted_est = MaxIterEstimator(max_iter=i + 1).fit()
        log_item = {
            "task_path": (0, i),
            f"estimator_name_depth_{depth}": est_name,
            f"task_name_depth_{depth}": "fit",
            f"task_id_depth_{depth}": 0,
            f"sequential_subtasks_depth_{depth}": True,
            f"estimator_name_depth_{depth + 1}": est_name,
            f"task_name_depth_{depth + 1}": f"iteration {i}",
            f"task_id_depth_{depth + 1}": i,
            f"sequential_subtasks_depth_{depth + 1}": True,
        }
        scores = scorer(fitted_est, X, y)
        if not isinstance(scores, dict):
            scores = {score_names[0]: scores}
        expected_log.append({**log_item, **scores})

    # fit root task
    log_item = {
        "task_path": (0,),
        f"estimator_name_depth_{depth}": est_name,
        f"task_name_depth_{depth}": "fit",
        f"task_id_depth_{depth}": 0,
        f"sequential_subtasks_depth_{depth}": True,
    }
    scores = scorer(fitted_est, X, y)
    if not isinstance(scores, dict):
        scores = {score_names[0]: scores}
    expected_log.append({**log_item, **scores})

    if as_frame:
        pd = pytest.importorskip("pandas")
        expected_log = pd.DataFrame(expected_log)

    return expected_log


def _make_expected_output_MetaEstimator(
    n_outer, n_inner, max_iter, scoring, as_frame, X, y
):
    """Generate the expected output of a ScoringMonitor on a MetaEstimator.

    The sub-estimators are expected to be MaxIterEstimator.
    """
    meta_est_name = MetaEstimator.__name__

    expected_log = []
    for i in range(n_outer):
        for j in range(n_inner):
            meta_est_context_levels = {
                "estimator_name_depth_0": meta_est_name,
                "task_name_depth_0": "fit",
                "task_id_depth_0": 0,
                "sequential_subtasks_depth_0": False,
                "estimator_name_depth_1": meta_est_name,
                "task_name_depth_1": "outer",
                "task_id_depth_1": i,
                "sequential_subtasks_depth_1": True,
            }

            estimator_log = _make_expected_output_MaxIterEstimator(
                max_iter, scoring, False, X, y, depth=2
            )
            for entry in estimator_log:
                # root task_id of the estimator is inherited from the inner task
                # adjust the task_id and rewrite the task_path
                entry["task_id_depth_2"] = j
                task_path = {"task_path": (0, i, j) + entry["task_path"][1:]}
                del entry["task_path"]
                expected_log.append({**task_path, **meta_est_context_levels, **entry})

    if as_frame:
        pd = pytest.importorskip("pandas")
        expected_log = pd.DataFrame(expected_log)

    return expected_log


def _custom_score_func(y_true, y_pred):
    """Custom score to test the ScoringMonitor with a callable."""
    return 0


custom_score = make_scorer(_custom_score_func)


def test_score_after_fit():
    """Check that the logged scores are the same as if the fit ended at each step."""
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    max_iter = 10
    callback = ScoringMonitor(scoring="r2")
    estimator = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    estimator.fit(X=X, y=y)

    log = callback.get_logs(as_frame=False, select="most_recent")["data"]

    scorer = check_scoring(None, "r2")
    for i in range(max_iter):
        estimator = MaxIterEstimator(max_iter=i + 1).fit(X, y)
        score = scorer(estimator, X, y)
        assert_allclose(score, log[i]["r2"])


@pytest.mark.parametrize(
    "scoring",
    ["neg_mean_squared_error", ("neg_mean_squared_error", "r2"), custom_score],
)
@pytest.mark.parametrize("as_frame", [True, False])
def test_logged_values(scoring, as_frame):
    """Test that the correct values are logged with a simple estimator."""
    if as_frame:
        pytest.importorskip("pandas")

    max_iter = 3
    callback = ScoringMonitor(scoring=scoring)
    estimator = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    estimator.fit(X=X, y=y)

    log = callback.get_logs(as_frame=as_frame, select="most_recent")["data"]
    expected_log = _make_expected_output_MaxIterEstimator(
        max_iter, scoring, as_frame, X, y
    )

    if as_frame:
        assert log.equals(expected_log)
    else:
        assert log == expected_log


@pytest.mark.parametrize("prefer", ["processes", "threads"])
@pytest.mark.parametrize(
    "scoring",
    ["neg_mean_squared_error", ("neg_mean_squared_error", "r2"), custom_score],
)
@pytest.mark.parametrize("as_frame", [True, False])
def test_logged_values_meta_estimator(prefer, scoring, as_frame):
    """Test that the correct values are logged with a meta-estimator."""
    if as_frame:
        pytest.importorskip("pandas")

    n_outer, n_inner, max_iter = 3, 2, 5
    callback = ScoringMonitor(scoring=scoring)
    est = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    meta_est = MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=2, prefer=prefer
    )
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    meta_est.fit(X=X, y=y)

    log = callback.get_logs(as_frame=as_frame, select="most_recent")["data"]
    expected_log = _make_expected_output_MetaEstimator(
        n_outer, n_inner, max_iter, scoring, as_frame, X, y
    )

    # The log items might not be in the same order as the expected log due to the
    # parallelization of the meta-estimator. hence we sort the log by id of the
    # parallel task, "outer", which at depth 1.
    if as_frame:
        log.sort_values(
            by=["task_id_depth_1"], inplace=True, ignore_index=True, kind="stable"
        )
        assert log.equals(expected_log)
    else:
        log.sort(key=lambda x: x["task_id_depth_1"])
        assert log == expected_log


@pytest.mark.parametrize("as_frame", [True, False])
def test_get_logs_output_type_no_fit(as_frame):
    """Check that get_logs return empty containers of the right type before fit."""
    if as_frame:
        pytest.importorskip("pandas")

    callback = ScoringMonitor(scoring="neg_mean_squared_error")

    # "all" logs is always a list of run dicts.
    logs_all = callback.get_logs(select="all", as_frame=as_frame)
    assert isinstance(logs_all, list)
    assert len(logs_all) == 0

    log_most_recent = callback.get_logs(select="most_recent", as_frame=as_frame)
    assert isinstance(log_most_recent, dict)
    assert len(log_most_recent) == 0


@pytest.mark.parametrize("as_frame", [True, False])
def test_get_logs_output_type(as_frame):
    """Test the type of the get_logs output."""
    if as_frame:
        pd = pytest.importorskip("pandas")

    callback = ScoringMonitor(scoring="neg_mean_squared_error")
    estimator = MaxIterEstimator().set_callbacks(callback)
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    # 2 runs
    estimator.fit(X, y)
    estimator.fit(X, y)

    logs_all = callback.get_logs(select="all", as_frame=as_frame)
    assert isinstance(logs_all, list)
    assert len(logs_all) == 2
    assert all(isinstance(log, dict) for log in logs_all)

    expected_data_type = pd.DataFrame if as_frame else list
    assert all(isinstance(log["data"], expected_data_type) for log in logs_all)

    log_most_recent = callback.get_logs(select="most_recent", as_frame=as_frame)
    assert isinstance(log_most_recent, dict)
    assert isinstance(log_most_recent["data"], expected_data_type)


def test_estimator_without_reconstruction_attributes():
    """Smoke test on an estimator which does not provide reconstruction_attributes."""
    callback = ScoringMonitor(scoring="r2")
    WhileEstimator().set_callbacks(callback).fit()
    assert len(callback.get_logs(select="all")) == 0
    assert len(callback.get_logs(select="most_recent")) == 0


def test_scoringmonitor_sample_weights():
    """Check that the ScoringMonitor works with sample weights."""
    rng = np.random.RandomState(0)
    X, y = make_regression(n_samples=100, n_features=2, random_state=rng)
    sample_weight = rng.randint(0, 5, size=X.shape[0])

    # no sample weights
    callback = ScoringMonitor(scoring="r2")
    MaxIterEstimator().set_callbacks(callback).fit(X=X, y=y)
    log_no_sw = callback.get_logs(as_frame=False, select="most_recent")["data"]

    # sample weights
    callback = ScoringMonitor(scoring="r2")
    MaxIterEstimator().set_callbacks(callback).fit(
        X=X, y=y, sample_weight=sample_weight
    )
    log_sw = callback.get_logs(as_frame=False, select="most_recent")["data"]

    assert log_no_sw != log_sw
