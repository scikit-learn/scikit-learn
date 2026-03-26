# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import re

import numpy as np
import pytest

from sklearn import config_context
from sklearn.callback import ScoringMonitor
from sklearn.callback.tests._utils import (
    MaxIterEstimator,
    MetaEstimator,
    WhileEstimator,
)
from sklearn.datasets import make_regression
from sklearn.metrics import check_scoring, make_scorer, r2_score
from sklearn.model_selection import train_test_split
from sklearn.utils._metadata_requests import UnsetMetadataPassedError


def _make_expected_output_MaxIterEstimator(
    max_iter, scoring, as_frame, X_train, y_train, X_val, y_val, depth=0
):
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
            f"estimator_name_depth_{depth}": est_name,
            f"task_name_depth_{depth}": "fit",
            f"task_id_depth_{depth}": 0,
            f"estimator_name_depth_{depth + 1}": est_name,
            f"task_name_depth_{depth + 1}": "iteration",
            f"task_id_depth_{depth + 1}": i,
        }
        for eval_on in ("train", "val"):
            X, y = (X_train, y_train) if eval_on == "train" else (X_val, y_val)
            scores = scorer(fitted_est, X, y)
            if isinstance(scores, dict):
                log_item.update({f"{eval_on}_{k}": v for k, v in scores.items()})
            else:
                log_item[f"{eval_on}_{score_names[0]}"] = scores
        expected_log.append(log_item)

    # fit root task
    log_item = {
        f"estimator_name_depth_{depth}": est_name,
        f"task_name_depth_{depth}": "fit",
        f"task_id_depth_{depth}": 0,
    }
    for eval_on in ("train", "val"):
        X, y = (X_train, y_train) if eval_on == "train" else (X_val, y_val)
        scores = scorer(fitted_est, X, y)
        if isinstance(scores, dict):
            log_item.update({f"{eval_on}_{k}": v for k, v in scores.items()})
        else:
            log_item[f"{eval_on}_{score_names[0]}"] = scores
    expected_log.append(log_item)

    if as_frame:
        pd = pytest.importorskip("pandas")
        expected_log = pd.DataFrame(expected_log)

    return expected_log


def _make_expected_output_MetaEstimator(
    n_outer, n_inner, max_iter, scoring, as_frame, X_train, y_train, X_val, y_val
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
                "estimator_name_depth_1": meta_est_name,
                "task_name_depth_1": "outer",
                "task_id_depth_1": i,
            }

            estimator_log = _make_expected_output_MaxIterEstimator(
                max_iter, scoring, False, X_train, y_train, X_val, y_val, depth=2
            )
            # root task_id of the estimator is inherited from the inner task
            for entry in estimator_log:
                entry["task_id_depth_2"] = j
                expected_log.append({**meta_est_context_levels, **entry})

    if as_frame:
        pd = pytest.importorskip("pandas")
        expected_log = pd.DataFrame(expected_log)

    return expected_log


def _custom_score_func(y_true, y_pred):
    """Custom score to test the ScoringMonitor with a callable."""
    return 0


custom_score = make_scorer(_custom_score_func)


@pytest.mark.parametrize("eval_on", ["train", "val", "both"])
def test_eval_on(eval_on):
    max_iter = 3
    callback = ScoringMonitor(eval_on=eval_on, scoring="r2")
    estimator = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)
    X, X_val, y, y_val = train_test_split(X, y, test_size=0.2, random_state=0)

    estimator.fit(X=X, y=y, X_val=X_val, y_val=y_val)
    log = callback.get_logs(as_frame=False, select="most_recent")["data"]

    # we expect one score for each iteration + the score at the end of fit
    assert len(log) == 1 + max_iter
    for row in log:
        if eval_on in ("train", "both"):
            assert "train_r2" in row
        if eval_on in ("val", "both"):
            assert "val_r2" in row


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
    callback = ScoringMonitor(eval_on="both", scoring=scoring)
    estimator = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)
    X, X_val, y, y_val = train_test_split(X, y, test_size=0.2, random_state=0)

    estimator.fit(X=X, y=y, X_val=X_val, y_val=y_val)

    log = callback.get_logs(as_frame=as_frame, select="most_recent")["data"]
    expected_log = _make_expected_output_MaxIterEstimator(
        max_iter, scoring, as_frame, X, y, X_val, y_val
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
    callback = ScoringMonitor(eval_on="both", scoring=scoring)
    est = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    meta_est = MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=2, prefer=prefer
    )
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)
    X, X_val, y, y_val = train_test_split(X, y, test_size=0.2, random_state=0)

    with config_context(enable_metadata_routing=True):
        est.set_fit_request(X_val=True, y_val=True)
        meta_est.fit(X=X, y=y, X_val=X_val, y_val=y_val)

    log = callback.get_logs(as_frame=as_frame, select="most_recent")["data"]
    expected_log = _make_expected_output_MetaEstimator(
        n_outer, n_inner, max_iter, scoring, as_frame, X, y, X_val, y_val
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
    callback = ScoringMonitor(eval_on="both", scoring="r2")
    WhileEstimator().set_callbacks(callback).fit()
    assert len(callback.get_logs(select="all")) == 0
    assert len(callback.get_logs(select="most_recent")) == 0


def test_sample_weights_and_metadata_routing():
    """Check that the ScoringMonitor works with sample weights and metadata-routing.

    - passing sample weights results in a different score than not passing them.
    - passing sample weights without metadata-routing enabled gives the same scores as
      passing them with metadata-routing enabled.
    - Not requesting sample weights gives an error if metadata-routing is enabled.
    """
    rng = np.random.RandomState(0)
    X, y = make_regression(n_samples=100, n_features=2, random_state=rng)
    sample_weight = rng.randint(0, 5, size=X.shape[0])

    # no sample weights
    callback = ScoringMonitor(eval_on="train", scoring="r2")
    MaxIterEstimator().set_callbacks(callback).fit(X=X, y=y)
    log_no_sw = callback.get_logs(as_frame=False, select="most_recent")["data"]

    # sample weights, no metadata-routing
    callback = ScoringMonitor(eval_on="train", scoring="r2")
    MaxIterEstimator().set_callbacks(callback).fit(
        X=X, y=y, sample_weight=sample_weight
    )
    log_sw_no_mr = callback.get_logs(as_frame=False, select="most_recent")["data"]

    # sample weights, metadata-routing
    with config_context(enable_metadata_routing=True):
        scorer = make_scorer(r2_score)
        scorer.set_score_request(sample_weight=True)
        callback = ScoringMonitor(eval_on="train", scoring={"r2": scorer})
        MaxIterEstimator().set_callbacks(callback).fit(
            X=X, y=y, sample_weight=sample_weight
        )
        log_sw_mr = callback.get_logs(as_frame=False, select="most_recent")["data"]

        # error if sample_weight not requested
        scorer = make_scorer(r2_score)
        callback = ScoringMonitor(eval_on="train", scoring={"r2": scorer})
        est = MaxIterEstimator().set_callbacks(callback)
        with pytest.raises(
            TypeError,
            match=re.escape("score got unexpected argument(s) {'sample_weight'}"),
        ):
            est.fit(X=X, y=y, sample_weight=sample_weight)

    assert log_no_sw != log_sw_no_mr
    assert log_sw_no_mr == log_sw_mr


def test_validation_set_metadata_routing():
    """Integration test for metadata-routing on the validation set.

    X_val and y_val must be requested for the MaxIterEstimator to be able to use them.
    """
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)
    X, X_val, y, y_val = train_test_split(X, y, test_size=0.2, random_state=0)

    n_outer, n_inner, max_iter = 2, 3, 10

    callback = ScoringMonitor(eval_on="both", scoring="r2")
    est = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)

    # Without metadata-routing enabled, passing X_val and y_val gives an error
    msg = re.escape(
        "[X_val, y_val] are passed but are not explicitly set as requested or not "
        "requested for MaxIterEstimator.fit"
    )
    with pytest.raises(UnsetMetadataPassedError, match=msg):
        MetaEstimator(est).fit(X=X, y=y, X_val=X_val, y_val=y_val)

    with config_context(enable_metadata_routing=True):
        # passing X_val and y_val without requesting them gives the same error
        with pytest.raises(UnsetMetadataPassedError, match=msg):
            MetaEstimator(est).fit(X=X, y=y, X_val=X_val, y_val=y_val)

        # with metadata-routing enabled and requested
        est.set_fit_request(X_val=True, y_val=True)
        MetaEstimator(est, n_outer=n_outer, n_inner=n_inner).fit(
            X=X, y=y, X_val=X_val, y_val=y_val
        )
        log = callback.get_logs(as_frame=False, select="most_recent")["data"]

        # 1 score for each iteration + 1 score at the end of fit in MaxIterEstimator
        assert len(log) == n_outer * n_inner * (1 + max_iter)

        # The scores on the train and validation sets are different
        assert any([entry["train_r2"] != entry["val_r2"] for entry in log])
