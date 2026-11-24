# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import gc
import re

import numpy as np
import pytest

from sklearn import config_context
from sklearn.base import clone
from sklearn.callback import ScoringMonitor
from sklearn.callback._scoring_monitor import ScoringMonitorLog
from sklearn.callback._testing.estimators import (
    MaxIterEstimator,
    MetaEstimator,
    WhileEstimator,
)
from sklearn.callback._transport import _listeners
from sklearn.datasets import make_regression
from sklearn.metrics import check_scoring, get_scorer, make_scorer
from sklearn.utils._testing import assert_allclose
from sklearn.utils.parallel import Parallel, delayed


def _get_score_names(scoring, scorer):
    if isinstance(scoring, str):
        return [scoring]
    elif callable(scoring):
        return ["score"]
    else:  # multi-scorer
        return list(scorer._scorers.keys())


def _make_expected_output_MaxIterEstimator(
    max_iter, scoring, as_pandas, X, y, sample_weight=None, root_task_id=0
):
    """Generate the expected output of a ScoringMonitor on a MaxIterEstimator."""
    scorer = (
        check_scoring(None, scoring)
        if scoring not in ("no_train_score", "no_val_score")
        else None
    )
    score_names = _get_score_names(scoring, scorer)

    est_name = MaxIterEstimator.__name__

    expected_log = []

    # fit root task
    fitted_est = MaxIterEstimator(max_iter=max_iter).fit()
    log_item = {
        "task_id_path": (root_task_id,),
        "parent_task_id_path": (),
        "estimator_name": est_name,
        "task_name": "fit",
        "task_id": root_task_id,
        "sequential_subtasks": True,
    }
    if scorer is not None:
        scores = (
            scorer(fitted_est, X, y, sample_weight=sample_weight)
            if scorer._accept_sample_weight()
            else scorer(fitted_est, X, y)
        )
        if not isinstance(scores, dict):
            scores = {score_names[0]: scores}
        expected_log.append({**log_item, **scores})

    # fit loop iterations
    for i in range(max_iter):
        fitted_est = MaxIterEstimator(max_iter=i + 1).fit()
        log_item = {
            "task_id_path": (root_task_id, i),
            "parent_task_id_path": (root_task_id,),
            "estimator_name": est_name,
            "task_name": f"iteration {i}",
            "task_id": i,
            "sequential_subtasks": True,
        }
        if scorer is not None:
            scores = (
                scorer(fitted_est, X, y, sample_weight=sample_weight)
                if scorer._accept_sample_weight()
                else scorer(fitted_est, X, y)
            )
            if not isinstance(scores, dict):
                scores = {score_names[0]: scores}
            expected_log.append({**log_item, **scores})

    if as_pandas:
        pd = pytest.importorskip("pandas")
        expected_log = pd.DataFrame(expected_log)

    return expected_log


def _make_expected_output_MetaEstimator(
    n_outer,
    n_inner,
    max_iter,
    scoring,
    as_pandas,
    include_lineage,
    X,
    y,
    sample_weight=None,
):
    """Generate the expected output of a ScoringMonitor on a MetaEstimator.

    The sub-estimators are expected to be MaxIterEstimator.
    """
    meta_est_name = MetaEstimator.__name__

    expected_log = []
    if scoring in ("no_train_score", "no_val_score"):
        if as_pandas:
            pd = pytest.importorskip("pandas")
            expected_log = pd.DataFrame(expected_log)
        return expected_log

    for i in range(n_outer):
        for j in range(n_inner):
            estimator_log = _make_expected_output_MaxIterEstimator(
                max_iter,
                scoring,
                False,
                X,
                y,
                root_task_id=j,
                sample_weight=sample_weight,
            )
            for row in estimator_log:
                row["task_id_path"] = (0, i) + row["task_id_path"]
                row["parent_task_id_path"] = (0, i) + row["parent_task_id_path"]
                expected_log.append(row)

    if include_lineage:
        # fit root task
        expected_log.append(
            {
                "task_id_path": (0,),
                "parent_task_id_path": (),
                "estimator_name": meta_est_name,
                "task_name": "fit",
                "task_id": 0,
                "sequential_subtasks": False,
            }
        )
        # outer tasks
        for i in range(n_outer):
            expected_log.append(
                {
                    "task_id_path": (0, i),
                    "parent_task_id_path": (0,),
                    "estimator_name": meta_est_name,
                    "task_name": "outer",
                    "task_id": i,
                    "sequential_subtasks": True,
                }
            )

    sorting_key = lambda x: (len(x["parent_task_id_path"]), x["parent_task_id_path"])
    expected_log.sort(key=sorting_key)

    if as_pandas:
        pd = pytest.importorskip("pandas")
        expected_log = pd.DataFrame(expected_log)

    return expected_log


def _custom_score_func(y_true, y_pred):
    """Custom score to test the ScoringMonitor with a scorer made by `make_scorer`."""
    return 0


custom_score = make_scorer(_custom_score_func)


def _fit_for_logged_values(estimator, scoring_val, with_sample_weight):
    """Generates the inputs and fits the estimator for testing the logged values.

    Helper function for test_logged_values and test_logged_values_meta_estimator.
    Returns the training data, the optional validation data, and the sample
    weights passed to `fit`, so the expected logs can be built from the same
    inputs.
    """
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)
    if scoring_val != "no_val_score":
        X_val, y_val = make_regression(n_samples=10, n_features=2, random_state=1)
    else:
        X_val, y_val = None, None

    rng = np.random.RandomState(0)
    sample_weight = rng.randint(1, 5, size=X.shape[0]) if with_sample_weight else None
    sample_weight_val = (
        rng.randint(1, 5, size=X_val.shape[0])
        if with_sample_weight and scoring_val != "no_val_score"
        else None
    )

    if scoring_val != "no_val_score":
        with config_context(enable_metadata_routing=True):
            if with_sample_weight:
                estimator.fit(
                    X=X,
                    y=y,
                    X_val=X_val,
                    y_val=y_val,
                    sample_weight=sample_weight,
                    sample_weight_val=sample_weight_val,
                )
            else:
                estimator.fit(X=X, y=y, X_val=X_val, y_val=y_val)
    else:
        if with_sample_weight:
            estimator.fit(X=X, y=y, sample_weight=sample_weight)
        else:
            estimator.fit(X=X, y=y)

    return X, y, X_val, y_val, sample_weight, sample_weight_val


def test_score_after_fit():
    """Check that the logged scores are the same as if the fit ended at each step."""
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    max_iter = 10
    callback = ScoringMonitor(scoring_train="r2")
    estimator = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    estimator.fit(X=X, y=y)

    log = callback.get_logs()
    # select the rows corresponding to the iteration tasks
    iter_log = [
        row for row in log.train_scores if row["task_name"].startswith("iteration")
    ]

    scorer = check_scoring(None, "r2")
    for i in range(max_iter):
        estimator = MaxIterEstimator(max_iter=i + 1).fit(X, y)
        score = scorer(estimator, X, y)
        assert_allclose(score, iter_log[i]["r2"])


@pytest.mark.parametrize(
    "scoring_train",
    [
        "neg_mean_squared_error",
        ("neg_mean_squared_error", "r2"),
        custom_score,
        "no_train_score",
    ],
)
@pytest.mark.parametrize("scoring_val", ["neg_mean_squared_error", "no_val_score"])
@pytest.mark.parametrize("with_sample_weight", [True, False])
@pytest.mark.parametrize("as_pandas", [True, False])
def test_logged_values(scoring_train, scoring_val, with_sample_weight, as_pandas):
    """Test that the correct values are logged with a simple estimator."""
    if as_pandas:
        pytest.importorskip("pandas")

    if scoring_train == "no_train_score" and scoring_val == "no_val_score":
        pytest.skip("At least one scorer must be set.")

    if scoring_train is custom_score and with_sample_weight:
        pytest.skip("Custom score does not accept sample weight.")

    max_iter = 3
    callback = ScoringMonitor(scoring_train=scoring_train, scoring_val=scoring_val)
    estimator = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    X, y, X_val, y_val, sample_weight, sample_weight_val = _fit_for_logged_values(
        estimator, scoring_val, with_sample_weight
    )

    log = callback.get_logs()
    log_train = getattr(
        log, "train_scores" if not as_pandas else "train_scores_as_pandas"
    )
    expected_log_train = _make_expected_output_MaxIterEstimator(
        max_iter, scoring_train, as_pandas, X, y, sample_weight
    )

    log_val = getattr(log, "val_scores" if not as_pandas else "val_scores_as_pandas")
    expected_log_val = _make_expected_output_MaxIterEstimator(
        max_iter, scoring_val, as_pandas, X_val, y_val, sample_weight_val
    )

    if as_pandas:
        assert log_train.equals(expected_log_train)
        assert log_val.equals(expected_log_val)
    else:
        assert log_train == expected_log_train
        assert log_val == expected_log_val


@pytest.mark.parametrize("prefer", ["processes", "threads"])
@pytest.mark.parametrize(
    "scoring_train",
    ["neg_mean_squared_error", ("neg_mean_squared_error", "r2"), custom_score],
)
@pytest.mark.parametrize("scoring_val", ["neg_mean_squared_error", "no_val_score"])
@pytest.mark.parametrize("as_pandas", [True, False])
@pytest.mark.parametrize("include_lineage", [True, False])
@pytest.mark.parametrize("with_sample_weight", [True, False])
def test_logged_values_meta_estimator(
    prefer, scoring_train, scoring_val, as_pandas, include_lineage, with_sample_weight
):
    """Test that the correct values are logged with a meta-estimator."""
    if as_pandas:
        pytest.importorskip("pandas")

    n_outer, n_inner, max_iter = 3, 2, 5
    callback = ScoringMonitor(scoring_train=scoring_train, scoring_val=scoring_val)
    est = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    meta_est = MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=2, prefer=prefer
    )
    X, y, X_val, y_val, sample_weight, sample_weight_val = _fit_for_logged_values(
        meta_est, scoring_val, with_sample_weight
    )

    log = callback.get_logs(include_lineage=include_lineage)
    log_train = log.train_scores_as_pandas if as_pandas else log.train_scores
    log_val = log.val_scores_as_pandas if as_pandas else log.val_scores
    expected_log_train = _make_expected_output_MetaEstimator(
        n_outer,
        n_inner,
        max_iter,
        scoring_train,
        as_pandas,
        include_lineage,
        X,
        y,
        sample_weight,
    )
    expected_log_val = _make_expected_output_MetaEstimator(
        n_outer,
        n_inner,
        max_iter,
        scoring_val,
        as_pandas,
        include_lineage,
        X_val,
        y_val,
        sample_weight_val,
    )

    if as_pandas:
        assert log_train.equals(expected_log_train)
        assert log_val.equals(expected_log_val)
    else:
        assert log_train == expected_log_train
        assert log_val == expected_log_val


@pytest.mark.parametrize("select", ["all", "most_recent"])
def test_get_logs_output_type_no_fit(select):
    """Check that get_logs raises an error before fit."""
    callback = ScoringMonitor(scoring_train="neg_mean_squared_error")

    with pytest.raises(ValueError, match="No logs to retrieve"):
        callback.get_logs(select=select)


@pytest.mark.parametrize("as_pandas", [True, False])
def test_get_logs_output_type(as_pandas):
    """Test the type of the get_logs output."""
    if as_pandas:
        pd = pytest.importorskip("pandas")

    callback = ScoringMonitor(scoring_train="neg_mean_squared_error")
    estimator = MaxIterEstimator().set_callbacks(callback)
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    # 2 runs
    estimator.fit(X, y)
    estimator.fit(X, y)

    logs_all = callback.get_logs(select="all")
    assert isinstance(logs_all, list)
    assert len(logs_all) == 2
    assert all(isinstance(log, ScoringMonitorLog) for log in logs_all)

    expected_data_type = pd.DataFrame if as_pandas else list
    attr = "train_scores" if not as_pandas else "train_scores_as_pandas"
    assert all(isinstance(getattr(log, attr), expected_data_type) for log in logs_all)

    log_most_recent = callback.get_logs(select="most_recent")
    assert isinstance(log_most_recent, ScoringMonitorLog)
    assert isinstance(getattr(log_most_recent, attr), expected_data_type)


@pytest.mark.parametrize("select", ["all", "most_recent"])
def test_estimator_without_reconstruction_attributes(select):
    """Smoke test on an estimator which does not provide reconstruction_attributes."""
    callback = ScoringMonitor(scoring_train="r2")
    WhileEstimator().set_callbacks(callback).fit()

    with pytest.raises(ValueError, match="No logs to retrieve"):
        callback.get_logs(select=select)


def test_scoring_monitor_no_metadata_routing_estimator():
    """Smoke test for estimators that don't implement metadata routing."""
    X, y = make_regression(n_samples=10, n_features=2, random_state=0)
    cb = ScoringMonitor(scoring_train="r2")
    est = WhileEstimator().set_callbacks(cb)
    est.fit(X, y)
    expected_log_length = 22  # WhileEstimator does 21 iterations + fit task
    assert len(cb.get_logs().train_scores) == expected_log_length


@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_scoringmonitor_metadata_routing():
    """Test the routing of metadata to the scorer."""
    score_func = lambda y_true, y_pred, req_arg: req_arg
    scorer_train = make_scorer(score_func).set_score_request(req_arg=True)
    scorer_val = make_scorer(score_func).set_score_request(req_arg="req_arg_val")
    cb = ScoringMonitor(scoring_train=scorer_train, scoring_val=scorer_val)
    X, y = make_regression(n_samples=10, n_features=2, random_state=0)
    MaxIterEstimator(max_iter=2).set_callbacks(cb).fit(
        X, y, X_val=X, y_val=y, req_arg="train_value", req_arg_val="val_value"
    )
    logs = cb.get_logs()
    assert logs.train_scores
    assert all([ll["score"] == "train_value" for ll in logs.train_scores])
    assert logs.val_scores
    assert all([ll["score"] == "val_value" for ll in logs.val_scores])


@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_scoringmonitor_metadata_routing_meta_estimator():
    """Test the routing of metadata to the scorer in a meta-estimator."""
    score_func = lambda y_true, y_pred, req_arg: req_arg
    scorer_train = make_scorer(score_func).set_score_request(req_arg=True)
    scorer_val = make_scorer(score_func).set_score_request(req_arg="req_arg_val")
    cb = ScoringMonitor(scoring_train=scorer_train, scoring_val=scorer_val)
    X, y = make_regression(n_samples=10, n_features=2, random_state=0)
    est = MaxIterEstimator().set_callbacks(cb)
    MetaEstimator(estimator=est).fit(
        X, y, X_val=X, y_val=y, req_arg="train_value", req_arg_val="val_value"
    )
    logs = cb.get_logs()
    assert logs.train_scores
    assert all([ll["score"] == "train_value" for ll in logs.train_scores])
    assert logs.val_scores
    assert all([ll["score"] == "val_value" for ll in logs.val_scores])


def _get_ancestors_info_path(log, parent_task_id_path):
    log = {row["task_id_path"]: row for row in log}
    ancestors_info_path = []
    while parent_task_id_path:
        task_info = log[parent_task_id_path]
        ancestors_info_path.append(
            {
                "estimator_name": task_info["estimator_name"],
                "task_name": task_info["task_name"],
                "task_id": task_info["task_id"],
            }
        )
        parent_task_id_path = task_info["parent_task_id_path"]
    ancestors_info_path.reverse()
    return ancestors_info_path


def _get_ancestors_info_path_pandas(log, parent_task_id_path):
    log = log.set_index("task_id_path", drop=False)
    ancestors_info_path = []
    while parent_task_id_path:
        task_info = log.loc[[parent_task_id_path]].iloc[0]
        ancestors_info_path.append(
            {
                "estimator_name": task_info["estimator_name"],
                "task_name": task_info["task_name"],
                "task_id": task_info["task_id"],
            }
        )
        parent_task_id_path = task_info["parent_task_id_path"]
    ancestors_info_path.reverse()
    return ancestors_info_path


@pytest.mark.parametrize("as_pandas", [True, False])
def test_get_logs_include_lineage_ancestor_retrieval(as_pandas):
    """Check that ancestor task info can be retrieved from a given task."""
    if as_pandas:
        pytest.importorskip("pandas")

    n_outer, n_inner, max_iter = 2, 3, 5
    callback = ScoringMonitor(scoring_train="r2")
    est = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    meta_est = MetaEstimator(est, n_outer=n_outer, n_inner=n_inner)
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    meta_est.fit(X=X, y=y)

    attr = "train_scores" if not as_pandas else "train_scores_as_pandas"
    log = callback.get_logs(include_lineage=True)
    log = getattr(log, attr)

    # pick an arbitrary parent of iteration tasks
    parent_path = (0, 1, 1)
    expected_ancestors_info_path = [
        {
            "estimator_name": "MetaEstimator",
            "task_name": "fit",
            "task_id": 0,
        },
        {
            "estimator_name": "MetaEstimator",
            "task_name": "outer",
            "task_id": 1,
        },
        {
            "estimator_name": "MaxIterEstimator",
            "task_name": "fit",
            "task_id": 1,
        },
    ]

    if as_pandas:
        sibling_rows = log.loc[log["parent_task_id_path"] == parent_path]
        ancestors_info_path = _get_ancestors_info_path_pandas(log, parent_path)
    else:
        sibling_rows = [row for row in log if row["parent_task_id_path"] == parent_path]
        ancestors_info_path = _get_ancestors_info_path(log, parent_path)

    assert len(sibling_rows) == max_iter
    assert ancestors_info_path == expected_ancestors_info_path


@pytest.mark.parametrize("backend", ["threading", "loky"])
def test_scoring_monitor_no_callback_support(backend):
    """Test ScoringMonitor within a parallelized function not supporting callbacks.

    The expected result is multiple separate logs, similar to the result of multiple
    sequential fit executions.
    """

    def clone_and_fit(estimator, X, y):
        clone(estimator).fit(X, y)

    def func(estimator, X, y, n_fits):
        Parallel(n_jobs=2, backend=backend)(
            delayed(clone_and_fit)(estimator, X, y) for _ in range(n_fits)
        )

    max_iter = 3
    n_fits = 4

    callback = ScoringMonitor(scoring_train="r2")
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)
    func(
        MaxIterEstimator(max_iter=max_iter).set_callbacks(callback), X, y, n_fits=n_fits
    )
    log_all = callback.get_logs(select="all")

    expected_log = _make_expected_output_MaxIterEstimator(
        max_iter, scoring="r2", as_pandas=False, X=X, y=y
    )

    assert len(log_all) == n_fits
    for log in log_all:
        assert log.train_scores == expected_log
        assert log.estimator_name == "MaxIterEstimator"


def test_scoring_monitor_listener_closed_on_gc():
    """Check that listener is closed when callback is garbage collected."""
    callback = ScoringMonitor(scoring_train="r2")
    listener_address = callback._channel._listener_address
    assert listener_address in _listeners

    MaxIterEstimator().set_callbacks(callback).fit()

    del callback
    gc.collect()

    assert listener_address not in _listeners


def test_no_scorer_error():
    """Test the error when not setting a scorer."""
    with pytest.raises(
        ValueError,
        match="`scoring_train='no_train_score'` and `scoring_val='no_val_score'`",
    ):
        ScoringMonitor()


def test_scoring_val_no_metadata_routing_error():
    """Test the error when using a validation scorer with metadata routing disabled."""
    est = MaxIterEstimator().set_callbacks(ScoringMonitor(scoring_val="r2"))
    with pytest.raises(
        ValueError,
        match=r"scorer on validation data .* only supported when metadata routing",
    ):
        est.fit()


def test_repr_scoringmonitorlog():
    """Test the repr of the ScoringMonitorLog."""
    cb = ScoringMonitor(scoring_train="r2")
    X, y = make_regression(n_samples=10, n_features=2)
    MaxIterEstimator().set_callbacks(cb).fit(X, y)
    pattern = (
        r"ScoringMonitorLog\(run_id=.*, estimator_name=MaxIterEstimator, timestamp=.*\)"
    )
    assert re.match(pattern, repr(cb.get_logs()))


def test_accept_sample_weight():
    """Test the _accept_sample_weight method."""
    cb = ScoringMonitor(scoring_train="r2")
    # True when the scorer accepts sample_weight.
    assert cb._accept_sample_weight("on_fit_task_end")
    # False on wrong hook.
    assert not cb._accept_sample_weight("on_fit_task_begin")

    # False with a warning when the scorer does not accept sample_weight.
    with pytest.warns(UserWarning, match="does not support sample_weight"):
        assert not ScoringMonitor(scoring_train=custom_score)._accept_sample_weight(
            "on_fit_task_end"
        )

    # True with a warning when a multiscorer with at least one scorer that accepts
    # sample_weiht but at least one that does not.
    with pytest.warns(UserWarning, match="does not support sample_weight"):
        assert ScoringMonitor(scoring_train=["r2", "rand_score"])._accept_sample_weight(
            "on_fit_task_end"
        )

    # False with a warning for a callable, even if it has a sample_weight parameter,
    # because it cannot be routed any metadata.
    scoring = lambda estimator, X, y, sample_weight=None: 1.0
    with pytest.warns(UserWarning, match="does not support sample_weight"):
        assert not ScoringMonitor(scoring_train=scoring)._accept_sample_weight(
            "on_fit_task_end"
        )


@pytest.mark.parametrize(
    "scoring, expected_scores",
    [
        (lambda estimator, X, y: 1.0, {"score": 1.0}),
        (lambda estimator, X, y: {"a": 1.0, "b": 2.0}, {"a": 1.0, "b": 2.0}),
    ],
)
@pytest.mark.parametrize("dataset", ["train", "val"])
def test_callable_scoring(scoring, expected_scores, dataset):
    """Check that a callable (not made by `make_scorer`) can be used as scorer.

    It can return either a single value or a dict of values, and be used on the
    training or the validation data.
    """
    X, y = make_regression(n_samples=10, n_features=2, random_state=0)
    fit_params = {}
    if dataset == "val":
        X_val, y_val = make_regression(n_samples=10, n_features=2, random_state=1)
        fit_params = {"X_val": X_val, "y_val": y_val}

    # Scoring on the validation data requires metadata routing.
    with config_context(
        enable_metadata_routing=dataset == "val", enable_metadata_auto_requests=True
    ):
        cb = ScoringMonitor(**{f"scoring_{dataset}": scoring})
        MaxIterEstimator().set_callbacks(cb).fit(X, y, **fit_params)

    scores = getattr(cb.get_logs(), f"{dataset}_scores")
    assert scores
    for row in scores:
        assert {k: row[k] for k in expected_scores} == expected_scores


@config_context(enable_metadata_routing=True)
def test_val_scorer_with_aliased_data_request():
    """Test that the validation scorer with aliased data request is working."""
    cb_no_alias = ScoringMonitor(scoring_val="r2")
    cb_alias = ScoringMonitor(scoring_val="r2")
    X, y = make_regression(n_samples=10, n_features=2)
    X_val_1, y_val_1 = make_regression(n_samples=10, n_features=2, random_state=0)
    X_val_2, y_val_2 = make_regression(n_samples=10, n_features=2, random_state=1)
    cb_alias.set_on_fit_task_end_request(X_val="X_val_alias", y_val="y_val_alias")
    MaxIterEstimator().set_callbacks(cb_no_alias, cb_alias).fit(
        X, y, X_val=X_val_1, y_val=y_val_1, X_val_alias=X_val_2, y_val_alias=y_val_2
    )
    alias_score_values = [log_item["r2"] for log_item in cb_alias.get_logs().val_scores]
    no_alias_score_values = [
        log_item["r2"] for log_item in cb_no_alias.get_logs().val_scores
    ]
    # Verify that the callback with alias request is not using the same data
    assert not np.allclose(alias_score_values, no_alias_score_values)

    cb_no_alias = ScoringMonitor(scoring_val="r2")
    MaxIterEstimator().set_callbacks(cb_no_alias).fit(
        X, y, X_val=X_val_2, y_val=y_val_2
    )
    no_alias_score_values = [
        log_item["r2"] for log_item in cb_no_alias.get_logs().val_scores
    ]
    # Verify that with or without alias give the same score with the same data
    assert_allclose(no_alias_score_values, alias_score_values)


@config_context(enable_metadata_routing=True)
def test_val_scorer_with_aliased_sample_weight_request():
    """Test that the validation scorer with aliased sample weight request is working."""
    cb_no_alias = ScoringMonitor(scoring_val="r2")
    scorer = get_scorer("r2").set_score_request(sample_weight="sample_weight_alias")
    cb_alias = ScoringMonitor(scoring_val=scorer)
    X, y = make_regression(n_samples=10, n_features=2)
    X_val, y_val = make_regression(n_samples=10, n_features=2, random_state=0)
    sample_weight = np.random.randint(0, 5, size=X_val.shape[0])

    MaxIterEstimator().set_callbacks(cb_no_alias, cb_alias).fit(
        X, y, X_val=X_val, y_val=y_val, sample_weight_alias=sample_weight
    )
    alias_score_values = [
        log_item["score"] for log_item in cb_alias.get_logs().val_scores
    ]
    no_alias_score_values = [
        log_item["r2"] for log_item in cb_no_alias.get_logs().val_scores
    ]
    # Verify that requesting the sample weights result in different scores
    assert not np.allclose(alias_score_values, no_alias_score_values)


@config_context(enable_metadata_routing=True, enable_metadata_auto_requests=True)
def test_scoring_monitor_does_not_modify_scorers():
    """Check that the scorers given to the callback are not modified.

    The validation scorer automatically requests `sample_weight_val`. This must not
    leak to the user's scorer object, nor to the training scorer when both are the same
    object.
    """

    def max_sample_weight(y_true, y_pred, sample_weight):
        """Score function revealing the sample weights it received."""
        return np.max(sample_weight)

    scorer = make_scorer(max_sample_weight)
    requests_before = dict(scorer.get_metadata_routing().score.requests)

    cb = ScoringMonitor(scoring_train=scorer, scoring_val=scorer)

    X, y = make_regression(n_samples=10, n_features=2, random_state=0)
    X_val, y_val = make_regression(n_samples=10, n_features=2, random_state=1)
    MaxIterEstimator().set_callbacks(cb).fit(
        X,
        y,
        sample_weight=1,
        X_val=X_val,
        y_val=y_val,
        sample_weight_val=2,
    )

    # The original scorer object is untouched.
    assert scorer.get_metadata_routing().score.requests == requests_before

    # The training scorer used `sample_weight` and the validation scorer used
    # `sample_weight_val`, although they come from the same scorer object.
    logs = cb.get_logs()
    assert all(row["score"] == 1.0 for row in logs.train_scores)
    assert all(row["score"] == 2.0 for row in logs.val_scores)


@config_context(enable_metadata_routing=True)
def test_val_scorer_only_X():
    """Smoke test for the validation scorer on a scorer that only requires X."""

    def score_func(y_true, y_pred):
        return y_pred.sum()

    cb = ScoringMonitor(scoring_val=make_scorer(score_func))

    X, _ = make_regression(n_samples=10, n_features=2)
    X_val, _ = make_regression(n_samples=10, n_features=2)

    # without alias
    MaxIterEstimator().set_callbacks(cb).fit(X, X_val=X_val)

    # with alias
    cb.set_on_fit_task_end_request(X_val="X_val_alias")
    MaxIterEstimator().set_callbacks(cb).fit(X, X_val_alias=X_val)
