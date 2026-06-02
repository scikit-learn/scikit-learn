# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest

from sklearn import set_config
from sklearn.base import clone
from sklearn.callback import ScoringMonitor
from sklearn.callback._scoring_monitor import ScoringMonitorLog
from sklearn.callback.tests._utils import (
    MaxIterEstimator,
    MetaEstimator,
    WhileEstimator,
)
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
    max_iter, scoring, as_pandas, X, y, root_task_id=0
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
        scores = scorer(fitted_est, X, y)
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
            scores = scorer(fitted_est, X, y)
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
):
    """Generate the expected output of a ScoringMonitor on a MetaEstimator.

    The sub-estimators are expected to be MaxIterEstimator.
    """
    meta_est_name = MetaEstimator.__name__

    expected_log = []
    for i in range(n_outer):
        for j in range(n_inner):
            estimator_log = _make_expected_output_MaxIterEstimator(
                max_iter,
                scoring,
                False,
                X,
                y,
                root_task_id=j,
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
    """Custom score to test the ScoringMonitor with a callable."""
    return 0


custom_score = make_scorer(_custom_score_func)


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
@pytest.mark.parametrize(
    "scoring_val",
    ["neg_mean_squared_error", "no_val_score"],
)
@pytest.mark.parametrize("as_pandas", [True, False])
def test_logged_values(scoring_train, scoring_val, as_pandas):
    """Test that the correct values are logged with a simple estimator."""
    if as_pandas:
        pytest.importorskip("pandas")

    if scoring_train == "no_train_score" and scoring_val == "no_val_score":
        pytest.skip("At least one scorer must be set.")

    if scoring_val != "no_val_score":
        set_config(enable_metadata_routing=True)
        X_val, y_val = make_regression(n_samples=10, n_features=2, random_state=1)
    else:
        X_val, y_val = None, None

    max_iter = 3
    callback = ScoringMonitor(scoring_train=scoring_train, scoring_val=scoring_val)
    estimator = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    if scoring_val != "no_val_score":
        estimator.fit(X=X, y=y, X_val=X_val, y_val=y_val)
    else:
        estimator.fit(X=X, y=y)

    log = callback.get_logs()
    log_train = getattr(
        log, "train_scores" if not as_pandas else "train_scores_as_pandas"
    )
    expected_log_train = _make_expected_output_MaxIterEstimator(
        max_iter, scoring_train, as_pandas, X, y
    )

    log_val = getattr(log, "val_scores" if not as_pandas else "val_scores_as_pandas")
    expected_log_val = _make_expected_output_MaxIterEstimator(
        max_iter, scoring_val, as_pandas, X_val, y_val
    )

    if as_pandas:
        assert log_train.equals(expected_log_train)
        assert log_val.equals(expected_log_val)
    else:
        assert log_train == expected_log_train
        assert log_val == expected_log_val
    set_config(enable_metadata_routing=False)


@pytest.mark.parametrize("prefer", ["processes", "threads"])
@pytest.mark.parametrize(
    "scoring_train",
    ["neg_mean_squared_error", ("neg_mean_squared_error", "r2"), custom_score],
)
@pytest.mark.parametrize("as_pandas", [True, False])
@pytest.mark.parametrize("include_lineage", [True, False])
def test_logged_values_meta_estimator_train(
    prefer, scoring_train, as_pandas, include_lineage
):
    """Test that the correct values are logged with a meta-estimator."""
    if as_pandas:
        pytest.importorskip("pandas")

    n_outer, n_inner, max_iter = 3, 2, 5
    callback = ScoringMonitor(scoring_train=scoring_train)
    est = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    meta_est = MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=2, prefer=prefer
    )
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    meta_est.fit(X=X, y=y)

    attr = "train_scores" if not as_pandas else "train_scores_as_pandas"
    log = callback.get_logs(include_lineage=include_lineage)
    log_train = getattr(log, attr)
    expected_log_train = _make_expected_output_MetaEstimator(
        n_outer,
        n_inner,
        max_iter,
        scoring_train,
        as_pandas,
        include_lineage,
        X,
        y,
    )

    if as_pandas:
        assert log_train.equals(expected_log_train)
    else:
        assert log_train == expected_log_train
    assert not log.val_scores


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


@pytest.mark.parametrize("enable_metadata_routing", [True, False])
def test_scoringmonitor_sample_weights_metadata(enable_metadata_routing):
    """Check that the ScoringMonitor works with sample weights."""
    rng = np.random.RandomState(0)
    X, y = make_regression(n_samples=100, n_features=2, random_state=rng)
    sample_weight = rng.randint(0, 5, size=X.shape[0])
    scorer_train = get_scorer("r2")
    scorer_val = "no_val_score"

    if enable_metadata_routing:
        set_config(enable_metadata_routing=True)
        X_val, y_val = make_regression(n_samples=10, n_features=2, random_state=rng)
        scorer_val = get_scorer("r2")

    # no sample weights
    callback = ScoringMonitor(scoring_train=scorer_train, scoring_val=scorer_val)
    if enable_metadata_routing:
        MaxIterEstimator().set_callbacks(callback).fit(
            X=X, y=y, X_val=X_val, y_val=y_val
        )
    else:
        MaxIterEstimator().set_callbacks(callback).fit(X=X, y=y)
    log_no_sw = callback.get_logs()

    # sample weights
    if enable_metadata_routing:
        scorer_train.set_score_request(sample_weight=True)
        scorer_val.set_score_request(sample_weight="sample_weight_val")
    callback = ScoringMonitor(scoring_train=scorer_train, scoring_val=scorer_val)
    if enable_metadata_routing:
        MaxIterEstimator().set_callbacks(callback).fit(
            X=X,
            y=y,
            X_val=X_val,
            y_val=y_val,
            sample_weight=sample_weight,
            sample_weight_val=sample_weight[:10],
        )
    else:
        MaxIterEstimator().set_callbacks(callback).fit(
            X=X, y=y, sample_weight=sample_weight
        )
    log_sw = callback.get_logs()

    assert log_no_sw.train_scores != log_sw.train_scores
    if enable_metadata_routing:
        assert log_no_sw.val_scores != log_sw.val_scores

    set_config(enable_metadata_routing=False)


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
        max_iter,
        scoring="r2",
        as_pandas=False,
        X=X,
        y=y,
    )

    assert len(log_all) == n_fits
    for log in log_all:
        assert log.train_scores == expected_log
        assert log.estimator_name == "MaxIterEstimator"


def test_scoring_monitor_listener_closed_on_gc():
    """Check that listener is closed when callback is garbage collected."""
    from sklearn.callback._transport import _listeners, _message_consumers

    callback = ScoringMonitor(scoring_train="r2")
    listener_address = callback._listener_handle.address
    assert listener_address in _listeners
    assert listener_address in _message_consumers

    del callback
    assert listener_address not in _listeners
    assert listener_address not in _message_consumers


def test_no_scorer_error():
    """Test the error when not setting a scorer."""
    with pytest.raises(
        ValueError,
        match="`scoring_train='no_train_score'` and `scoring_val='no_val_score'`",
    ):
        ScoringMonitor()


def test_scoring_val_no_metadata_routing_error():
    """Test the error when using a validation scorer with metadata routing disabled."""

    with pytest.raises(
        ValueError,
        match="a scorer on validation data .* supported if enable_metadata_routing",
    ):
        ScoringMonitor(scoring_val="r2")
