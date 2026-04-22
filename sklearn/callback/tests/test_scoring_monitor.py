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


def _make_expected_output_MaxIterEstimator(
    max_iter, scoring, as_frame, X, y, root_task_id=0
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

    # fit root task
    fitted_est = MaxIterEstimator(max_iter=max_iter).fit()
    log_item = {
        "parent_task_id_path": (),
        "estimator_name": est_name,
        "task_name": "fit",
        "task_id": root_task_id,
        "sequential_subtasks": True,
    }
    scores = scorer(fitted_est, X, y)
    if not isinstance(scores, dict):
        scores = {score_names[0]: scores}
    expected_log.append({**log_item, **scores})

    # fit loop iterations
    for i in range(max_iter):
        fitted_est = MaxIterEstimator(max_iter=i + 1).fit()
        log_item = {
            "parent_task_id_path": (root_task_id,),
            "estimator_name": est_name,
            "task_name": f"iteration {i}",
            "task_id": i,
            "sequential_subtasks": True,
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
    n_outer, n_inner, max_iter, scoring, as_frame, include_lineage, X, y
):
    """Generate the expected output of a ScoringMonitor on a MetaEstimator.

    The sub-estimators are expected to be MaxIterEstimator.
    """
    meta_est_name = MetaEstimator.__name__

    expected_log = []
    for i in range(n_outer):
        for j in range(n_inner):
            estimator_log = _make_expected_output_MaxIterEstimator(
                max_iter, scoring, False, X, y, root_task_id=j
            )
            for row in estimator_log:
                row["parent_task_id_path"] = (0, i) + row["parent_task_id_path"]
                expected_log.append(row)

    if include_lineage:
        # fit root task
        expected_log.append(
            {
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
                    "parent_task_id_path": (0,),
                    "estimator_name": meta_est_name,
                    "task_name": "outer",
                    "task_id": i,
                    "sequential_subtasks": True,
                }
            )
        # add task_id_path column
        for i, row in enumerate(expected_log):
            task_id_path = row["parent_task_id_path"] + (row["task_id"],)
            expected_log[i] = {"task_id_path": task_id_path, **row}

    sorting_key = lambda x: (len(x["parent_task_id_path"]), x["parent_task_id_path"])
    expected_log.sort(key=sorting_key)

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
    log_iterations = [row for row in log if row["parent_task_id_path"] == (0,)]

    scorer = check_scoring(None, "r2")
    for i in range(max_iter):
        estimator = MaxIterEstimator(max_iter=i + 1).fit(X, y)
        score = scorer(estimator, X, y)
        assert_allclose(score, log_iterations[i]["r2"])


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
@pytest.mark.parametrize("include_lineage", [True, False])
def test_logged_values_meta_estimator(prefer, scoring, as_frame, include_lineage):
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

    log = callback.get_logs(
        as_frame=as_frame, select="most_recent", include_lineage=include_lineage
    )["data"]
    expected_log = _make_expected_output_MetaEstimator(
        n_outer, n_inner, max_iter, scoring, as_frame, include_lineage, X, y
    )

    if as_frame:
        assert log.equals(expected_log)
    else:
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


@pytest.mark.parametrize("as_frame", [True, False])
def test_get_logs_include_lineage_ancestor_retrieval(as_frame):
    """Check that ancestor task info can be retrieved from a given task."""
    if as_frame:
        pytest.importorskip("pandas")

    n_outer, n_inner, max_iter = 2, 3, 5
    callback = ScoringMonitor(scoring="r2")
    est = MaxIterEstimator(max_iter=max_iter).set_callbacks(callback)
    meta_est = MetaEstimator(est, n_outer=n_outer, n_inner=n_inner)
    X, y = make_regression(n_samples=100, n_features=2, random_state=0)

    meta_est.fit(X=X, y=y)

    log = callback.get_logs(as_frame=as_frame, include_lineage=True)["data"]

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

    if as_frame:
        sibling_rows = log.loc[log["parent_task_id_path"] == parent_path]
        ancestors_info_path = _get_ancestors_info_path_pandas(log, parent_path)
    else:
        sibling_rows = [row for row in log if row["parent_task_id_path"] == parent_path]
        ancestors_info_path = _get_ancestors_info_path(log, parent_path)

    assert len(sibling_rows) == max_iter
    assert ancestors_info_path == expected_ancestors_info_path
