# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import re
import textwrap
from unittest import mock

import pytest

from sklearn.base import clone
from sklearn.callback import ProgressBar
from sklearn.callback.tests._utils import (
    HeterogeneousMetaEstimator,
    MaxIterEstimator,
    MetaEstimator,
    NoSubtaskEstimator,
    WhileEstimator,
)
from sklearn.utils._optional_dependencies import check_rich_support
from sklearn.utils._testing import (
    assert_allclose,
    assert_run_python_script_without_output,
)
from sklearn.utils.parallel import Parallel, delayed


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("prefer", ["threads", "processes"])
@pytest.mark.parametrize(
    "InnerEstimator", [MaxIterEstimator, WhileEstimator, NoSubtaskEstimator]
)
@pytest.mark.parametrize("max_propagation_depth", [1, 2, None])
def test_progressbar(n_jobs, prefer, InnerEstimator, max_propagation_depth, capsys):
    """Check the output of the progress bars and their completion."""
    pytest.importorskip("rich")

    n_inner = 2
    n_outer = 3

    est = InnerEstimator()
    meta_est = MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=n_jobs, prefer=prefer
    )
    meta_est.set_callbacks(ProgressBar(max_propagation_depth=max_propagation_depth))
    meta_est.fit()

    captured = capsys.readouterr()

    assert re.search(r"MetaEstimator - fit", captured.out)
    for i in range(n_outer):
        assert re.search(rf"MetaEstimator - outer #{i}", captured.out)

    # Progress bars of inner estimators are displayed only if max_propagation_depth > 1
    # (or None, which means all levels are displayed)
    if max_propagation_depth is None or max_propagation_depth > 1:
        for i in range(n_inner):
            assert re.search(
                rf"MetaEstimator - inner \| {est.__class__.__name__} - fit #{i}",
                captured.out,
            )

    # Check that all bars are 100% complete
    assert re.search(r"100%", captured.out)
    assert not re.search(r"\b[0-9]{1,2}%", captured.out)


def test_progressbar_requires_rich_error():
    """Check that we raise an informative error when rich is not installed."""
    try:
        check_rich_support("test_progressbar_requires_rich_error")
        pytest.skip("This test requires rich to not be installed.")
    except ImportError:
        err_msg = "Progressbar requires rich"
        with pytest.raises(ImportError, match=err_msg):
            ProgressBar()


def test_clone_after_fit():
    """Cloning an estimator with a ProgressBar preserves the callback reference.

    clone copies _skl_callbacks by reference so that a single callback instance can
    track every clone.
    """
    pytest.importorskip("rich")
    pb = ProgressBar()
    est = MaxIterEstimator().set_callbacks(pb).fit()
    cloned = clone(est)
    assert cloned._skl_callbacks[0] is pb


@pytest.mark.parametrize("backend", ["threading", "loky"])
def test_progressbar_no_callback_support(backend):
    """Sanity check for ProgressBar within function not supporting callbacks.

    It's hard to check the output from sub-processes so this test only checks that it
    doesn't crash and that there are no threads leftover running.
    """
    pytest.importorskip("rich")

    def clone_and_fit(estimator):
        clone(estimator).fit()

    def func(estimator, n_fits):
        Parallel(n_jobs=2, backend=backend)(
            delayed(clone_and_fit)(estimator) for _ in range(n_fits)
        )

    progressbar = ProgressBar()
    n_fits = 4
    func(MaxIterEstimator().set_callbacks(progressbar), n_fits=n_fits)

    from sklearn.callback._progressbar import _run_monitors, _run_queues

    # After fit completes, `teardown` has popped the per-fit queue and monitor from
    # the module-level registry and the per-fit listener handle from the instance.
    assert progressbar._listener_handles == {}
    assert _run_queues == {}
    assert _run_monitors == {}


@pytest.mark.parametrize("prefer", ["threads", "processes"])
def test_progressbar_outside_main_module(prefer):
    """Check that ProgressBar does not trigger spawn errors outside `__main__`."""
    pytest.importorskip("rich")
    code = f"""
    from sklearn.callback import ProgressBar
    from sklearn.callback.tests._utils import MaxIterEstimator, MetaEstimator

    est = MaxIterEstimator()
    meta_est = MetaEstimator(
        est, n_outer=2, n_inner=1, n_jobs=2, prefer='{prefer}'
    )
    meta_est.set_callbacks(ProgressBar())
    meta_est.fit()
    """

    pattern = "An attempt has been made to start a new process"
    assert_run_python_script_without_output(
        textwrap.dedent(code), pattern=pattern, timeout=120
    )


def test_progress_during_fit():
    """Check that the completion of a bottom-level progressbar increments linearly."""
    pytest.importorskip("rich")
    from sklearn.callback._progressbar import RichProgressMonitor

    records = []
    orig_on_task_end = RichProgressMonitor._on_task_end

    def recording_on_task_end(self, task_info):
        orig_on_task_end(self, task_info)
        records.append(self.root_rich_task.progress)

    max_iter = 7
    with mock.patch.object(RichProgressMonitor, "_on_task_end", recording_on_task_end):
        MaxIterEstimator(max_iter=max_iter).set_callbacks(ProgressBar()).fit()

    # progress after each iteration + 100% at the end of fit
    expected = [i / max_iter for i in range(1, max_iter + 1)] + [1.0]
    assert_allclose(records, expected)


@pytest.mark.parametrize(
    "meta_estimator",
    [
        MetaEstimator(MaxIterEstimator(max_iter=5), n_outer=4),
        HeterogeneousMetaEstimator([MaxIterEstimator(max_iter=5), None, None, None]),
        HeterogeneousMetaEstimator([None, None, MaxIterEstimator(max_iter=5), None]),
        HeterogeneousMetaEstimator([None, None, None, MaxIterEstimator(max_iter=5)]),
    ],
)
def test_progress_during_fit_composition(meta_estimator):
    """Check the recursive computation of the progress of nested estimators."""
    pytest.importorskip("rich")
    from sklearn.callback._progressbar import RichProgressMonitor

    records = []
    orig_on_task_end = RichProgressMonitor._on_task_end

    def check_progress(task):
        # Check recursively that the completion of each progress bar is the average
        # completion of its children.
        if not task.children:
            expected = 1.0
        else:
            expected = sum(c.progress for c in task.children.values()) / task.total
            for child in task.children.values():
                check_progress(child)
        assert_allclose(task.progress, expected)

    def recording_on_task_end(self, task_info):
        path = task_info["path"]
        orig_on_task_end(self, task_info)
        records.append([path, self.root_rich_task.progress])
        check_progress(self.root_rich_task)

    with mock.patch.object(RichProgressMonitor, "_on_task_end", recording_on_task_end):
        meta = meta_estimator.set_callbacks(ProgressBar())
        meta.fit()

    # MetaEstimator's fit has 4 subtasks (n_outer=4). The top-level progress bar should
    # be at 25%, 50%, 75% and 100% at the end of the subtasks.
    end_of_outer_subtasks = [progress for path, progress in records if len(path) == 2]
    expected_progress = [0.25, 0.5, 0.75, 1.0]
    assert_allclose(end_of_outer_subtasks, expected_progress)


def test_estimator_without_subtasks(capsys):
    """Check that a progress bar is displayed for an estimator without subtasks."""
    pytest.importorskip("rich")

    NoSubtaskEstimator().set_callbacks(ProgressBar()).fit()
    captured = capsys.readouterr()
    assert re.search(r"NoSubtaskEstimator - fit", captured.out)
    assert re.search(r"100%", captured.out)
