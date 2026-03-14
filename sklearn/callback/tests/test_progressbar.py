# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import re
import sys

import pytest

from sklearn.base import clone
from sklearn.callback import ProgressBar
from sklearn.callback.tests._utils import (
    MaxIterEstimator,
    MetaEstimator,
    WhileEstimator,
)
from sklearn.utils._optional_dependencies import check_rich_support
from sklearn.utils.parallel import Parallel, delayed


@pytest.mark.skipif(
    sys.version_info < (3, 12, 8),
    reason="Race conditions can appear because of multiprocessing issues for python"
    " < 3.12.8.",
)
@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("prefer", ["threads", "processes"])
@pytest.mark.parametrize("InnerEstimator", [MaxIterEstimator, WhileEstimator])
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
    """Smoke test for cloning after fit with a progressbar attached.

    Initialized `ProgressBar` instances use a multiprocessing.Manager.Queue instance
    that cannot be deepcopied. This test is there to ensure that future changes
    in clone will not make it attempt to naively call copy.deepcopy on the
    _skl_callbacks attribute of the estimator.
    """
    pytest.importorskip("rich")
    est = MaxIterEstimator().set_callbacks(ProgressBar()).fit()
    clone(est)


@pytest.mark.skipif(
    sys.version_info < (3, 12, 8),
    reason="Race conditions can appear because of multiprocessing issues for python"
    " < 3.12.8.",
)
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

    if backend == "loky":
        # Since ProgressBar is pickled in different subprocesses and managers are not
        # picklable, a new manager is created for each subprocess and the queues are
        # effectively process-local.
        assert len(progressbar._run_queues) == 0
        # The monitors are process-local by construction.
        assert len(progressbar._run_monitors) == 0
    else:  # "threading"
        # The state is shared across threads so we expect one queue and monitor per fit.
        # in the shared state.
        assert len(progressbar._run_queues) == n_fits
        assert len(progressbar._run_monitors) == n_fits
        # All monitor threads are finished.
        assert not any(mon.is_alive() for mon in progressbar._run_monitors.values())
        # All queues are empty.
        assert all(queue.empty() for queue in progressbar._run_queues.values())
