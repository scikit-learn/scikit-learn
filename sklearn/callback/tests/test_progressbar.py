# License: BSD 3 clause
# Authors: the scikit-learn developers

import re

import pytest

from sklearn.callback import ProgressBar
from sklearn.utils._optional_dependencies import check_rich_support
from sklearn.utils._testing import SkipTest

from ._utils import Estimator, MetaEstimator, WhileEstimator


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("prefer", ["threads", "processes"])
@pytest.mark.parametrize("InnerEstimator", [Estimator, WhileEstimator])
@pytest.mark.parametrize("max_estimator_depth", [1, 2, None])
def test_progressbar(n_jobs, prefer, InnerEstimator, max_estimator_depth, capsys):
    """Check the output of the progress bars and their completion."""
    pytest.importorskip("rich")

    n_inner = 2
    n_outer = 3

    est = InnerEstimator()
    meta_est = MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=n_jobs, prefer=prefer
    )
    meta_est.set_callbacks(ProgressBar(max_estimator_depth=max_estimator_depth))
    meta_est.fit()

    captured = capsys.readouterr()

    assert re.search(r"MetaEstimator - fit", captured.out)
    for i in range(n_outer):
        assert re.search(rf"MetaEstimator - outer #{i}", captured.out)

    # Progress bars of inner estimators are displayed only if max_estimator_depth > 1
    # (or None, which means all levels are displayed)
    if max_estimator_depth is None or max_estimator_depth > 1:
        for i in range(n_inner):
            assert re.search(
                rf"MetaEstimator - inner #{i} | {est.__class__.__name__} - fit",
                captured.out,
            )

    # Check that all bars are 100% complete
    assert re.search(r"100%", captured.out)
    assert not re.search(r"[1-9]%", captured.out)


def test_progressbar_requires_rich_error():
    """Check that we raise an informative error when rich is not installed."""
    try:
        check_rich_support("test_progressbar_requires_rich_error")
    except ImportError:
        err_msg = "Progressbar requires rich"
        with pytest.raises(ImportError, match=err_msg):
            ProgressBar()
    else:
        raise SkipTest("This test requires rich to not be installed.")
