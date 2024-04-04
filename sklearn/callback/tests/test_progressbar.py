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
def test_progressbar(n_jobs, prefer, InnerEstimator, capsys):
    """Check the output of the progress bars and their completion."""
    pytest.importorskip("rich")

    est = InnerEstimator()
    meta_est = MetaEstimator(est, n_jobs=n_jobs, prefer=prefer)
    meta_est._set_callbacks(ProgressBar())
    meta_est.fit()

    captured = capsys.readouterr()

    assert re.search(r"MetaEstimator - fit", captured.out)
    for i in range(4):
        assert re.search(rf"MetaEstimator - outer #{i}", captured.out)
    for i in range(3):
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
