# License: BSD 3 clause
# Authors: the scikit-learn developers

import re

import pytest

from sklearn.callback import ProgressBar

from ._utils import Estimator, MetaEstimator


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("prefer", ["threads", "processes"])
def test_progressbar(n_jobs, prefer, capsys):
    """Check the output of the progress bars and their completion"""
    pytest.importorskip("rich")

    est = Estimator()
    meta_est = MetaEstimator(est, n_jobs=n_jobs, prefer=prefer)
    meta_est._set_callbacks(ProgressBar())
    meta_est.fit(None, None)

    captured = capsys.readouterr()

    assert re.search(r"MetaEstimator - fit", captured.out)
    for i in range(4):
        assert re.search(rf"MetaEstimator - outer #{i}", captured.out)
    for i in range(3):
        assert re.search(rf"MetaEstimator - inner #{i} | Estimator - fit", captured.out)

    # Check that all bars are 100% complete
    assert re.search(r"100%", captured.out)
    assert not re.search(r"[1-9]%", captured.out)
