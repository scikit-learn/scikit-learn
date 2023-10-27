import textwrap

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

    expected_output = """\
        MetaEstimator - fit                            ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
          MetaEstimator - outer #0                     ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #0 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #1 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #2 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
          MetaEstimator - outer #1                     ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #0 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #1 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #2 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
          MetaEstimator - outer #2                     ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #0 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #1 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #2 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
          MetaEstimator - outer #3                     ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #0 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #1 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
            MetaEstimator - inner #2 | Estimator - fit ━━━━━━━━━━━━━━━━━━━━ 100% 0:00:00
        """

    assert captured.out == textwrap.dedent(expected_output)
