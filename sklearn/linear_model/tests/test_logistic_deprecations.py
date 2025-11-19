import pytest

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV


def test_logistic_n_jobs_deprecated_init():
    with pytest.warns(FutureWarning, match="n_jobs"):
        LogisticRegression(n_jobs=2)


def test_logisticcv_n_jobs_deprecated_init():
    with pytest.warns(FutureWarning, match="n_jobs"):
        LogisticRegressionCV(n_jobs=2)
