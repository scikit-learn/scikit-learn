"""Regression test for thread-safety issues in liblinear-backed estimators."""

import numpy as np
import pytest

from sklearn.datasets import make_regression
from sklearn.exceptions import ConvergenceWarning
from sklearn.svm import LinearSVR
from sklearn.utils.parallel import Parallel, delayed


@pytest.mark.xfail(
    reason=(
        "LinearSVR/liblinear is not thread-safe when fitting in threads; see "
        "https://github.com/scikit-learn/scikit-learn/issues/31883"
    ),
    strict=False,
)
def test_linear_svr_fit_threading_backend_is_deterministic():
    # Use a fixed dataset and random_state to make sequential fits deterministic.
    X, y = make_regression(
        n_samples=200,
        n_features=30,
        noise=0.1,
        random_state=42,
    )

    C_range = np.logspace(-6, 6, 13)

    def fit_one(C):
        est = LinearSVR(C=C, random_state=0, max_iter=5000)
        with pytest.warns(ConvergenceWarning):
            est.fit(X, y)
        return est.coef_.copy()

    sequential = [fit_one(C) for C in C_range]

    # Repeat a few times to increase the chance of triggering the race.
    for _ in range(3):
        parallel = Parallel(n_jobs=4, backend="threading")(
            delayed(fit_one)(C) for C in C_range
        )
        np.testing.assert_array_equal(sequential, parallel)
