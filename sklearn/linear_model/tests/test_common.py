# Author: Maria Telenczuk <https://github.com/maikia>
#
# License: BSD 3 clause

import pytest

import numpy as np

from sklearn.base import is_classifier
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import RidgeClassifierCV

from sklearn.utils import check_random_state


@pytest.mark.parametrize(
    'normalize, n_warnings, warning',
    [(True, 1, FutureWarning),
     (False, 1, FutureWarning),
     ("deprecated", 0, None)]
)
@pytest.mark.parametrize(
    "estimator",
    [LinearRegression, Ridge, RidgeCV, RidgeClassifier, RidgeClassifierCV]
)
# FIXME remove test in 1.2
def test_linear_model_normalize_deprecation_message(
     estimator,
     normalize, n_warnings, warning
):
    # check that we issue a FutureWarning when normalize was set in
    # linear model
    rng = check_random_state(0)
    n_samples = 200
    n_features = 2
    X = rng.randn(n_samples, n_features)
    X[X < 0.1] = 0.0
    y = rng.rand(n_samples)
    if is_classifier(estimator):
        y = np.sign(y)

    model = estimator(normalize=normalize)
    with pytest.warns(warning) as record:
        model.fit(X, y)
    assert len(record) == n_warnings
    if n_warnings:
        assert "'normalize' was deprecated" in str(record[0].message)
