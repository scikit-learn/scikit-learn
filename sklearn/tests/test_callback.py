# License: BSD 3 clause

import numpy as np
import pytest

from sklearn.datasets import load_iris
from sklearn._callbacks import Callback
from sklearn.utils.testing import all_estimators
from sklearn.utils.estimator_checks import _construct_instance, _safe_tags


@pytest.fixture(scope="module")
def iris():
    X, y = load_iris(return_X_y=True)
    return X, y.reshape(-1, 1)


@pytest.mark.parametrize(
    "name, Estimator",
    [
        (name, Estimator)
        for (name, Estimator) in all_estimators()
        if not name.startswith("_")
    ],
)
def test_callback(name, Estimator, iris):

    if name in [
        "SparseRandomProjection",
        "SelectKBest",
        "MultiLabelBinarizer",
        "LabelEncoder",
        "LabelBinarizer",
        "KernelCenterer",
        "IsotonicRegression",
        "GaussianRandomProjection",
        "DictVectorizer",
    ]:
        # TODO: make estimator.fit(X, y) work for these with
        # some appropriate data / parameters
        pytest.skip("Callbacks not yet supported")

    if name in [
        "DictionaryLearning",
        "MiniBatchDictionaryLearning",
        "TheilSenRegressor",
    ]:
        # TODO: make these tests faster
        pytest.skip("Slow tests, skipping for now")

    estimator = _construct_instance(Estimator)

    tags = _safe_tags(estimator)

    class CheckCallback(Callback):
        def __init__(self):
            self.n_calls = 0

        def fit(self, X, y):
            pass

        def __call__(self, **kwargs):
            self.n_calls += 1

    callback = CheckCallback()
    estimator._set_callbacks([callback])
    if tags.get("X_types", []) == ["string"]:
        X = ["some document", "another document"]
        y = None
    else:
        X, y = iris
    assert callback.n_calls == 0
    estimator.fit(X, y)
    assert callback.n_calls > 0
