import pytest

import numpy as np

from sklearn.impute._base import _BaseImputer


@pytest.fixture
def data():
    X = np.random.randn(10, 2)
    X[::2] = np.nan
    return X


class NoFitIndicatorImputer(_BaseImputer):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return self._concatenate_indicator(X, self._transform_indicator(X))


class NoTransformIndicatorImputer(_BaseImputer):
    def fit(self, X, y=None):
        super()._fit_indicator(X)
        return self

    def transform(self, X, y=None):
        return self._concatenate_indicator(X, None)


def test_base_imputer_not_fit(data):
    imputer = NoFitIndicatorImputer(add_indicator=True)
    err_msg = "Make sure to call _fit_indicator before _transform_indicator"
    with pytest.raises(ValueError, match=err_msg):
        imputer.fit(data).transform(data)
    with pytest.raises(ValueError, match=err_msg):
        imputer.fit_transform(data)


def test_base_imputer_not_transform(data):
    imputer = NoTransformIndicatorImputer(add_indicator=True)
    err_msg = ("Call _fit_indicator and _transform_indicator in the "
               "imputer implementation")
    with pytest.raises(ValueError, match=err_msg):
        imputer.fit(data).transform(data)
    with pytest.raises(ValueError, match=err_msg):
        imputer.fit_transform(data)
