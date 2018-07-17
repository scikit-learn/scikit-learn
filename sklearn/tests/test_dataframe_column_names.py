# Author: Nicolas Hug
# License: BSD 3 clause

import numpy as np
import pytest

from sklearn.base import BaseEstimator
from sklearn.preprocessing import StandardScaler
from sklearn.utils.testing import SkipTest
try:
    import pandas as pd
except ImportError:
    raise SkipTest("Pandas not found")


def test_check_column_names():

    class CustomEstimator(BaseEstimator):

        def fit(self, X, y=None):
            # forgot to call _check_column_names
            return self

        def transform(self, X):

            self._check_column_names(X)
            return X

        def predict(self, X):

            self._check_column_names(X)
            return np.zeros(shape=X.shape[0])

    df = pd.DataFrame({'a': np.arange(-1, 1, .1),
                       'b': np.arange(-1, 1, .1)})

    est = CustomEstimator()
    est.fit(df)
    with pytest.raises(RuntimeError):
        est.transform(df)
    with pytest.raises(RuntimeError):
        est.predict(df)


def test_column_names_standard_scaler():
    df = pd.DataFrame({'a': np.arange(-1, 1, .1),
                       'b': np.arange(-1, 1, .1)})
    df2 = pd.DataFrame({'c': np.arange(-1, 1, .1),
                        'd': np.arange(-1, 1, .1)})

    ss = StandardScaler()

    for function in (ss.fit, ss.partial_fit):
        function(df)
        ss.transform(df)  # all is fine
        ss.inverse_transform(df)  # all is fine

        # different column order
        with pytest.raises(ValueError):
            ss.transform(df[['b', 'a']])
        with pytest.raises(ValueError):
            ss.inverse_transform(df[['b', 'a']])

        # completely different names
        with pytest.raises(ValueError):
            ss.transform(df2)
        with pytest.raises(ValueError):
            ss.inverse_transform(df2)

        # column order OK but unknown columns
        with pytest.raises(ValueError):
            ss.transform(pd.concat([df, df2]))
        with pytest.raises(ValueError):
            ss.inverse_transform(pd.concat([df, df2]))
