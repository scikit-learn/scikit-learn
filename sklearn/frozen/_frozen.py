# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ..base import BaseEstimator
from ..utils._available_if import available_if


def _estimator_has(attr):
    """Check that final_estimator has `attr`.

    Used together with `available_if`.
    """

    def check(self):
        # raise original `AttributeError` if `attr` does not exist
        getattr(self.estimator, attr)
        return True

    return check


class FrozenEstimator(BaseEstimator):
    def __init__(self, estimator):
        self.estimator = estimator

    @available_if(_estimator_has("__getitem__"))
    def __getitem__(self, *args, **kwargs):
        return self.estimator.__getitem__(*args, **kwargs)

    def __getattr__(self, name):
        # `estimator`'s attributes are now accessible
        return getattr(self.estimator, name)

    def __sklearn_clone__(self):
        return self

    def fit(self, X, y, *args, **kwargs):
        # Fitting does not change the state of the estimator
        return self

    @available_if(_estimator_has("fit_transform"))
    def fit_transform(self, X, y=None, *args, **kwargs):
        # fit_transform only transforms the data
        return self.estimator.transform(X, y, *args, **kwargs)
