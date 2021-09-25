# Author: Lars Buitinck
# License: 3-clause BSD

import numpy as np
from ..base import BaseEstimator
from ._base import SelectorMixin
from ..utils.validation import check_is_fitted


class SelectIndex(SelectorMixin, BaseEstimator):
    """TODO"""

    def __init__(self, idx):
        self.idx = np.array(idx)

    def fit(self, X, y=None):
        """TODO"""

        X = self._validate_data(
            X,
            accept_sparse=("csr", "csc"),
            dtype=np.float64,
            force_all_finite="allow-nan",
        )

        n_columns = X.shape[1]

        support_mask = np.zeros(n_columns, dtype=bool)
        support_mask[self.idx] = True

        self.support_mask = support_mask

        return self

    def _get_support_mask(self):
        check_is_fitted(self)

        return self.support_mask

    def _more_tags(self):
        return {"allow_nan": True}
