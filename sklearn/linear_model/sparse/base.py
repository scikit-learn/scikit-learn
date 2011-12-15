
# Author: Mathieu Blondel <mathieu@mblondel.org>
#
# License: BSD Style.

import numpy as np
from ...base import TransformerMixin


class CoefSelectTransformerMixin(TransformerMixin):
    """Mixin for linear models that can find sparse solutions.
    """

    def transform(self, X, threshold=1e-10):
        import scipy.sparse as sp
        X = sp.csc_matrix(X)
        ind = np.arange(X.shape[1])

        if len(self.coef_.shape) == 1 or self.coef_.shape[1] == 1:
            # 2-class case
            coef = np.ravel(self.coef_)
        else:
            # multi-class case
            coef = np.mean(self.coef_, axis=0)

        return X[:, ind[coef > threshold]]
