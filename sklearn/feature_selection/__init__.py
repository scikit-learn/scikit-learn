"""
The :mod:`sklearn.feature_selection` module implements feature selection
algorithms. It currently includes univariate filter selection methods and the
recursive feature elimination algorithm.
"""

from .univariate_selection import chi2
from .univariate_selection import f_classif
from .univariate_selection import f_regression
from .univariate_selection import SelectPercentile
from .univariate_selection import SelectKBest
from .univariate_selection import SelectFpr
from .univariate_selection import SelectFdr
from .univariate_selection import SelectFwe
from .univariate_selection import GenericUnivariateSelect

from .rfe import RFE
from .rfe import RFECV

import numpy as np

from ..base import TransformerMixin


class SelectorMixin(TransformerMixin):
    """"Transformer mixin selecting features based on importance weights.

    This implementation can be mixin on any estimator that exposes a
    ``feature_importances_`` or ``coef_`` attribute to evaluate the relative
    importance of individual features for feature selection.
    """
    def transform(self, X, threshold=None):
        """Reduce X to its most important features.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        threshold : string, float or None, optional (default=None)
            The threshold value to use for feature selection. Features whose
            importance is greater or equal are kept while the others are
            discarded. If "median", then the threshold value is the median of
            the feature importances. If "mean", then the threshold value is the
            mean of the feature importances. If None and if available, the
            object attribute ``threshold`` is used. Otherwise, "mean" is used
            by default.

        Returns
        -------
        X_r : array of shape [n_samples, n_selected_features]
            The input samples with only the selected features.
        """
        # Retrieve importance vector
        if hasattr(self, "feature_importances_"):
            importances = self.feature_importances_

        elif hasattr(self, "coef_"):
            if self.coef_.ndim == 1:
                importances = self.coef_ ** 2

            else:
                importances = np.sum(self.coef_ ** 2, axis=0)

        else:
            raise ValueError("Missing `feature_importances_` or `coef_`"
                             " attribute.")

        # Retrieve threshold
        if threshold is None:
            threshold = getattr(self, "threshold", "mean")

        if threshold == "median":
            threshold = np.median(importances)

        elif threshold == "mean":
            threshold = np.mean(importances)

        else:
            threshold = float(threshold)

        # Selection
        return X[:, importances >= threshold]
