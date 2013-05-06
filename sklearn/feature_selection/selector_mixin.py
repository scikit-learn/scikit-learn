# Authors: Gilles Louppe, Mathieu Blondel
# License: BSD 3 clause

import numpy as np

from ..base import TransformerMixin
<<<<<<< HEAD
from ..externals import six
from ..utils import safe_mask, atleast2d_or_csc

=======
from .etc import SelectBetween
>>>>>>> Reimplement SelectorMixin using SelectBetween

class SelectorMixin(TransformerMixin):
    """Transformer mixin selecting features based on importance weights.

    This implementation can be mixin on any estimator that exposes a
    ``feature_importances_`` or ``coef_`` attribute to evaluate the relative
    importance of individual features for feature selection.
    """

    def make_transformer(self):
        if getattr(self, 'penalty') == 'l1':
            # the natural default threshold is 0 when l1 penalty was used
            default_threshold = 1e-5
        else:
            default_threshold = 'mean'
        return SelectBetween(self.__get_feature_importances,
                             minimum=default_threshold)

    def __get_feature_importances(self, X=None, y=None):
        if hasattr(self, "feature_importances_"):
            importances = self.feature_importances_
            if importances is None:
                raise ValueError("Importance weights not computed. Please set"
                                 " the compute_importances parameter before "
                                 "fit.")

        elif hasattr(self, "coef_"):
            if self.coef_.ndim == 1:
                importances = np.abs(self.coef_)

            else:
                importances = np.sum(np.abs(self.coef_), axis=0)

        else:
            raise ValueError("Missing `feature_importances_` or `coef_`"
                             " attribute, did you forget to set the "
                             "estimator's parameter to compute it?")
        return importances

    def transform(self, X, threshold=None):
        """Reduce X to its most important features.

        Parameters
        ----------
        X : array or scipy sparse matrix of shape [n_samples, n_features]
            The input samples.

        threshold : string, float or None, optional (default=None)
            The threshold value to use for feature selection. Features whose
            importance is greater or equal are kept while the others are
            discarded. If "median" (resp. "mean"), then the threshold value is
            the median (resp. the mean) of the feature importances. A scaling
            factor (e.g., "1.25*mean") may also be used. If None and if
            available, the object attribute ``threshold`` is used. Otherwise,
            "mean" is used by default.

        Returns
        -------
        X_r : array of shape [n_samples, n_selected_features]
            The input samples with only the selected features.
        """
        transformer = self.make_transformer()
        if threshold is not None:
            transformer.minimum = threshold
        return transformer.fit_transform(X)
