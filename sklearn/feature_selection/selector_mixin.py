# Authors: Gilles Louppe, Mathieu Blondel
# License: BSD

import numpy as np

from ..base import TransformerMixin
from ..utils import safe_mask, atleast2d_or_csr


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
        X = atleast2d_or_csr(X)
        # Retrieve importance vector
        if hasattr(self, "feature_importances_"):
            importances = self.feature_importances_
            if importances is None:
                raise ValueError("Importance weights not computed. Please  "
                    "Set the compute_importances parameter before fit.")

        elif hasattr(self, "coef_"):
            if self.coef_.ndim == 1:
                importances = np.abs(self.coef_)

            else:
                importances = np.sum(np.abs(self.coef_), axis=0)

        else:
            raise ValueError("Missing `feature_importances_` or `coef_`"
                             " attribute, did you forget to set the "
                             "estimator's parameter to compute it?")
        if len(importances) != X.shape[1]:
            raise ValueError("X has different number of features than"
                    " during model fitting.")

        # Retrieve threshold
        if threshold is None:
            if hasattr(self, "penalty") and self.penalty == "l1":
                # the natural default threshold is 0 when l1 penalty was used
                threshold = getattr(self, "threshold", 1e-5)
            else:
                threshold = getattr(self, "threshold", "mean")

        if isinstance(threshold, basestring):
            if "*" in threshold:
                scale, reference = threshold.split("*")
                scale = float(scale.strip())
                reference = reference.strip()

                if reference == "median":
                    reference = np.median(importances)
                elif reference == "mean":
                    reference = np.mean(importances)
                else:
                    raise ValueError("Unknown reference: " + reference)

                threshold = scale * reference

            elif threshold == "median":
                threshold = np.median(importances)

            elif threshold == "mean":
                threshold = np.mean(importances)

        else:
            threshold = float(threshold)

        # Selection
        mask = importances >= threshold

        if np.any(mask):
            mask = safe_mask(X, mask)
            return X[:, mask]
        else:
            raise ValueError("Invalid threshold: all features are discarded.")
