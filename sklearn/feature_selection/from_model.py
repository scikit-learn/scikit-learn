# Authors: Gilles Louppe, Mathieu Blondel
# License: BSD 3 clause

import numpy as np

from ..base import TransformerMixin
from ..externals import six
from ..utils import safe_mask, check_array
from ..utils.validation import NotFittedError, check_is_fitted


class _LearntSelectorMixin(object):
    # Note because of the extra threshold parameter in transform, this does
    # not naturally extend from SelectorMixin
    """Transformer mixin selecting features based on importance weights.

    This implementation can be mixin on any estimator that exposes a
    ``feature_importances_`` or ``coef_`` attribute to evaluate the relative
    importance of individual features for feature selection.
    """
    def transform(self, X, threshold=None):
        """Reduce X to its most important features.

        Uses ``coef_`` or ``feature_importances_`` to determine the most
        important features.  For models with a ``coef_`` for each class, the
        absolute sum over the classes is used.

        Parameters
        ----------
        X : array or scipy sparse matrix of shape (n_samples, n_features)
            The input samples.

        threshold : string or float, optional
            The threshold value to use for feature selection. Features
            whose importance is greater or equal are kept while the others
            are discarded.

            1. If the penalty used is "l1", then the threshold value
               is 1e-5 times the maximum coefficient.
            2. The default threshold used otherwise is that which is
               provided in the object, if not present than the mean is used.
            3. The threshold provided can be a float, "mean", "median"
               or multiplied by a scaling factor, such as ("1.25*mean").

        Returns
        -------
        X_r : array of shape (n_samples, n_selected_features)
            The input samples with only the selected features.
        """
        check_is_fitted(self, ('coef_', 'feature_importances_'), 
                        all_or_any=any)

        X = check_array(X, 'csc')
        # Retrieve importance vector
        if hasattr(self, "feature_importances_"):
            importances = self.feature_importances_

        elif hasattr(self, "coef_"):
            if self.coef_ is None:
                msg = "This model is not fitted yet. Please call fit() first" 
                raise NotFittedError(msg)

            if self.coef_.ndim == 1:
                importances = np.abs(self.coef_)
            else:
                importances = np.sum(np.abs(self.coef_), axis=0)

        if len(importances) != X.shape[1]:
            raise ValueError("X has different number of features than"
                             " during model fitting.")

        # Retrieve threshold
        if threshold is None:
            threshold = getattr(self, "threshold", None)

        if threshold is None:
            # Lasso has a l1 penalty but no penalty param.
            if (hasattr(self, "penalty") and self.penalty == "l1" or
                'Lasso' in self.__class__.__name__):
                # the natural default threshold is 0 when l1 penalty was used
                threshold = 1e-5 * np.max(importances)
            else:
                threshold = "mean"

        if isinstance(threshold, six.string_types):
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
        try:
            mask = importances >= threshold
        except TypeError:
            # Fails in Python 3.x when threshold is str;
            # result is array of True
            raise ValueError("Invalid threshold: all features are discarded.")

        if np.any(mask):
            mask = safe_mask(X, mask)
            return X[:, mask]
        else:
            raise ValueError("Invalid threshold: all features are discarded.")
