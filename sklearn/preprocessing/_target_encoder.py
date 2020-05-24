import numpy as np

from ..utils.validation import check_is_fitted
from ..utils.validation import column_or_1d
from ._encoders import _BaseEncoder


class TargetRegressorEncoder(_BaseEncoder):
    """Target Encoder for Regression Targets.

    Each category is encoded based its effect on the target variable. The
    encoding scheme takes a weighted average estimated by a multilevel
    linear model.

    Read more in the :ref:`User Guide <target_regressor_encoder>`.

    Parameters
    ----------
    categories : 'auto' or a list of array-like, default='auto'
        Categories (unique values) per feature:

        - 'auto' : Determine categories automatically from the training data.
        - list : `categories[i]` holds the categories expected in the ith
          column. The passed categories should not mix strings and numeric
          values within a single feature, and should be sorted in case of
          numeric values.

        The used categories can be found in the `categories_` attribute.

    Attributes
    ----------
    cat_encodings_ : list of ndarray
        The encoding corresponding to the categories in the `categories_`
        attribute.

    categories_ : list of ndarray
        The categories of each feature determined during fitting
        (in order of the features in X and corresponding with the output
        of :meth:`transform`).

    y_mean_ : float
        The overall mean of the target.

    See Also
    --------
    sklearn.preprocessing.OrdinalEncoder : Performs an ordinal (integer)
      encoding of the categorical features.
    sklearn.preprocessing.OneHotEncoder : Performs a one-hot encoding of
      categorical features.
    """
    def __init__(self, categories='auto'):
        self.categories = categories

    def fit(self, X, y):
        """Fit the TargetRegressorEncoder to X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the categories of each feature.

        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        self : object
        """
        self._fit(X, y)
        return self

    def transform(self, X, y=None):
        """Encodes X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to encode.

        Returns
        -------
        X_out : ndarray of shape (n_samples, n_features)
            Transformed input.
        """
        check_is_fitted(self)
        X_int, X_known = self._transform(X, handle_unknown='ignore')
        return self._transform_X_int(X_int, X_known)

    def fit_transform(self, X, y):
        """Fit the encoder and encodes `X`.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to encode.

        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        X_out : ndarray of shape (n_samples, n_features)
            Transformed input.
        """
        X_int, X_known = self._fit(X, y)
        return self._transform_X_int(X_int, X_known)

    def _transform_X_int(self, X_int, X_known):
        """Transform integer encoded X. `X_int` and `X_known` are returned
        by `_transform`.
        """
        out = np.empty_like(X_int, dtype=float)
        for i, cat_encoding in enumerate(self.cat_encodings_):
            np.take(cat_encoding, X_int[:, i], out=out[:, i])
            out[~X_known[:, i], i] = self.y_mean_
        return out

    def _fit(self, X, y):
        """Fit the encoder"""
        y = column_or_1d(y, warn=True)
        fit_results = super()._fit(X, handle_unknown='ignore',
                                   return_counts=True)
        X_int, X_known = self._transform(X, handle_unknown='ignore')
        # Makes sure unknown categories are not used fot fitting
        X_int[~X_known] = -1
        self.y_mean_ = y_mean = np.mean(y)

        # y is constant the encoding will be the constant
        if np.ptp(y) == 0.0:
            self.cat_encodings_ = [
                np.full(len(cat), fill_value=y_mean, dtype=float)
                for cat in self.categories_]
            return X_int, X_known

        y_variance = np.var(y)

        n_samples, n_features = X_int.shape
        # temporary location to store masks when encoding each category
        tmp_mask = np.zeros(n_samples, dtype=bool)

        cat_encodings = []
        category_counts = fit_results['category_counts']

        for i in range(n_features):
            n_cats = len(self.categories_[i])
            cat_means = np.zeros(n_cats, dtype=float)
            cat_var_ratio = np.ones(n_cats, dtype=float)

            for encoding in range(n_cats):
                np.equal(X_int[:, i], encoding, out=tmp_mask)
                y_tmp = y[tmp_mask]
                if y_tmp.size:
                    cat_means[encoding] = np.mean(y_tmp)
                    cat_var_ratio[encoding] = np.var(y_tmp)

            # partial-pooling estimates
            cat_counts = category_counts[i]
            cat_var_ratio /= y_variance

            cat_encoded = cat_counts * cat_means + cat_var_ratio * y_mean
            cat_encoded /= cat_counts + cat_var_ratio
            cat_encodings.append(cat_encoded)

        self.cat_encodings_ = cat_encodings
        return X_int, X_known
