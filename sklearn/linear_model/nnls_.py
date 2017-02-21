import numpy as np
from scipy import optimize
from sklearn.linear_model import base
from sklearn.linear_model.base import _rescale_data
from sklearn.utils import check_X_y


class NNLS(base.LinearModel):
    def __init__(self, fit_intercept=True, normalize=False, copy_X=True):
        """
        Non negative least squares

        Parameters
        ----------
        fit_intercept : bool, optional, default : True
            whether to calculate the intercept for this model. If set
            to false, no intercept will be used in calculations
            (e.g. data is expected to be already centered).

        normalize : boolean, optional, default False
            If True, the regressors X will be normalized before regression.
            This parameter is ignored when `fit_intercept` is set to False.
            When the regressors are normalized, note that this makes the
            hyperparameters learnt more robust and almost independent of the
            number of samples. The same property is not valid for standardized
            data. However, if you wish to standardize, please use
            `preprocessing.StandardScaler` before calling `fit` on an
            estimator with `normalize=False`.

        copy_X : boolean, optional, default True
            If True, X will be copied; else, it may be overwritten.

        Attributes
        ----------
        coef_ : array, shape (n_features, ) or (n_targets, n_features)
            Estimated coefficients for the linear regression problem.
            If multiple targets are passed during the fit (y 2D), this
            is a 2D array of shape (n_targets, n_features), while if only
            one target is passed, this is a 1D array of length n_features.

        residues_ : array, shape (n_targets,) or (1,) or empty
            Sum of residuals. Squared Euclidean 2-norm for each target passed
            during the fit. If the linear regression problem is
            under-determined (the number of linearly independent rows of the
            training matrix is less than its number of linearly independent
            columns), this is an empty array. If the target vector passed
            during the fit is 1-dimensional, this is a (1,) shape array.

        """
        self.fit_intercept = fit_intercept
        self.normalize = normalize
        self.copy_X = copy_X

    def fit(self, X, y, sample_weight=None):
        """
        Fit linear model.

        Parameters
        ----------
        X : numpy array [n_samples,n_features]
            Training data

        y : numpy array of shape [n_samples]
            Target values

        sample_weight : numpy array of shape [n_samples]
            Individual weights for each sample


        Returns
        -------
        self : returns an instance of self.
        """
        X, y = check_X_y(X, y,
                         y_numeric=True, multi_output=True)

        if (sample_weight is not None) and \
           np.atleast_1d(sample_weight).ndim > 1:
            raise ValueError("Sample weights must be 1D array or scalar")

        X, y, X_offset, y_offset, X_scale = self._preprocess_data(
            X, y, fit_intercept=self.fit_intercept, normalize=self.normalize,
            copy=self.copy_X, sample_weight=sample_weight)

        if sample_weight is not None:
            # Sample weight can be implemented via a simple rescaling.
            X, y = _rescale_data(X, y, sample_weight)

        self.coef_, self.residues_ = optimize.nnls(X, y)
        self._set_intercept(X_offset, y_offset, X_scale)

        return self
