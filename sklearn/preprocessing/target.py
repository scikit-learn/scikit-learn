# Authors: Guillaume Lemaitre <guillaume.lemaitre@inria.fr>
#          Andreas Mueller <amueller@ais.uni-bonn.de>
# License: BSD 3 clause

import numpy as np

from ..base import BaseEstimator, RegressorMixin, clone, is_regressor
from ..linear_model import LinearRegression
from ..utils.fixes import signature
from ..utils.validation import check_is_fitted, check_array
from ._function_transformer import FunctionTransformer

__all__ = ['TransformedTargetRegressor']


class TransformedTargetRegressor(BaseEstimator, RegressorMixin):
    """Meta-estimator to regress on a transformed target.

    Useful for applying a non-linear transformation in regression
    problems. This transformation can be given as a Transformer such as the
    QuantileTransformer or as a function and its inverse such as ``np.log`` and
    ``np.exp``.

    The computation during ``fit`` is::

        regressor.fit(X, func(y))

    or::

        regressor.fit(X, transformer.transform(y))

    The computation during ``predict`` is::

        inverse_func(regressor.predict(X))

    or::

        transformer.inverse_transform(regressor.predict(X))

    Parameters
    ----------
    regressor : object, (default=LinearRegression())
        Regressor object such as derived from ``RegressorMixin``.

    transformer : object, (default=None)
        Estimator object such as derived from ``TransformerMixin``. Cannot be
        set at the same time as ``func`` and ``inverse_func``. If ``None`` and
        ``func`` and ``inverse_func`` are ``None`` as well, the transformer
        will be an identity transformer.

    func : function, optional
        Function to apply to ``y`` before passing to ``fit``. Cannot be set at
        the same time than ``transformer``. If ``None`` and ``transformer`` is
        ``None`` as well, the function used will be the identity function.

    inverse_func : function, optional
        Function to apply to the prediction of the regressor. Cannot be set at
        the same time than ``transformer``. If ``None`` and ``transformer`` as
        well, the function used will be the identity function. The inverse
        function is used to return to the same space of the original training
        labels during prediction.

    check_inverse : bool, (default=True)
        Whether to check that ``transform`` followed by ``inverse_transform``
        or ``func`` followed by ``inverse_func`` leads to the original data.

    Attributes
    ----------
    regressor_ : object
        Fitted regressor.

    transformer_ : object
        Used transformer in ``fit`` and ``predict``.

    y_ndim_ : int
        Number of targets.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.preprocessing import TransformedTargetRegressor
    >>> tt = TransformedTargetRegressor(regressor=LinearRegression(),
    ...                                 func=np.log, inverse_func=np.exp)
    >>> X = np.arange(4).reshape(-1, 1)
    >>> y = np.exp(2 * X).ravel()
    >>> tt.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    TransformedTargetRegressor(check_inverse=True,
        regressor=LinearRegression(copy_X=True, fit_intercept=True, n_jobs=1,
                                   normalize=False),
        func=<ufunc 'log'>, inverse_func=<ufunc 'exp'>, transformer=None)
    >>> tt.score(X, y)
    1.0
    >>> tt.regressor_.coef_
    array([[ 2.]])

    """
    def __init__(self, regressor=LinearRegression(), transformer=None,
                 func=None, inverse_func=None, check_inverse=True):
        self.regressor = regressor
        self.transformer = transformer
        self.func = func
        self.inverse_func = inverse_func
        self.check_inverse = check_inverse
        # we probably need to change this ones we have tags
        self._estimator_type = regressor._estimator_type

    def _fit_transformer(self, y, sample_weight):
        if (self.transformer is not None and
                (self.func is not None or self.inverse_func is not None)):
            raise ValueError("Both 'transformer' and functions 'func'/"
                             "'inverse_func' cannot be set at the same time.")
        elif self.transformer is not None:
            self.transformer_ = clone(self.transformer)
        else:
            self.transformer_ = FunctionTransformer(
                func=self.func, inverse_func=self.inverse_func, validate=False)
        fit_parameters = signature(self.transformer_.fit).parameters
        if "sample_weight" in fit_parameters:
            self.transformer_.fit(y, sample_weight=sample_weight)
        else:
            self.transformer_.fit(y)
        self.transformer_.fit(y)
        if self.check_inverse:
            n_subsample = min(1000, y.shape[0])
            subsample_idx = np.random.choice(range(y.shape[0]),
                                             size=n_subsample, replace=False)
            if not np.allclose(
                    y[subsample_idx],
                    self.transformer_.inverse_transform(
                        self.transformer_.transform(y[subsample_idx]))):
                raise ValueError("The provided functions or transformer are"
                                 " not strictly inverse of each other.")

    def fit(self, X, y, sample_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape (n_samples,)
            Target vector relative to X.

        sample_weight : array-like, shape (n_samples,) optional
            Array of weights that are assigned to individual samples.
            If not provided, then each sample is given unit weight.

        Returns
        -------
        self : object
            Returns self.
        """
        y = check_array(y, ensure_2d=False)
        # memorize if y should be a multi-output
        self.y_ndim_ = y.ndim
        if y.ndim == 1 and self.func is None:
            y_2d = y.reshape(-1, 1)
        else:
            y_2d = y
        self._fit_transformer(y_2d, sample_weight)
        self.regressor_ = clone(self.regressor)
        if sample_weight is not None:
            self.regressor_.fit(X, self.transformer_.transform(y_2d),
                                sample_weight=sample_weight)
        else:
            self.regressor_.fit(X, self.transformer_.transform(y_2d))
        return self

    def predict(self, X):
        """Predict using the base regressor, applying inverse.

        The regressor is used to predict and the ``inverse_func`` or
        ``inverse_transform`` is applied before returning the prediction.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            Samples.

        Returns
        -------
        y_hat : array, shape = (n_samples,)
            Predicted values.

        """
        check_is_fitted(self, "regressor_")
        pred = self.transformer_.inverse_transform(self.regressor_.predict(X))
        # if y is not a multi-output, it should be ravel
        if self.y_ndim_ == 1 and self.func is None:
            return pred.ravel()
        else:
            return pred

    def score(self, X, y, sample_weight=None):
        """Returns the coefficient of determination R^2 of the prediction.

        The coefficient R^2 is defined as (1 - u/v), where u is the regression
        sum of squares ((y_true - y_pred) ** 2).sum() and v is the residual sum
        of squares ((y_true - y_true.mean()) ** 2).sum().  Best possible score
        is 1.0 and it can be negative (because the model can be arbitrarily
        worse). A constant model that always predicts the expected value of y,
        disregarding the input features, would get a R^2 score of 0.0. Note
        that the score is computed in the original space using the
        ``inverse_transform`` of ``transformer_``.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Test samples.

        y : array-like, shape = (n_samples) or (n_samples, n_outputs)
            True values for X.

        sample_weight : array-like, shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        score : float
            R^2 of self.predict(X) wrt. y.

        """

        check_is_fitted(self, "regressor_")
        if not is_regressor(self.regressor_):
            if not hasattr(self.regressor_, "_estimator_type"):
                err = "regressor has declared no _estimator_type."
            else:
                err = "regressor has _estimator_type {}".format(
                    self.regressor_._estimator_type)
            raise NotImplementedError("TransformedTargetRegressor should be a"
                                      " regressor. This " + err)
        else:
            return super(TransformedTargetRegressor, self).score(X, y,
                                                                 sample_weight)
