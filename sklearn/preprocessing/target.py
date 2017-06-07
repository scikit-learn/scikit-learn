# Authors: Guillaume Lemaitre <guillaume.lemaitre@inria.fr>
#          Andreas Mueller <amueller@ais.uni-bonn.de>
# License: BSD 3 clause

import numpy as np

from ..base import BaseEstimator, RegressorMixin, clone, is_regressor
from ..linear_model import LinearRegression
from ..utils.validation import check_is_fitted
from ._function_transformer import FunctionTransformer

__all__ = ['TransformedTargetRegressor']


class TransformedTargetRegressor(BaseEstimator, RegressorMixin):
    """Meta-estimator to apply a transformation to the target before fitting.

    Useful for applying a non-linear transformation in regression
    problems. This transformation can be given as a Transformer such as the
    QuantileTransformer or as a function and its inverse such as ``np.log`` and
    ``np.exp``.

    The computation during ``fit`` is::

        estimator.fit(X, func(y))

    or::

        estimator.fit(X, transformer.transform(y))

    The computation during ``predict`` is::

        inverse_func(estimator.predict(X))

    or::

        transformer.inverse_transform(estimator.predict(X))

    Parameters
    ----------
    estimator : object, (default=LinearRegression())
        Estimator object derived from ``RegressorMixin``.

    transformer : object, (default=None)
        Estimator object derived from ``TransformerMixin``. Cannot be set at
        the same time as ``func`` and ``inverse_func``. If ``None`` and
        ``func`` and ``inverse_func`` are ``None`` as well, the transformer
        will be an identity transformer.

    func : function, (default=None)
        Function to apply to ``y`` before passing to ``fit``. Cannot be set at
        the same time than ``transformer``. If ``None`` and ``transformer`` is
        ``None`` as well, the function used will be the identity function.

    inverse_func : function, (default=None)
        Function apply to the prediction of the estimator. Cannot be set at
        the same time than ``transformer``. If ``None`` and ``transformer`` as
        well, the function used will be the identity function.

    check_invertible : bool, (default=True)
        Whether to check that ``transform`` followed by ``inverse_transform``
        or ``func`` followed by ``inverse_func`` lead to the original data.

    Attributes
    ----------
    estimator_ : object
        Fitted estimator.

    transformer_ : object
        Used transformer in ``fit`` and ``predict``.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.preprocessing import TransformedTargetRegressor
    >>> tt = TransformedTargetRegressor(estimator=LinearRegression(),
    ...                                 func=np.log, inverse_func=np.exp)
    >>> X = np.arange(4).reshape(-1, 1)
    >>> y = np.exp(2 * X).ravel()
    >>> tt.fit(X, y)
    ... #doctest: +NORMALIZE_WHITESPACE
    TransformedTargetRegressor(check_invertible=True,
        estimator=LinearRegression(copy_X=True, fit_intercept=True, n_jobs=1,
                                   normalize=False),
        func=<ufunc 'log'>, inverse_func=<ufunc 'exp'>, transformer=None)
    >>> tt.score(X, y)
    1.0
    >>> tt.estimator_.coef_
    array([[ 2.]])

    """
    def __init__(self, estimator=LinearRegression(), transformer=None,
                 func=None, inverse_func=None, check_invertible=True):
        self.estimator = estimator
        self.transformer = transformer
        self.func = func
        self.inverse_func = inverse_func
        self.check_invertible = check_invertible
        # we probably need to change this ones we have tags
        self._estimator_type = estimator._estimator_type

    def _validate_transformer(self, y):
        if (self.transformer is not None and
                (self.func is not None or self.inverse_func is not None)):
            raise ValueError("Both 'transformer' and functions 'func'/"
                             "'inverse_func' cannot be set at the same time.")
        elif self.transformer is not None:
            self.transformer_ = clone(self.transformer)
        else:
            self.transformer_ = FunctionTransformer(
                func=self.func, inverse_func=self.inverse_func, validate=False)
        self.transformer_.fit(y)
        if self.check_invertible:
            n_subsample = min(1000, y.shape[0])
            subsample_idx = np.random.choice(range(y.shape[0]),
                                             size=n_subsample, replace=False)
            diff = np.abs((y[subsample_idx] -
                           self.transformer_.inverse_transform(
                               self.transformer_.transform(y[subsample_idx]))))
            if np.sum(diff) > 1e-7:
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
        # memorize if y should be a multi-output
        self.y_ndim_ = y.ndim
        if y.ndim == 1:
            y_2d = y.reshape(-1, 1)
        else:
            y_2d = y
        self._validate_transformer(y_2d)
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X, self.transformer_.transform(y_2d),
                            sample_weight=sample_weight)
        return self

    def predict(self, X):
        """Predict using the base estimator, applying inverse.

        The estimator is used to predict and the ``inverse_func`` or
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
        check_is_fitted(self, "estimator_")
        pred = self.transformer_.inverse_transform(self.estimator_.predict(X))
        # if y is not a multi-output, it should be ravel
        if self.y_ndim_ == 1:
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
        disregarding the input features, would get a R^2 score of 0.0.

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

        check_is_fitted(self, "estimator_")
        if not is_regressor(self.estimator_):
            if not hasattr(self.estimator_, "_estimator_type"):
                err = "estimator has declared no _estimator_type."
            else:
                err = "estimator has _estimator_type {}".format(
                    self.estimator_._estimator_type)
            raise NotImplementedError("TransformedTargetRegressor should be a"
                                      " regressor. This " + err)
        else:
            return super(TransformedTargetRegressor, self).score(X, y,
                                                                 sample_weight)
