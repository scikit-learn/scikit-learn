# Authors: Andreas Mueller < andreas.mueller@columbia.edu>
#          Guillaume Lemaitre <guillaume.lemaitre@inria.fr>
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
        Regressor object such as derived from ``RegressorMixin``. This
        regressor will be cloned during fitting.

    transformer : object, (default=None)
        Estimator object such as derived from ``TransformerMixin``. Cannot be
        set at the same time as ``func`` and ``inverse_func``. If ``None`` and
        ``func`` and ``inverse_func`` are ``None`` as well, the transformer
        will be an identity transformer. The transformer will be cloned during
        fitting.

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
                               func=<ufunc 'log'>,
                               inverse_func=<ufunc 'exp'>,
                               regressor=LinearRegression(copy_X=True,
                                                          fit_intercept=True,
                                                          n_jobs=1,
                                                          normalize=False),
                               transformer=None)
    >>> tt.score(X, y)
    1.0
    >>> tt.regressor_.coef_
    array([ 2.])

    """
    def __init__(self, regressor=None, transformer=None,
                 func=None, inverse_func=None, check_inverse=True):
        self.regressor = regressor
        self.transformer = transformer
        self.func = func
        self.inverse_func = inverse_func
        self.check_inverse = check_inverse

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
        if self.check_inverse:
            n_subsample = min(10, y.shape[0])
            subsample_idx = np.random.choice(range(y.shape[0]),
                                             size=n_subsample, replace=False)
            if not np.allclose(
                    y[subsample_idx],
                    self.transformer_.inverse_transform(
                        self.transformer_.transform(y[subsample_idx])),
                    atol=1e-4):
                raise ValueError("The provided functions or transformer are"
                                 " not strictly inverse of each other. If"
                                 " you are sure you want to proceed regardless"
                                 ", set 'check_inverse=False'")

    def fit(self, X, y, sample_weight=None):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape (n_samples,)
            Target values.

        sample_weight : array-like, shape (n_samples,) optional
            Array of weights that are assigned to individual samples.
            If not provided, then each sample is given unit weight.

        Returns
        -------
        self : object
            Returns self.
        """
        y = check_array(y, ensure_2d=False)
        self.y_ndim_ = y.ndim
        if y.ndim == 1 and self.func is None:
            y_2d = y.reshape(-1, 1)
        else:
            y_2d = y
        self._fit_transformer(y_2d, sample_weight)
        if self.regressor is None:
            self.regressor_ = LinearRegression()
        else:
            self.regressor_ = clone(self.regressor)
        if sample_weight is not None:
            self.regressor_.fit(X, self.transformer_.fit_transform(y_2d),
                                sample_weight=sample_weight)
        else:
            self.regressor_.fit(X, self.transformer_.fit_transform(y_2d))
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
        if self.y_ndim_ == 1 and self.func is None:
            return pred.ravel()
        else:
            return pred
