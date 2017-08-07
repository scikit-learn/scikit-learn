"""Stacking API"""

# Author: Caio Oliveira <caioaao@gmail.com>
# License: BSD 3 clause

from ..base import (BaseEstimator, TransformerMixin, MetaEstimatorMixin)
from ..model_selection import cross_val_predict
from ..pipeline import (_name_estimators, FeatureUnion)
from ..preprocessing import FunctionTransformer


class StackingTransformer(BaseEstimator, MetaEstimatorMixin, TransformerMixin):
    """Transformer to turn estimators into meta-estimators for model stacking

    In stacked generalization, meta estimators are combined in layers to
    improve the final result. To prevent data leaks between layers, a procedure
    similar to cross validation is adopted, where the model is trained in one
    part of the set and predicts the other part. In `StackingTransformer`, it
    happens during `fit_transform`, as the result of this procedure is what
    should be used by the next layers.

    Parameters
    ----------
    base_estimator : the estimator to be blended.

    cv : cv to be used, optional (default=3)
        Will be passed to `cross_val_predict` during `fit_transform`.

    method : string, optional (default='auto')
        Invokes the passed method name of the passed estimator. If the method
        is `auto`, will try to invoke `predict_proba` or `predict` in that
        order.

    n_jobs : int, optional (default=1)
        Number of jobs to be passed to `cross_val_predict` during
        `fit_transform`.

    """
    def __init__(self, base_estimator, cv=3, method='auto', n_jobs=1):
        self.base_estimator = base_estimator
        self.cv = cv
        self.method = method
        self.n_jobs = n_jobs

    def fit(self, X, y, **fit_params):
        """Fit the estimator.

        This should only be used in special situations. See :ref:`user guide`
        for more info.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        **fit_params : parameters to be passed to the base estimator.

        Returns
        -------
        self : object

        """
        self.base_estimator.fit(X, y, **fit_params)
        return self

    def _method_name(self):
        if self.method == 'auto':
            if getattr(self.base_estimator, 'predict_proba', None):
                method = 'predict_proba'
            else:
                method = 'predict'
        else:
            method = self.method

        return method

    def transform(self, *args, **kwargs):
        """Transform dataset.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)
            Input data to be transformed. Use ``dtype=np.float32`` for maximum
            efficiency. Sparse matrices are also supported, use sparse
            ``csr_matrix`` for maximum efficiency.

        Returns
        -------
        X_transformed : sparse matrix, shape=(n_samples, n_out)
            Transformed dataset.
        """
        t = getattr(self.base_estimator, self._method_name())
        preds = t(*args, **kwargs)

        if preds.ndim == 1:
            preds = preds.reshape(-1, 1)

        return preds

    def fit_transform(self, X, y, **fit_params):
        """Fit estimator and transform dataset.

        Parameters
        ----------
        X : array-like or sparse matrix, shape=(n_samples, n_features)
            Input data used to build forests. Use ``dtype=np.float32`` for
            maximum efficiency.

        y : array-like, shape = [n_samples]
            Target values.

        **fit_params : parameters to be passed to the base estimator.

        Returns
        -------
        X_transformed : sparse matrix, shape=(n_samples, n_out)
            Transformed dataset.
        """
        preds = cross_val_predict(self.base_estimator, X, y, cv=self.cv,
                                  method=self._method_name(),
                                  n_jobs=self.n_jobs, fit_params=fit_params)

        self.base_estimator.fit(X, y, **fit_params)

        if preds.ndim == 1:
            preds = preds.reshape(-1, 1)

        return preds


def _identity(x):
    return x


def _identity_transformer():
    """Contructs a transformer that returns its input unchanged"""
    return FunctionTransformer(_identity)


class StackLayer(FeatureUnion):
    """ Single layer for model stacking

    Parameters
    ----------
    base_estimators : list of estimators to be used in stacking

    restack : bool, optional (default=False)
        Whether input should be concatenated to the transformation.

    cv : cv to be used, optional (default=3)
        Will be passed to `StackingTransformer` for each base estimator.

    method : string, optional (default='auto')
        Invokes the passed method name of the estimators. If the method is
        `auto`, will try to invoke `predict_proba` or `predict` in that
        order.

    n_jobs : int, optional (default=1)
        Number of jobs to be passed to `cross_val_predict` during
        `fit_transform`.

    transformer_weights : dict, optional (default=None)
        Multiplicative weights for features per transformer.
        Keys are transformer names, values the weights.

    Example
    -------
    >>> from sklearn.ensemble import StackLayer
    >>> from sklearn.linear_model import Ridge
    >>> from sklearn.svm import SVC
    >>> l = StackLayer([('svc', SVC()), ('ridge', Ridge())])
    >>> l # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    StackLayer(base_estimators=[('svc', SVC(C=1.0, ...)),
                                ('ridge', Ridge(alpha=1.0, ...))],
               cv=3, method='auto', n_jobs=1, restack=False,
               transformer_weights=None)
    """
    def __init__(self, base_estimators=[], restack=False, cv=3, method='auto',
                 n_jobs=1, transformer_weights=None):
        self.base_estimators = base_estimators
        self.restack = restack
        self.cv = cv
        self.method = method
        self.n_jobs = n_jobs
        self.transformer_weights = transformer_weights
        self._update_layer()

    def _wrap_estimator(self, estimator):
        return StackingTransformer(estimator, cv=self.cv, method=self.method,
                                   n_jobs=self.n_jobs)

    def _update_layer(self):
        self.transformer_list = [(name, self._wrap_estimator(x))
                                 for name, x in self.base_estimators]
        if self.restack:
            self.transformer_list.append(('restacker',
                                          _identity_transformer()))

    def get_params(self, deep=True):
        """Get parameters for this estimator.

        Parameters
        ----------
        deep : boolean, optional
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        return self._get_params('base_estimators', deep=deep)

    def set_params(self, **kwargs):
        """Set the parameters of this estimator.

        Valid parameter keys can be listed with ``get_params()``.

        Returns
        -------
        self
        """
        self._set_params('base_estimators', **kwargs)
        self._update_layer()
        return self


def make_stack_layer(*base_estimators, **kwargs):
    """Construct a single layer for a stacked model.

    This is a shorthand for the StackLayer constructor to automatically name
    the estimators.

    Parameters
    ----------
    *base_estimators : list of base estimators.

    **kwargs : Keyword arguments to be passed to `StackLayer`.

    Returns
    -------
    StackLayer

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.ensemble import make_stack_layer
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> from sklearn.svm import SVC
    >>> t = make_stack_layer(KNeighborsClassifier(n_neighbors=3), SVC())
    >>> X = np.array([[1, 3], [.12, 1], [.5, -2], [1, -1], [-2, .1], [7, -84]])
    >>> y = np.array([1, 0, 0, 1, 0, 1])
    >>> t.fit_transform(X, y) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0.66666667,  0.33333333,  0.        ],
           [ 0.66666667,  0.33333333,  1.        ],
           [ 0.66666667,  0.33333333,  1.        ],
           [ 0.66666667,  0.33333333,  0.        ],
           [ 0.66666667,  0.33333333,  0.        ],
           [ 0.66666667,  0.33333333,  0.        ]])

    """
    return StackLayer(_name_estimators(base_estimators), **kwargs)
