"""Stacked ensembles"""

# Author: Caio Oliveira <caioaao@gmail.com>
# License: BSD 3 clause

from ..base import (BaseEstimator, TransformerMixin, MetaEstimatorMixin)
from ..model_selection import cross_val_predict
from ..pipeline import FeatureUnion
from ..preprocessing import FunctionTransformer


class StackingTransformer(BaseEstimator, MetaEstimatorMixin, TransformerMixin):
    """Transformer to turn estimators into meta-estimators for model stacking

    In stacked generalization, meta estimators are combined in layers to
    improve the final result. To prevent data leaks between layers, a procedure
    similar to cross validation is adopted, where the model is trained in one
    part of the set and predicts the other part. In ``StackingTransformer``, it
    happens during ``fit_transform``, as the result of this procedure is what
    should be used by the next layers. Note that this behavior is different
    from ``fit().transform()``. Read more in the :ref:`User Guide <ensemble>`.

    Parameters
    ----------
    base_estimator : the estimator to be blended.

    cv : cv to be used, optional (default=3)
        Will be passed to ``cross_val_predict`` during ``fit_transform``.

    method : string, optional (default='auto')
        Invokes the passed method name of the passed estimator. If the method
        is ``auto``, will try to invoke, for each estimator, ``predict_proba``,
        ``decision_function`` or ``predict`` in that order.

    n_jobs : int, optional (default=1)
        Number of jobs to be passed to ``cross_val_predict`` during
        ``fit_transform``.

    """
    def __init__(self, base_estimator, cv=3, method='auto', n_jobs=1):
        self.base_estimator = base_estimator
        self.cv = cv
        self.method = method
        self.n_jobs = n_jobs

    def fit(self, X, y, **fit_params):
        """Fit the estimator.

        This should only be used in special situations. Read more in the
        :ref:`User Guide <ensemble>`.

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
            elif getattr(self.base_estimator, 'decision_function', None):
                method = 'decision_function'
            else:
                method = 'predict'
        else:
            method = self.method

        return method

    def transform(self, *args, **kwargs):
        """Transform dataset.

        Note that, unlike ``fit_transform()``, this won't return the cross
        validation predictions. Read more in the :ref:`User Guide <ensemble>`.

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

        Note that this behavior is different from ``fit().transform()`` as it
        will return the cross validation predictions instead.. Read more in the
        :ref:`User Guide <ensemble>`.

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
    return FunctionTransformer(_identity, accept_sparse=True)


def make_stack_layer(base_estimators, restack=False, cv=3, method='auto',
                     n_jobs=1, n_cv_jobs=1, transformer_weights=None):
    """ Construct single layer for model stacking

    Parameters
    ----------
    base_estimators : list of estimators to be used in stacking

    restack : bool, optional (default=False)
        Whether input should be concatenated to the transformation.

    cv : cv to be used, optional (default=3)
        Will be passed to ``StackingTransformer`` for each base estimator.

    method : string, optional (default='auto')
        Invokes the passed method name of the estimators. If the method is
        ``auto``, will try to invoke ``predict_proba`` or ``predict`` in that
        order.

    n_jobs : int, optional (default=1)
        Number of jobs to run in parallel. Each job will be assigned to a base
        estimator.

    n_cv_jobs: int, optional (default=1)
        Number of jobs to be passed to each base estimator's
        ``cross_val_predict`` during ``fit_transform``.
        If ``n_jobs != 1``, ``n_cv_jobs`` must be 1 and vice-versa, since
        nested parallelism is not supported and will likely break in future
        versions.

    transformer_weights : dict, optional (default=None)
        Multiplicative weights for features per transformer.
        Keys are transformer names, values the weights.

    Returns
    -------
    FeatureUnion

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.ensemble import make_stack_layer
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> from sklearn.svm import SVC
    >>> t = make_stack_layer([('knn', KNeighborsClassifier(n_neighbors=3)),
    ...                       ('svc', SVC())])
    >>> X = np.array([[1, 3], [.12, 1], [.5, -2], [1, -1], [-2, .1], [7, -84]])
    >>> y = np.array([1, 0, 0, 1, 0, 1])
    >>> t.fit_transform(X, y) # doctest: +NORMALIZE_WHITESPACE
    array([[  6.66666667e-01,   3.33333333e-01,  -4.44715665e-04],
           [  6.66666667e-01,   3.33333333e-01,   1.04443234e-02],
           [  6.66666667e-01,   3.33333333e-01,   2.01594603e-02],
           [  6.66666667e-01,   3.33333333e-01,  -6.22390628e-02],
           [  6.66666667e-01,   3.33333333e-01,  -1.10049720e-01],
           [  6.66666667e-01,   3.33333333e-01,  -4.09417930e-02]])
    """
    transformer_list = [(name, StackingTransformer(estimator, cv=cv,
                                                   method=method,
                                                   n_jobs=n_cv_jobs))
                        for name, estimator in base_estimators]
    if restack:
        transformer_list.append(('restacker', _identity_transformer()))

    return FeatureUnion(transformer_list, n_jobs=n_jobs,
                        transformer_weights=transformer_weights)
