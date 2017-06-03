"""Stacking API"""

# Author: Caio Oliveira <caioaao@gmail.com>
# License: BSD 3 clause

from ..base import (BaseEstimator, TransformerMixin, MetaEstimatorMixin)
from ..model_selection import cross_val_predict
from ..pipeline import (_name_estimators, FeatureUnion, Pipeline)
from ..preprocessing import FunctionTransformer


class StackMetaEstimator(BaseEstimator, MetaEstimatorMixin, TransformerMixin):
    """Transformer to turn estimators into meta-estimators for model stacking

    In stacked generalization, meta estimators are combined in layers to
    improve the final result. To prevent data leaks between layers, a procedure
    similar to cross validation is adopted, where the model is trained in one
    part of the set and predicts the other part. In `StackMetaEstimator`, it
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

    def fit(self, *args, **kwargs):
        raise NotImplementedError("Use `fit_transform` instead")

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
        t = getattr(self.base_estimator, self._method_name())
        preds = t(*args, **kwargs)

        if preds.ndim == 1:
            preds = preds.reshape(-1, 1)

        return preds

    def fit_transform(self, X, y, **fit_params):
        preds = cross_val_predict(self.base_estimator, X, y, cv=self.cv,
                                  method=self._method_name(),
                                  n_jobs=self.n_jobs, fit_params=fit_params)

        self.base_estimator.fit(X, y, **fit_params)

        if preds.ndim == 1:
            preds = preds.reshape(-1, 1)

        return preds


def _identity_transformer():
    """Contructs a transformer that does nothing"""
    return FunctionTransformer(lambda x: x)


def make_stack_layer(*base_estimators, **kwargs):
    """Construct a single layer for a stacked model.

    This is a wrapper around pipelines to provide a more convenient API for
    stacking models.

    Parameters
    ----------
    *base_estimators : list of base estimators.

    restacking : bool, optional (default=False)
        When true, the transformer will return the input.

    transformer_weights : dict, optional (default=None)
        Multiplicative weights for features per transformer.
        Keys are transformer names, values the weights.

    **kwargs : Keyword arguments to be passed to `StackMetaEstimator`.

    Returns
    -------
    f : FeatureUnion with every base estimator wrapped in a
         `StackMetaEstimator`.

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
    restacking = kwargs.pop('restacking', False)
    transformer_weights = kwargs.pop('transformer_weights', None)

    named_estimators = [(n, StackMetaEstimator(estimator, **kwargs))
                        for n, estimator in _name_estimators(base_estimators)]
    if restacking:
        named_estimators.extend(_name_estimators([_identity_transformer()]))

    return FeatureUnion(named_estimators,
                        transformer_weights=transformer_weights)

def stack_estimators(*estimators, **kwargs):
    """Construct a stacked estimator

    This is a wrapper around pipelines to provide a more convenient API for
    stacking models. Estimators in `estimators_matrix` are wrapped in
    `StackMetaEstimator`.

    Parameters
    ----------
    *estimators : Estimators for stacking. Every param but the last one must
    be a list of estimators that will be wrapped with `StackMetaEstimator`.
    The last argument on the list of estimators should be a single estimator
    and will be used as the last estimator of the stack (the combiner).

    meta_estimator : Estimator that will stay on top of the stack.

    **kwargs : Keyword arguments to be passed to `make_stacked_layer`.

    Returns
    -------
    Pipeline

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.ensemble import stack_estimators
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> from sklearn.svm import SVC
    >>> from sklearn.linear_model import LogisticRegression
    >>> eclf = stack_estimators([KNeighborsClassifier(n_neighbors=2), SVC()],
    ...                         [KNeighborsClassifier(n_neighbors=3)],
    ...                         LogisticRegression())
    >>> X = np.array([[1, 3], [.12, 1], [.5, -2], [1, -1], [-2, .1], [7, -84]])
    >>> y = np.array([1, 0, 0, 1, 0, 1])
    >>> eclf.fit(X, y).predict(X)
    array([0, 1, 1, 0, 1, 0])
    """
    if len(estimators) < 1:
        raise ValueError("Invalid number of layers.")

    meta_estimators = []
    for i, stack_row in enumerate(estimators[:-1]):
        meta_estimators.append(("layer%d" % i,
                                make_stack_layer(*stack_row, **kwargs)))

    meta_estimators.append(("combiner", estimators[-1]))

    return Pipeline(meta_estimators)
