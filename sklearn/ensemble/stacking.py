"""Stacking API"""

# Author: Caio Oliveira <caioaao@gmail.com>
# License: BSD 3 clause

from ..base import BaseEstimator, TransformerMixin
from ..model_selection import cross_val_predict
from ..pipeline import make_union, make_pipeline


class BlendedClassifierTransformer(BaseEstimator, TransformerMixin):
    """ Transformer to turn estimators into blended estimators

    This is used for stacking models. Blending a classifier prevents data leaks
    between stacks. Blending will happen only when calling `fit_transform`, as
    it's the only stage where this makes sense.

    This is probably useless by itself. See `make_stack_layer` for more info on
    how to use this API for stacking.


    Parameters
    ----------
    clf: the estimator to be blended.

    cv: cv to be used. Will be passed to `cross_val_predict`

    method: string, optional, default: 'auto'
        Invokes the passed method name of the passed estimator. If the method is
        `auto`, will try to invoke `predict_proba` or `predict` in that order.
    """
    def __init__(self, clf, cv=3, method='auto'):
        self.clf = clf
        self.cv = cv
        self.method = method

    def fit(self, *args, **kwargs):
        self.clf = self.clf.fit(*args, **kwargs)
        return self

    def _method_name(self):
        if self.method == 'auto':
            if getattr(self.clf, 'predict_proba', None):
                method = 'predict_proba'
            else:
                method = 'predict'
        else:
            method = self.method

        return method

    def transform(self, *args, **kwargs):
        return getattr(self.clf, self._method_name())(*args, **kwargs)

    def fit_transform(self, X, y):
        preds = cross_val_predict(self.clf, X, y, cv=self.cv,
                                  method=self._method_name())

        self.clf.fit(X, y)

        return preds


def make_stack_layer(*clfs, **kwargs):
    """ Construct a single layer for a stacked model.

    This is a wrapper around pipeline API to make creating a stacked classifier
    more convenient.

    Parameters
    ----------
    *clfs: list of base estimators.

    **kwargs: Keyword arguments to be passed to `make_union`.

    Returns
    -------
    f : FeatureUnion with every base estimator wrapped in a
         `BlendedClassifierTransformer`.

    Examples
    --------
    >>> from sklearn.ensemble import make_stack_layer
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> from sklearn.svm import SVC
    >>> make_stack_layer([KNeighborsClassifier(), SVC()])  # doctest:  +NORMALIZE_WHITESPACE
    FeatureUnion(n_jobs=1,
           transformer_list=[('blendedclassifiertransformer', BlendedClassifierTransformer(clf=[KNeighborsClassifier(algorithm='auto', leaf_size=30, metric='minkowski',
               metric_params=None, n_jobs=1, n_neighbors=5, p=2,
               weights='uniform'), SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
      decision_function_shape='ovr', degree=3, gamma='auto', kernel='rbf',
      max_iter=-1, probability=False, random_state=None, shrinking=True,
      tol=0.001, verbose=False)],
                   cv=3))],
           transformer_weights=None)
    """
    return make_union(*[BlendedClassifierTransformer(clf) for clf in clfs],
                      **kwargs)


def make_stacked_classifier(clf_matrix, meta_estimator, **kwargs):
    """ Construct a stacked classifier

    This wraps `pipeline` API to make creating a stacked classifier more
    convenient.

    Parameters
    ----------
    clf_matrix: 2D matrix with base classifiers. Each row will be turned into a
        layer in the stack.

    meta_estimator: Estimator that will stay on top of the stack.

    **kwargs: Keyword arguments to be passed to `make_stacked_layer.

    Returns
    -------
    p: Pipeline

    Examples
    --------
    >>> from sklearn.ensemble import make_stacked_classifier
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> from sklearn.svm import SVC
    >>> from sklearn.linear_model import LogisticRegression
    >>> clf = make_stacked_classifier([[KNeighborsClassifier(), SVC()], [KNeighborsClassifier(), SVC()]], LogisticRegression())
    >>> print map(lambda x: x[0], clf.steps)
    ['featureunion-1', 'featureunion-2', 'logisticregression']
    """
    clfs = [make_stack_layer(*row, **kwargs)
            for row in clf_matrix]
    clfs.append(meta_estimator)
    return make_pipeline(*clfs)
