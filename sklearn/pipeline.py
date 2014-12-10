"""
The :mod:`sklearn.pipeline` module implements utilities to build a composite
estimator, as a chain of transforms and estimators.
"""
# Author: Edouard Duchesnay
#         Gael Varoquaux
#         Virgile Fritsch
#         Alexandre Gramfort
#         Lars Buitinck
# Licence: BSD

from collections import defaultdict

import numpy as np
from scipy import sparse

from .base import BaseEstimator, TransformerMixin
from .externals.joblib import Parallel, delayed, Memory
from .externals import six
from .utils import tosequence
from .externals.six import iteritems

__all__ = ['Pipeline', 'FeatureUnion']

# joblib memory object for caching functions in this module
_memory = Memory('../.Pipeline', verbose=0)


# One round of beers on me if someone finds out why the backslash
# is needed in the Attributes section so as not to upset sphinx.

class Pipeline(BaseEstimator):
    """Pipeline of transforms with a final estimator.

    Sequentially apply a list of transforms and a final estimator.
    Intermediate steps of the pipeline must be 'transforms', that is, they
    must implement fit and transform methods.
    The final estimator only needs to implement fit.

    The purpose of the pipeline is to assemble several steps that can be
    cross-validated together while setting different parameters.
    For this, it enables setting parameters of the various steps using their
    names and the parameter name separated by a '__', as in the example below.

    Parameters
    ----------
    steps: list
        List of (name, transform) tuples (implementing fit/transform) that are
        chained, in the order in which they are chained, with the last object
        an estimator.

    cache: boolean
        Whether results from steps in the pipeline should be cached for future
        computations. This will result in stages in the pipeline NOT BEING
        CALLED if they have previously been passed the same input (even across
        multiple program invokations). This can save significant computation
        by only computing expensive tranforms one time, but will IGNORE changes
        to the workings of any stage. It is therefore recommended only if
        elements of the pipeline are stable code (e.g., classes in sklearn)
        with functionality that will not change during development of your
        program.

    Examples
    --------
    >>> from sklearn import svm
    >>> from sklearn.datasets import samples_generator
    >>> from sklearn.feature_selection import SelectKBest
    >>> from sklearn.feature_selection import f_regression
    >>> from sklearn.pipeline import Pipeline
    >>> # generate some data to play with
    >>> X, y = samples_generator.make_classification(
    ...     n_informative=5, n_redundant=0, random_state=42)
    >>> # ANOVA SVM-C
    >>> anova_filter = SelectKBest(f_regression, k=5)
    >>> clf = svm.SVC(kernel='linear')
    >>> steps = (('anova', anova_filter), ('svc', clf))
    >>> anova_svm = Pipeline(steps)
    >>> # You can set the parameters using the names issued
    >>> # For instance, fit using a k of 10 in the SelectKBest
    >>> # and a parameter 'C' of the svm
    >>> anova_svm.set_params(anova__k=10, svc__C=.1).fit(X, y)
    ...                        # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Pipeline(cache=False, steps=[...])
    >>> prediction = anova_svm.predict(X)
    >>> anova_svm.score(X, y)                        # doctest: +ELLIPSIS
    0.77...
    >>> # we may want these results again in the future--let's
    >>> # persistently cache them
    >>> caching_anova_svm = Pipeline(steps, cache=False)
    >>> caching_anova_svm.set_params(anova__k=10, svc__C=.1).fit(X, y)
    ...                        # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Pipeline(cache=False, steps=[...])
    >>> prediction2 = caching_anova_svm.predict(X)
    >>> all(prediction == prediction2)
    True
    >>> # now suppose we need to compute this somewhere else
    >>> # in our code, or are running this program a separate time
    >>> prediction3 = caching_anova_svm.predict(X)  # identical, but faster
    >>> all(prediction == prediction3)
    True
    >>> # the cached results will be retrieved even if we construct a
    >>> # new object, as long as the steps, parameters, and arguments
    >>> # are the same
    >>> caching_anova_svm2 = Pipeline(steps, cache=False)
    >>> caching_anova_svm2.set_params(anova__k=10, svc__C=.1).fit(X, y)
    ...                        # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Pipeline(cache=False, steps=[...])
    >>> prediction4 = caching_anova_svm2.predict(X)
    >>> all(prediction == prediction4)
    True
    """

    # BaseEstimator interface

    def __init__(self, steps, cache=False):
        self.cache = cache

        self.named_steps = dict(steps)
        names, estimators = zip(*steps)
        if len(self.named_steps) != len(steps):
            raise ValueError("Names provided are not unique: %s" % (names,))

        # shallow copy of steps
        self.steps = tosequence(zip(names, estimators))
        transforms = estimators[:-1]
        estimator = estimators[-1]

        for t in transforms:
            if (not (hasattr(t, "fit") or hasattr(t, "fit_transform")) or not
                    hasattr(t, "transform")):
                raise TypeError("All intermediate steps a the chain should "
                                "be transforms and implement fit and transform"
                                " '%s' (type %s) doesn't)" % (t, type(t)))

        if not hasattr(estimator, "fit"):
            raise TypeError("Last step of chain should implement fit "
                            "'%s' (type %s) doesn't)"
                            % (estimator, type(estimator)))

    def get_params(self, deep=True):
        if not deep:
            return super(Pipeline, self).get_params(deep=False)
        else:
            out = self.named_steps.copy()
            for name, step in six.iteritems(self.named_steps):
                for key, value in six.iteritems(step.get_params(deep=True)):
                    out['%s__%s' % (name, key)] = value
            return out

    # Estimator interface

    def _pre_transform(self, X, y=None, **fit_params):
        fit_params_steps = dict((step, {}) for step, _ in self.steps)
        for pname, pval in six.iteritems(fit_params):
            step, param = pname.split('__', 1)
            fit_params_steps[step][param] = pval

        Xt = X
        for i, step in enumerate(self.steps[:-1]):
            name, transform = step
            fit_params = fit_params_steps[name]
            Xt, fitted = self._fit_transform(transform, Xt, y, **fit_params)

            # store fitted estimator
            self.steps[i] = (name, fitted)
            self.named_steps[name] = fitted

        return Xt, fit_params_steps[self.steps[-1][0]]

    def fit(self, X, y=None, **fit_params):
        """Fit all the transforms one after the other and transform the
        data, then fit the transformed data using the final estimator.
        """
        Xt, fit_params = self._pre_transform(X, y, **fit_params)
        lastEstimator = self._get_last_estimator()
        estimator = self._fit(lastEstimator, Xt, y, **fit_params)
        self._set_last_estimator(estimator)
        return self

    def fit_transform(self, X, y=None, **fit_params):
        """Fit all the transforms one after the other and transform the
        data, then use fit_transform on transformed data using the final
        estimator."""
        Xt, fit_params = self._pre_transform(X, y, **fit_params)
        lastEstimator = self._get_last_estimator()
        Xt, estimator = self._fit_transform(lastEstimator, Xt, y, **fit_params)
        self._set_last_estimator(estimator)
        return Xt

    def predict(self, X):
        """Applies transforms to the data, and the predict method of the
        final estimator. Valid only if the final estimator implements
        predict."""
        Xt = self._transform_except_last_estimator(X)
        return self._get_last_estimator().predict(Xt)

    def predict_proba(self, X):
        """Applies transforms to the data, and the predict_proba method of the
        final estimator. Valid only if the final estimator implements
        predict_proba."""
        Xt = self._transform_except_last_estimator(X)
        return self._get_last_estimator().predict_proba(Xt)

    def decision_function(self, X):
        """Applies transforms to the data, and the decision_function method of
        the final estimator. Valid only if the final estimator implements
        decision_function."""
        Xt = self._transform_except_last_estimator(X)
        return self._get_last_estimator().decision_function(Xt)

    def predict_log_proba(self, X):
        Xt = self._transform_except_last_estimator(X)
        return self._get_last_estimator().predict_log_proba(Xt)

    def transform(self, X):
        """Applies transforms to the data, and the transform method of the
        final estimator. Valid only if the final estimator implements
        transform."""
        return self._apply_transform_steps(X, self.steps)

    def inverse_transform(self, X):
        if X.ndim == 1:
            X = X[None, :]
        Xt = X
        for name, step in self.steps[::-1]:
            Xt = step.inverse_transform(Xt)
        return Xt

    def score(self, X, y=None):
        """Applies transforms to the data, and the score method of the
        final estimator. Valid only if the final estimator implements
        score."""
        Xt = self._transform_except_last_estimator(X)
        return self._get_last_estimator().score(Xt, y)

    def _apply_transform_steps(self, X, steps):
        if not steps:
            return X
        if len(steps) == 1:
            # ensure that we have an iterable of transforms
            _, transform = steps[0]
            transforms = [transform]
        else:
            _, transforms = zip(*steps)
        if self.cache:
            return _cached_apply_transforms(X, transforms)
        return _apply_transforms(X, transforms)

    def _transform_except_last_estimator(self, X):
        return self._apply_transform_steps(X, self.steps[:-1])

    def _fit(self, estimator, X, y, **fit_params):
        if self.cache:
            return _cached_fit(estimator, X, y, **fit_params)
        return _fit(estimator, X, y, **fit_params)

    def _fit_transform(self, transform, X, y, **fit_params):
        if self.cache:
            return _cached_fit_transform(transform, X, y, **fit_params)
        return _fit_transform(transform, X, y, **fit_params)

    def _get_last_estimator(self):
        # each step is a (name, estimator) pair; we want the
        # last element of the last step
        return self.steps[-1][-1]

    def _set_last_estimator(self, estimator):
        name, _ = self.steps[-1]
        self.steps[-1] = (name, estimator)
        self.named_steps[name] = estimator

    @property
    def classes_(self):
        return self._get_last_estimator().classes_

    @property
    def _pairwise(self):
        # check if first estimator expects pairwise input
        return getattr(self.steps[0][1], '_pairwise', False)


def _name_estimators(estimators):
    """Generate names for estimators."""

    names = [type(estimator).__name__.lower() for estimator in estimators]
    namecount = defaultdict(int)
    for est, name in zip(estimators, names):
        namecount[name] += 1

    for k, v in list(six.iteritems(namecount)):
        if v == 1:
            del namecount[k]

    for i in reversed(range(len(estimators))):
        name = names[i]
        if name in namecount:
            names[i] += "-%d" % namecount[name]
            namecount[name] -= 1

    return list(zip(names, estimators))


def make_pipeline(*steps):
    """Construct a Pipeline from the given estimators.

    This is a shorthand for the Pipeline constructor; it does not require, and
    does not permit, naming the estimators. Instead, they will be given names
    automatically based on their types.

    Examples
    --------
    >>> from sklearn.naive_bayes import GaussianNB
    >>> from sklearn.preprocessing import StandardScaler
    >>> make_pipeline(StandardScaler(), GaussianNB())
    ...                     # doctest: +NORMALIZE_WHITESPACE
    Pipeline(cache=False,
         steps=[('standardscaler', StandardScaler(copy=True,
            with_mean=True, with_std=True)),
            ('gaussiannb', GaussianNB())])

    Returns
    -------
    p : Pipeline
    """
    return Pipeline(_name_estimators(steps))


def _fit(estimator, X, y, **fit_params):
    estimator.fit(X, y, **fit_params)
    return estimator


@_memory.cache
def _cached_fit(estimator, X, y, **fit_params):
    estimator.fit(X, y, **fit_params)
    return estimator


def _fit_transform(transformer, X, y, **fit_params):
    if hasattr(transformer, 'fit_transform'):
        X_transformed = transformer.fit_transform(X, y, **fit_params)
        return X_transformed, transformer
    else:
        X_transformed = transformer.fit(X, y, **fit_params).transform(X)
        return X_transformed, transformer


@_memory.cache
def _cached_fit_transform(transformer, X, y, **fit_params):
    if hasattr(transformer, 'fit_transform'):
        X_transformed = transformer.fit_transform(X, y, **fit_params)
        return X_transformed, transformer
    else:
        X_transformed = transformer.fit(X, y, **fit_params).transform(X)
        return X_transformed, transformer


def _apply_transforms(X, transforms):
    if not transforms:
        return X
    Xt = X
    for t in transforms:
        Xt = t.transform(Xt)
    return Xt


@_memory.cache
def _cached_apply_transforms(X, transforms):
    # if transforms is None or empty, just return X
    if not transforms:
        return X
    # otherwise, recurse; implementing this recursively
    # allows us to load cached results from as many
    # previous stages as possible
    Xt = _cached_apply_transforms(X, transforms[:-1])
    return transforms[-1].transform(Xt)


def _fit_transform_one(transformer, name, X, y, transformer_weights,
                       **fit_params):
    Xt, fitted_transformer = _fit_transform(transformer, X, y, **fit_params)
    if transformer_weights is not None and name in transformer_weights:
        return Xt * transformer_weights[name], fitted_transformer
    return Xt, fitted_transformer


def _fit_one_transformer(transformer, X, y):
    return transformer.fit(X, y)


def _transform_one(transformer, name, X, transformer_weights):
    if transformer_weights is not None and name in transformer_weights:
        # if we have a weight for this transformer, muliply output
        return transformer.transform(X) * transformer_weights[name]
    return transformer.transform(X)


class FeatureUnion(BaseEstimator, TransformerMixin):
    """Concatenates results of multiple transformer objects.

    This estimator applies a list of transformer objects in parallel to the
    input data, then concatenates the results. This is useful to combine
    several feature extraction mechanisms into a single transformer.

    Parameters
    ----------
    transformer_list: list of (string, transformer) tuples
        List of transformer objects to be applied to the data. The first
        half of each tuple is the name of the transformer.

    n_jobs: int, optional
        Number of jobs to run in parallel (default 1).

    transformer_weights: dict, optional
        Multiplicative weights for features per transformer.
        Keys are transformer names, values the weights.

    """
    def __init__(self, transformer_list, n_jobs=1, transformer_weights=None):
        self.transformer_list = transformer_list
        self.n_jobs = n_jobs
        self.transformer_weights = transformer_weights

    def get_feature_names(self):
        """Get feature names from all transformers.

        Returns
        -------
        feature_names : list of strings
            Names of the features produced by transform.
        """
        feature_names = []
        for name, trans in self.transformer_list:
            if not hasattr(trans, 'get_feature_names'):
                raise AttributeError("Transformer %s does not provide"
                                     " get_feature_names." % str(name))
            feature_names.extend([name + "__" + f for f in
                                  trans.get_feature_names()])
        return feature_names

    def fit(self, X, y=None):
        """Fit all transformers using X.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            Input data, used to fit transformers.
        """
        transformers = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_one_transformer)(trans, X, y)
            for name, trans in self.transformer_list)
        self._update_transformer_list(transformers)
        return self

    def fit_transform(self, X, y=None, **fit_params):
        """Fit all transformers using X, transform the data and concatenate
        results.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            Input data to be transformed.

        Returns
        -------
        X_t : array-like or sparse matrix, shape (n_samples, sum_n_components)
            hstack of results of transformers. sum_n_components is the
            sum of n_components (output dimension) over transformers.
        """
        result = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_transform_one)(trans, name, X, y,
                                        self.transformer_weights, **fit_params)
            for name, trans in self.transformer_list)

        Xs, transformers = zip(*result)
        self._update_transformer_list(transformers)
        if any(sparse.issparse(f) for f in Xs):
            Xs = sparse.hstack(Xs).tocsr()
        else:
            Xs = np.hstack(Xs)
        return Xs

    def transform(self, X):
        """Transform X separately by each transformer, concatenate results.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            Input data to be transformed.

        Returns
        -------
        X_t : array-like or sparse matrix, shape (n_samples, sum_n_components)
            hstack of results of transformers. sum_n_components is the
            sum of n_components (output dimension) over transformers.
        """
        Xs = Parallel(n_jobs=self.n_jobs)(
            delayed(_transform_one)(trans, name, X, self.transformer_weights)
            for name, trans in self.transformer_list)
        if any(sparse.issparse(f) for f in Xs):
            Xs = sparse.hstack(Xs).tocsr()
        else:
            Xs = np.hstack(Xs)
        return Xs

    def get_params(self, deep=True):
        if not deep:
            return super(FeatureUnion, self).get_params(deep=False)
        else:
            out = dict(self.transformer_list)
            for name, trans in self.transformer_list:
                for key, value in iteritems(trans.get_params(deep=True)):
                    out['%s__%s' % (name, key)] = value
            return out

    def _update_transformer_list(self, transformers):
        self.transformer_list[:] = [
            (name, new)
            for ((name, old), new) in zip(self.transformer_list, transformers)
        ]


# XXX it would be nice to have a keyword-only n_jobs argument to this function,
# but that's not allowed in Python 2.x.
def make_union(*transformers):
    """Construct a FeatureUnion from the given transformers.

    This is a shorthand for the FeatureUnion constructor; it does not require,
    and does not permit, naming the transformers. Instead, they will be given
    names automatically based on their types. It also does not allow weighting.

    Examples
    --------
    >>> from sklearn.decomposition import PCA, TruncatedSVD
    >>> make_union(PCA(), TruncatedSVD())    # doctest: +NORMALIZE_WHITESPACE
    FeatureUnion(n_jobs=1,
                 transformer_list=[('pca', PCA(copy=True, n_components=None,
                                               whiten=False)),
                                   ('truncatedsvd',
                                    TruncatedSVD(algorithm='randomized',
                                                 n_components=2, n_iter=5,
                                                 random_state=None, tol=0.0))],
                 transformer_weights=None)

    Returns
    -------
    f : FeatureUnion
    """
    return FeatureUnion(_name_estimators(transformers))
