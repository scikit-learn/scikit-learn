"""
The :mod:`sklearn.pipeline` module implements utilities to build a composite
estimator, as a chain of transforms and estimators.
"""
# Author: Edouard Duchesnay
#         Gael Varoquaux
#         Virgile Fritsch
#         Alexandre Gramfort
# Licence: BSD

from collections import defaultdict, namedtuple

import numpy as np
from scipy import sparse

from .base import BaseEstimator, TransformerMixin
from .externals.joblib import Parallel, delayed, Memory
from .externals import six
from .utils import tosequence
from .externals.six import iteritems

__all__ = ['Pipeline', 'FeatureUnion']


# One round of beers on me if someone finds out why the backslash
# is needed in the Attributes section so as not to upset sphinx.

class Pipeline(BaseEstimator):
    """Pipeline of transforms with a final estimator.

    Sequentially apply a list of transforms and a final estimator.
    Intermediate steps of the pipeline must be 'transforms', that is, they
    must implements fit and transform methods.
    The final estimator needs only implements fit.

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

    memory : string or joblib Memory instance (optional)
        A string value should point to a path to be used for caching.  This may
        add substantial overhead, but will benefit when searching over a
        parameter space in later steps of the pipeline and earlier steps are
        not cheap.

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
    >>> anova_svm = Pipeline([('anova', anova_filter), ('svc', clf)])

    >>> # You can set the parameters using the names issued
    >>> # For instance, fit using a k of 10 in the SelectKBest
    >>> # and a parameter 'C' of the svn
    >>> anova_svm.set_params(anova__k=10, svc__C=.1).fit(X, y)
    ...                                              # doctest: +ELLIPSIS
    Pipeline(steps=[...])

    >>> prediction = anova_svm.predict(X)
    >>> anova_svm.score(X, y)
    0.75
    """

    # BaseEstimator interface

    def __init__(self, steps, memory=None):
        self.memory = memory

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

    def _cached_run(self, method, X, *other_args):
        parameters = self.get_params(True)
        step_params = defaultdict(dict)
        for k, v in list(iteritems(parameters)):
            try:
                step_name, suffix = k.split('__', 1)
            except ValueError:
                # TODO: do something
                continue
            step_params[step_name][suffix] = v
        step_states = [(name, type(estimator), step_params[name])
                       for name, estimator in self.steps]
        if method in ('fit', 'fit_transform'):
            self._fit_args = (X,) + tuple(other_args)
        elif not hasattr(self, '_fit_args'):
            raise 'Pipeline is not fit when calling %r' % method

        if hasattr(self.memory, 'cache'):
            cache = self.memory.cache
        else:
            cache = Memory(self.memory).cache

        @cache
        def _fit(fit_args, step_states):
            """Fit the steps covered by steps_state"""
            X = fit_args[0]
            other_args = fit_args[1:]
            if len(step_states) > 1:
                self.steps[:len(step_states) - 1] = _fit(fit_args,
                                                         step_states[:-1])
                X = _transform(X, other_args, fit_args, step_states[:-1])
            step_name, step_estimator = self.steps[len(step_states) - 1]
            step_estimator.fit(X, *other_args)
            return self.steps[:len(step_states)]

        @cache
        def _transform(X, other_args, fit_args, step_states,
                       method='transform'):
            if len(step_states) > 1:
                X = _transform(X, other_args, fit_args, step_states[:-1])
            _, step_estimator = self.steps[len(step_states) - 1]
            if method == 'transform':
                other_args = ()
            return getattr(step_estimator, method)(X, *other_args)

        if method in ('fit', 'fit_transform'):
            self.steps = _fit(self._fit_args, step_states)
            if method == 'fit':
                return self
            method = 'transform'
        res = _transform(X, other_args, self._fit_args, step_states, method)
        return res

    # Estimator interface

    def fit(self, X, y=None, **fit_params):
        """Fit all the transforms one after the other and transform the
        data, then fit the transformed data using the final estimator.
        """
        return self._cached_run('fit', X, y)

    def fit_transform(self, X, y=None, **fit_params):
        """Fit all the transforms one after the other and transform the
        data, then use fit_transform on transformed data using the final
        estimator."""
        res = self._cached_run('fit_transform', X, y)
        return res

    def predict(self, X):
        """Applies transforms to the data, and the predict method of the
        final estimator. Valid only if the final estimator implements
        predict."""
        return self._cached_run('predict', X)

    def predict_proba(self, X):
        """Applies transforms to the data, and the predict_proba method of the
        final estimator. Valid only if the final estimator implements
        predict_proba."""
        return self._cached_run('predict_proba', X)

    def decision_function(self, X):
        """Applies transforms to the data, and the decision_function method of
        the final estimator. Valid only if the final estimator implements
        decision_function."""
        return self._cached_run('decision_function', X)

    def predict_log_proba(self, X):
        return self._cached_run('predict_log_proba', X)

    def transform(self, X):
        """Applies transforms to the data, and the transform method of the
        final estimator. Valid only if the final estimator implements
        transform."""
        return self._cached_run('transform', X)

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
        return self._cached_run('score', X, y)

    @property
    def _pairwise(self):
        # check if first estimator expects pairwise input
        return getattr(self.steps[0][1], '_pairwise', False)


def _fit_one_transformer(transformer, X, y):
    return transformer.fit(X, y)


def _transform_one(transformer, name, X, transformer_weights):
    if transformer_weights is not None and name in transformer_weights:
        # if we have a weight for this transformer, muliply output
        return transformer.transform(X) * transformer_weights[name]
    return transformer.transform(X)


def _fit_transform_one(transformer, name, X, y, transformer_weights,
                       **fit_params):
    if transformer_weights is not None and name in transformer_weights:
        # if we have a weight for this transformer, muliply output
        if hasattr(transformer, 'fit_transform'):
            X_transformed = transformer.fit_transform(X, y, **fit_params)
            return X_transformed * transformer_weights[name], transformer
        else:
            X_transformed = transformer.fit(X, y, **fit_params).transform(X)
            return X_transformed * transformer_weights[name], transformer
    if hasattr(transformer, 'fit_transform'):
        X_transformed = transformer.fit_transform(X, y, **fit_params)
        return X_transformed, transformer
    else:
        X_transformed = transformer.fit(X, y, **fit_params).transform(X)
        return X_transformed, transformer


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

