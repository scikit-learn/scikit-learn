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
from .externals.joblib import Parallel, delayed
from .externals import six
from .utils import tosequence
from .utils.metaestimators import if_delegate_has_method
from .externals.six import iteritems

__all__ = ['Pipeline', 'FeatureUnion']


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
    >>> # and a parameter 'C' of the svm
    >>> anova_svm.set_params(anova__k=10, svc__C=.1).fit(X, y)
    ...                                              # doctest: +ELLIPSIS
    Pipeline(steps=[...])
    >>> prediction = anova_svm.predict(X)
    >>> anova_svm.score(X, y)                        # doctest: +ELLIPSIS
    0.77...
    """

    # BaseEstimator interface

    def __init__(self, steps):
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
                raise TypeError("All intermediate steps of the chain should "
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

    @property
    def _final_estimator(self):
        return self.steps[-1][1]

    # Estimator interface

    def _pre_transform(self, X, y=None, **fit_params):
        fit_params_steps = dict((step, {}) for step, _ in self.steps)
        for pname, pval in six.iteritems(fit_params):
            step, param = pname.split('__', 1)
            fit_params_steps[step][param] = pval
        Xt = X
        for name, transform in self.steps[:-1]:
            if hasattr(transform, "fit_transform"):
                Xt = transform.fit_transform(Xt, y, **fit_params_steps[name])
            else:
                Xt = transform.fit(Xt, y, **fit_params_steps[name]) \
                              .transform(Xt)
        return Xt, fit_params_steps[self.steps[-1][0]]

    def fit(self, X, y=None, **fit_params):
        """Fit all the transforms one after the other and transform the
        data, then fit the transformed data using the final estimator.

        Parameters
        ----------
        X : iterable
            Training data. Must fulfill input requirements of first step of the
            pipeline.
        y : iterable, default=None
            Training targets. Must fulfill label requirements for all steps of
            the pipeline.
        """
        Xt, fit_params = self._pre_transform(X, y, **fit_params)
        self.steps[-1][-1].fit(Xt, y, **fit_params)
        return self

    def fit_transform(self, X, y=None, **fit_params):
        """Fit all the transforms one after the other and transform the
        data, then use fit_transform on transformed data using the final
        estimator.

        Parameters
        ----------
        X : iterable
            Training data. Must fulfill input requirements of first step of the
            pipeline.

        y : iterable, default=None
            Training targets. Must fulfill label requirements for all steps of
            the pipeline.
        """
        Xt, fit_params = self._pre_transform(X, y, **fit_params)
        if hasattr(self.steps[-1][-1], 'fit_transform'):
            return self.steps[-1][-1].fit_transform(Xt, y, **fit_params)
        else:
            return self.steps[-1][-1].fit(Xt, y, **fit_params).transform(Xt)

    @if_delegate_has_method(delegate='_final_estimator')
    def predict(self, X):
        """Applies transforms to the data, and the predict method of the
        final estimator. Valid only if the final estimator implements
        predict.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step of
            the pipeline.
        """
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict(Xt)

    @if_delegate_has_method(delegate='_final_estimator')
    def predict_proba(self, X):
        """Applies transforms to the data, and the predict_proba method of the
        final estimator. Valid only if the final estimator implements
        predict_proba.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step of
            the pipeline.
        """
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict_proba(Xt)

    @if_delegate_has_method(delegate='_final_estimator')
    def decision_function(self, X):
        """Applies transforms to the data, and the decision_function method of
        the final estimator. Valid only if the final estimator implements
        decision_function.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step of
            the pipeline.
        """
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].decision_function(Xt)

    @if_delegate_has_method(delegate='_final_estimator')
    def predict_log_proba(self, X):
        """Applies transforms to the data, and the predict_log_proba method of
        the final estimator. Valid only if the final estimator implements
        predict_log_proba.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step of
            the pipeline.
        """
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict_log_proba(Xt)

    @if_delegate_has_method(delegate='_final_estimator')
    def transform(self, X):
        """Applies transforms to the data, and the transform method of the
        final estimator. Valid only if the final estimator implements
        transform.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step of
            the pipeline.
        """
        Xt = X
        for name, transform in self.steps:
            Xt = transform.transform(Xt)
        return Xt

    @if_delegate_has_method(delegate='_final_estimator')
    def inverse_transform(self, X):
        """Applies inverse transform to the data.
        Starts with the last step of the pipeline and applies ``inverse_transform`` in
        inverse order of the pipeline steps.
        Valid only if all steps of the pipeline implement inverse_transform.

        Parameters
        ----------
        X : iterable
            Data to inverse transform. Must fulfill output requirements of the
            last step of the pipeline.
        """
        if X.ndim == 1:
            X = X[None, :]
        Xt = X
        for name, step in self.steps[::-1]:
            Xt = step.inverse_transform(Xt)
        return Xt

    @if_delegate_has_method(delegate='_final_estimator')
    def score(self, X, y=None):
        """Applies transforms to the data, and the score method of the
        final estimator. Valid only if the final estimator implements
        score.

        Parameters
        ----------
        X : iterable
            Data to score. Must fulfill input requirements of first step of the
            pipeline.

        y : iterable, default=None
            Targets used for scoring. Must fulfill label requirements for all steps of
            the pipeline.
        """
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].score(Xt, y)

    @property
    def classes_(self):
        return self.steps[-1][-1].classes_

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
    >>> make_pipeline(StandardScaler(), GaussianNB())    # doctest: +NORMALIZE_WHITESPACE
    Pipeline(steps=[('standardscaler',
                     StandardScaler(copy=True, with_mean=True, with_std=True)),
                    ('gaussiannb', GaussianNB())])

    Returns
    -------
    p : Pipeline
    """
    return Pipeline(_name_estimators(steps))


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
