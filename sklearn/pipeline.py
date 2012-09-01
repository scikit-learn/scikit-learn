"""
The :mod:`sklearn.pipeline` module implements utilites to build a composite
estimator, as a chain of transforms and estimators.
"""
# Author: Edouard Duchesnay
#         Gael Varoquaux
#         Virgile Fritsch
#         Alexandre Gramfort
# Licence: BSD

from .base import BaseEstimator

__all__ = ['Pipeline']


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

    Attributes
    ----------
    `steps` : list of (name, object)
        List of the named object that compose the pipeline, in the \
        order that they are applied on the data.

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

    def __init__(self, steps):
        self.named_steps = dict(steps)
        names, estimators = zip(*steps)
        if len(self.named_steps) != len(steps):
            raise ValueError("Names provided are not unique: %s" % names)

        self.steps = zip(names, estimators)     # shallow copy of steps
        transforms = estimators[:-1]
        estimator = estimators[-1]

        for t in transforms:
            if not (hasattr(t, "fit") or hasattr(t, "fit_transform")) \
              or not hasattr(t, "transform"):
                raise TypeError("All intermediate steps a the chain should "
                        "be transforms and implement fit and transform"
                        "'%s' (type %s) doesn't)" % (t, type(t)))

        if not hasattr(estimator, "fit"):
            raise TypeError("Last step of chain should implement fit "
                "'%s' (type %s) doesn't)" % (estimator, type(estimator)))

    def get_params(self, deep=True):
        if not deep:
            return super(Pipeline, self).get_params(deep=False)
        else:
            out = self.named_steps.copy()
            for name, step in self.named_steps.iteritems():
                for key, value in step.get_params(deep=True).iteritems():
                    out['%s__%s' % (name, key)] = value
            return out

    # Estimator interface

    def _pre_transform(self, X, y=None, **fit_params):
        fit_params_steps = dict((step, {}) for step, _ in self.steps)
        for pname, pval in fit_params.iteritems():
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
        """
        Xt, fit_params = self._pre_transform(X, y, **fit_params)
        self.steps[-1][-1].fit(Xt, y, **fit_params)
        return self

    def fit_transform(self, X, y=None, **fit_params):
        """Fit all the transforms one after the other and transform the
        data, then use fit_transform on transformed data using the final
        estimator. Valid only if the final estimator implements
        fit_transform."""
        Xt, fit_params = self._pre_transform(X, y, **fit_params)
        return self.steps[-1][-1].fit_transform(Xt, y, **fit_params)

    def predict(self, X):
        """Applies transforms to the data, and the predict method of the
        final estimator. Valid only if the final estimator implements
        predict."""
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict(Xt)

    def predict_proba(self, X):
        """Applies transforms to the data, and the predict_proba method of the
        final estimator. Valid only if the final estimator implements
        predict_proba."""
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict_proba(Xt)

    def decision_function(self, X):
        """Applies transforms to the data, and the decision_function method of
        the final estimator. Valid only if the final estimator implements
        decision_function."""
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].decision_function(Xt)

    def predict_log_proba(self, X):
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict_log_proba(Xt)

    def transform(self, X):
        """Applies transforms to the data, and the transform method of the
        final estimator. Valid only if the final estimator implements
        transform."""
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].transform(Xt)

    def inverse_transform(self, X):
        if X.ndim == 1:
            X = X[None, :]
        Xt = X
        for name, step in self.steps[:-1][::-1]:
            Xt = step.inverse_transform(Xt)
        return Xt

    def score(self, X, y=None):
        """Applies transforms to the data, and the score method of the
        final estimator. Valid only if the final estimator implements
        score."""
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].score(Xt, y)

    @property
    def _pairwise(self):
        # check if first estimator expects pairwise input
        return getattr(self.steps[0][1], '_pairwise', False)
