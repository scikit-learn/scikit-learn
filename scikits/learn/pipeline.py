"""
Pipeline: chain transforms and estimators to build a composite estimator.
"""
# Author: Edouard Duchesnay
#         Gael Varoquaux
#         Virgile Fritsch
#         Alexandre Gramfort
# Licence: BSD

from .base import BaseEstimator


class Pipeline(BaseEstimator):
    """ Pipeline of transforms with a final estimator

        Sequentialy apply a list of transforms and a final estimator
        Intermediate steps of the pipeline must be 'transforms', that
        is that they must implements fit & transform methods
        The final estimator need only implements fit.

        The purpose of the pipeline is to assemble several steps that can
        be cross-validated together while setting different parameters.
        For this, it enables to setting parameters of the various steps
        using their names and the parameter name separated by a '__',
        as in the example below.

        Attributes
        ===========
        steps: list of (names, object)
            List of the named object that compose the pipeline, in the
            order that they are applied on the data.

        Methods
        =======
        fit:
            Fit all the transforms one after the other and transform the
            data, then fit the transformed data using the final estimator
        fit_transform:
            Fit all the transforms one after the other and transform the
            data, then use fit_transform on transformed data using the final
            estimator. Valid only if the final estimator implements
            fit_transform.
        predict:
            Applies transforms to the data, and the predict method of the
            final estimator. Valid only if the final estimator implements
            predict.
        transform:
            Applies transforms to the data, and the transform method of the
            final estimator. Valid only if the final estimator implements
            transform.
        score:
            Applies transforms to the data, and the score method of the
            final estimator. Valid only if the final estimator implements
            score.


        Example
        =======

        >>> from scikits.learn import svm
        >>> from scikits.learn.datasets import samples_generator
        >>> from scikits.learn.feature_selection import SelectKBest, f_regression
        >>> from scikits.learn.pipeline import Pipeline

        >>> # generate some data to play with
        >>> X, y = samples_generator.test_dataset_classif(k=5)

        >>> # ANOVA SVM-C
        >>> anova_filter = SelectKBest(f_regression, k=5)
        >>> clf = svm.SVC(kernel='linear')
        >>> anova_svm = Pipeline([('anova', anova_filter), ('svc', clf)])

        >>> # You can set the parameters using the names issued
        >>> # For instance, fit using a k of 10 in the SelectKBest
        >>> # and a parameter 'C' of the svn
        >>> anova_svm.fit(X, y, anova__k=10, svc__C=.1) #doctest: +ELLIPSIS
        Pipeline(steps=[('anova', SelectKBest(k=10, score_func=<function f_regression at ...>)), ('svc', SVC(kernel='linear', C=0.1, probability=False, degree=3, coef0=0.0, tol=0.001,
          cache_size=100.0, shrinking=True, gamma=0.0))])

        >>> prediction = anova_svm.predict(X)
        >>> score = anova_svm.score(X)
    """

    #---------------------------------------------------------------------------
    # BaseEstimator interface
    #---------------------------------------------------------------------------

    def __init__(self, steps):
        """
        Parameters
        ==========
        steps: list
            List of (name, transform) object (implementing
            fit/transform) that are chained, in the order in which
            they are chained, with the last object an estimator.
        """
        self.named_steps = dict(steps)
        names, estimators = zip(*steps)
        self.steps = steps
        assert len(self.named_steps) == len(steps), ("Names provided are "
            "not unique: %s" % names)
        transforms = estimators[:-1]
        estimator = estimators[-1]
        for t in  transforms:
            assert (hasattr(t, "fit") or hasattr(t, "fit_transform")) and \
                    hasattr(t, "transform"), ValueError(
                "All intermediate steps a the chain should be transforms "
                "and implement fit and transform",
                "'%s' (type %s) doesn't)" % (t, type(t))
            )
        assert hasattr(estimator, "fit"), \
            ("Last step of chain should implement fit",
                "'%s' (type %s) doesn't)" % (estimator, type(estimator))
            )

    def _get_params(self, deep=True):
        if not deep:
            return super(Pipeline, self)._get_params(deep=False)
        else:
            out = self.named_steps.copy()
            for name, step in self.named_steps.iteritems():
                for key, value in step._get_params(deep=True).iteritems():
                    out['%s__%s' % (name, key)] = value
            return out

    #---------------------------------------------------------------------------
    # Estimator interface
    #---------------------------------------------------------------------------


    def _pre_transform(self, X, y=None, **params):
        self._set_params(**params)
        Xt = X
        for name, transform in self.steps[:-1]:
            if hasattr(transform, "fit_transform"):
                Xt = transform.fit_transform(Xt, y)
            else:
                Xt = transform.fit(Xt, y).transform(Xt)
        return Xt

    def fit(self, X, y=None, **params):
        Xt = self._pre_transform(X, y, **params)
        self.steps[-1][-1].fit(Xt, y)
        return self

    def fit_transform(self, X, y=None, **params):
        Xt = self._pre_transform(X, y, **params)
        return self.steps[-1][-1].fit_transform(Xt, y)

    def predict(self, X):
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict(Xt)

    def predict_proba(self, X):
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict_proba(Xt)

    def predict_log_proba(self, X):
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict_log_proba(Xt)

    def transform(self, X):
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
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].score(Xt, y)
