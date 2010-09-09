"""
Pipeline: chain transforms and estimators to build a composite estimator.
"""
# Author: Edouard Duchesnay
#         Gael Varoquaux
#         Virgile Fritsch
# Licence: BSD

from .base import BaseEstimator

class Pipeline(BaseEstimator):
    """ Pipeline of transforms with a final estimator 

        Sequentialy apply a list of transforms and a final estimator 
        A transform implements fit & transform methods
        A estimator implements fit & predict methods

        Example
        =======

        >>> from scikits.learn import svm, datasets
        >>> from scikits.learn.datasets import samples_generator
        >>> from scikits.learn.feature_selection import SelectKBest, f_regression
        >>> from scikits.learn.pipeline import Pipeline

        >>> # generate some data to play with
        >>> X, y = samples_generator.test_dataset_classif(k=5)

        >>> # ANOVA SVM-C
        >>> anova_filter = SelectKBest(f_regression, k=5)
        >>> clf = svm.SVC(kernel='linear')

        >>> anova_svm = Pipeline([anova_filter], clf)
        >>> _ = anova_svm.fit(X,y)

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
        self._named_steps = dict(steps)
        names, estimators = zip(*steps)
        self.steps = steps
        assert len(self._named_steps) == len(steps), ("Names provided are "
            "not unique: %s" % names)
        transforms = estimators[:-1]
        estimator = estimators[-1]
        for t in  transforms:
            assert hasattr(t, "fit") and hasattr(t, "transform"), ValueError(
                "All intermediate steps a the chain should be transforms "
                "and implement fit and transform",
                "'%s' (type %s) doesn't)" % (t, type(t))
            )
        assert hasattr(estimator, "fit"), \
            ("Last step of chain should implement fit",
                "'%s' (type %s) doesn't)" % (estimator, type(estimator))
            )

    def _get_params(self, deep=False):
        if not deep:
            return super(Pipeline, self)._get_params(deep=False)
        else:
            out = self._named_steps.copy()
            for name, step in self._named_steps.iteritems():
                for key, value in step._get_params(deep=True).iteritems():
                    out['%s__%s' % (name, key)] = value
        return out
    
    #---------------------------------------------------------------------------
    # Estimator interface
    #---------------------------------------------------------------------------

    def fit(self, X, y=None):
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.fit(Xt, y).transform(Xt)
        self.steps[-1][-1].fit(Xt, y)
        return self

    def predict(self, X):
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict(Xt)

    def score(self, X, y=None):
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].score(Xt, y)

