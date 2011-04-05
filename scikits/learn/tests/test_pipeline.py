"""
Test the pipeline module.
"""

from nose.tools import assert_raises, assert_equal, assert_false

import numpy as np

from ..base import BaseEstimator, clone
from ..pipeline import Pipeline
from ..svm import SVC
from ..linear_model import LogisticRegression
from ..feature_selection import SelectKBest, f_classif
from ..pca import PCA, RandomizedPCA
from ..datasets import load_iris


class IncorrectT(BaseEstimator):
    """Small class to test parameter dispatching.
    """

    def __init__(self, a=None, b=None):
        self.a = a
        self.b = b


class T(IncorrectT):

    def fit(self, X, y):
        return self


def test_pipeline_init():
    """ Test the various init parameters of the pipeline.
    """
    assert_raises(TypeError, Pipeline)
    # Check that we can't instantiate pipelines with objects without fit
    # method
    pipe = assert_raises(AssertionError, Pipeline,
                        [('svc', IncorrectT)])
    # Smoke test with only an estimator
    clf = T()
    pipe = Pipeline([('svc', clf)])
    assert_equal(pipe._get_params(deep=True),
                 dict(svc__a=None, svc__b=None, svc=clf))

    # Check that params are set
    pipe._set_params(svc__a=0.1)
    assert_equal(clf.a, 0.1)
    # Smoke test the repr:
    repr(pipe)

    # Test with two objects
    clf = SVC()
    filter1 = SelectKBest(f_classif)
    pipe = Pipeline([('anova', filter1), ('svc', clf)])

    # Check that params are set
    pipe._set_params(svc__C=0.1)
    assert_equal(clf.C, 0.1)
    # Smoke test the repr:
    repr(pipe)

    # Check that params are not set when naming them wrong
    assert_raises(AssertionError, pipe._set_params, anova__C=0.1)

    # Test clone
    pipe2 = clone(pipe)
    assert_false(pipe.named_steps['svc'] is pipe2.named_steps['svc'])

    # Check that appart from estimators, the parameters are the same
    params = pipe._get_params()
    params2 = pipe2._get_params()
    # Remove estimators that where copied
    params.pop('svc')
    params.pop('anova')
    params2.pop('svc')
    params2.pop('anova')
    assert_equal(params, params2)


def test_pipeline_methods_anova():
    """ Test the various methods of the pipeline.
    """
    iris = load_iris()
    X = iris.data
    y = iris.target
    # Test with Anova + LogisticRegression
    clf = LogisticRegression()
    filter1 = SelectKBest(f_classif, k=2)
    pipe = Pipeline([('anova', filter1), ('logistic', clf)])
    pipe.fit(X, y)
    pipe.predict(X)
    pipe.predict_proba(X)
    pipe.predict_log_proba(X)
    pipe.score(X, y)


def test_pipeline_methods_pca_svm():
    """Test the various methods of the pipeline."""
    iris = load_iris()
    X = iris.data
    y = iris.target
    # Test with PCA + SVC
    clf = SVC(probability=True)
    pca = PCA(n_components='mle', whiten=True)
    pipe = Pipeline([('pca', pca), ('svc', clf)])
    pipe.fit(X, y)
    pipe.predict(X)
    pipe.predict_proba(X)
    pipe.predict_log_proba(X)
    pipe.score(X, y)


def test_pipeline_methods_randomized_pca_svm():
    """Test the various methods of the pipeline."""
    iris = load_iris()
    X = iris.data
    y = iris.target
    # Test with PCA + SVC
    clf = SVC(probability=True)
    pca = RandomizedPCA(n_components=2, whiten=True)
    pipe = Pipeline([('pca', pca), ('svc', clf)])
    pipe.fit(X, y)
    pipe.predict(X)
    pipe.predict_proba(X)
    pipe.predict_log_proba(X)
    pipe.score(X, y)
