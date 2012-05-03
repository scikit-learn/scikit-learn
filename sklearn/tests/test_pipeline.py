"""
Test the pipeline module.
"""
import numpy as np

from nose.tools import assert_raises, assert_equal, assert_false, assert_true

from sklearn.base import BaseEstimator, clone
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.decomposition.pca import PCA, RandomizedPCA
from sklearn.datasets import load_iris
from sklearn.preprocessing import Scaler


class IncorrectT(BaseEstimator):
    """Small class to test parameter dispatching.
    """

    def __init__(self, a=None, b=None):
        self.a = a
        self.b = b


class T(IncorrectT):

    def fit(self, X, y):
        return self


class TransfT(T):

    def transform(self, X, y=None):
        return X


class FitParamT(BaseEstimator):
    """Mock classifier
    """

    def __init__(self):
        self.successful = False
        pass

    def fit(self, X, y, should_succeed=False):
        self.successful = should_succeed

    def predict(self, X):
        return self.successful


def test_pipeline_init():
    """ Test the various init parameters of the pipeline.
    """
    assert_raises(TypeError, Pipeline)
    # Check that we can't instantiate pipelines with objects without fit
    # method
    pipe = assert_raises(TypeError, Pipeline,
                        [('svc', IncorrectT)])
    # Smoke test with only an estimator
    clf = T()
    pipe = Pipeline([('svc', clf)])
    assert_equal(pipe.get_params(deep=True),
                 dict(svc__a=None, svc__b=None, svc=clf))

    # Check that params are set
    pipe.set_params(svc__a=0.1)
    assert_equal(clf.a, 0.1)
    # Smoke test the repr:
    repr(pipe)

    # Test with two objects
    clf = SVC()
    filter1 = SelectKBest(f_classif)
    pipe = Pipeline([('anova', filter1), ('svc', clf)])

    # Check that params are set
    pipe.set_params(svc__C=0.1)
    assert_equal(clf.C, 0.1)
    # Smoke test the repr:
    repr(pipe)

    # Check that params are not set when naming them wrong
    assert_raises(ValueError, pipe.set_params, anova__C=0.1)

    # Test clone
    pipe2 = clone(pipe)
    assert_false(pipe.named_steps['svc'] is pipe2.named_steps['svc'])

    # Check that appart from estimators, the parameters are the same
    params = pipe.get_params()
    params2 = pipe2.get_params()
    # Remove estimators that where copied
    params.pop('svc')
    params.pop('anova')
    params2.pop('svc')
    params2.pop('anova')
    assert_equal(params, params2)


def test_pipeline_methods_anova():
    """ Test the various methods of the pipeline (anova).
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


def test_pipeline_fit_params():
    """Test that the pipeline can take fit parameters
    """
    pipe = Pipeline([('transf', TransfT()), ('clf', FitParamT())])
    pipe.fit(X=None, y=None, clf__should_succeed=True)
    # classifier should return True
    assert_true(pipe.predict(None))
    # and transformer params should not be changed
    assert_true(pipe.named_steps['transf'].a is None)
    assert_true(pipe.named_steps['transf'].b is None)


def test_pipeline_methods_pca_svm():
    """Test the various methods of the pipeline (pca + svm)."""
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


def test_pipeline_methods_preprocessing_svm():
    """Test the various methods of the pipeline (preprocessing + svm)."""
    iris = load_iris()
    X = iris.data
    y = iris.target
    n_samples = X.shape[0]
    n_classes = len(np.unique(y))
    scaler = Scaler()
    pca = RandomizedPCA(n_components=2, whiten=True)
    clf = SVC(probability=True)

    for preprocessing in [scaler, pca]:
        pipe = Pipeline([('scaler', scaler), ('svc', clf)])
        pipe.fit(X, y)

        # check shapes of various prediction functions
        predict = pipe.predict(X)
        assert_equal(predict.shape, (n_samples,))

        proba = pipe.predict_proba(X)
        assert_equal(proba.shape, (n_samples, n_classes))

        log_proba = pipe.predict_log_proba(X)
        assert_equal(log_proba.shape, (n_samples, n_classes))

        decision_function = pipe.decision_function(X)
        assert_equal(decision_function.shape, (n_samples, n_classes))

        pipe.score(X, y)
