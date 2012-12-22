"""
Test the pipeline module.
"""
import numpy as np
from scipy import sparse

from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal

from sklearn.base import BaseEstimator, clone
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.decomposition.pca import PCA, RandomizedPCA
from sklearn.datasets import load_iris
from sklearn.preprocessing import StandardScaler
from sklearn.feature_extraction.text import CountVectorizer


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
    pipe = assert_raises(TypeError, Pipeline, [('svc', IncorrectT)])
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
    scaler = StandardScaler()
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


def test_feature_union():
    # basic sanity check for feature union
    iris = load_iris()
    X = iris.data
    X -= X.mean(axis=0)
    y = iris.target
    pca = RandomizedPCA(n_components=2, random_state=0)
    select = SelectKBest(k=1)
    fs = FeatureUnion([("pca", pca), ("select", select)])
    fs.fit(X, y)
    X_transformed = fs.transform(X)
    assert_equal(X_transformed.shape, (X.shape[0], 3))

    # check if it does the expected thing
    assert_array_almost_equal(X_transformed[:, :-1], pca.fit_transform(X))
    assert_array_equal(X_transformed[:, -1],
                       select.fit_transform(X, y).ravel())

    # test if it also works for sparse input
    # We use a different pca object to control the random_state stream
    fs = FeatureUnion([("pca", pca), ("select", select)])
    X_sp = sparse.csr_matrix(X)
    X_sp_transformed = fs.fit_transform(X_sp, y)
    assert_array_almost_equal(X_transformed, X_sp_transformed.toarray())

    # test setting parameters
    fs.set_params(select__k=2)
    assert_equal(fs.fit_transform(X, y).shape, (X.shape[0], 4))

    # test it works with transformers missing fit_transform
    fs = FeatureUnion([("mock", TransfT()), ("pca", pca), ("select", select)])
    X_transformed = fs.fit_transform(X, y)
    assert_equal(X_transformed.shape, (X.shape[0], 8))


def test_pipeline_transform():
    # Test whether pipeline works with a transformer at the end.
    # Also test pipline.transform and pipeline.inverse_transform
    iris = load_iris()
    X = iris.data
    pca = PCA(n_components=2)
    pipeline = Pipeline([('pca', pca)])

    # test transform and fit_transform:
    X_trans = pipeline.fit(X).transform(X)
    X_trans2 = pipeline.fit_transform(X)
    X_trans3 = pca.fit_transform(X)
    assert_array_almost_equal(X_trans, X_trans2)
    assert_array_almost_equal(X_trans, X_trans3)

    X_back = pipeline.inverse_transform(X_trans)
    X_back2 = pca.inverse_transform(X_trans)
    assert_array_almost_equal(X_back, X_back2)


def test_pipeline_fit_transform():
    # Test whether pipeline works with a transformer missing fit_transform
    iris = load_iris()
    X = iris.data
    y = iris.target
    transft = TransfT()
    pipeline = Pipeline([('mock', transft)])

    # test fit_transform:
    X_trans = pipeline.fit_transform(X, y)
    X_trans2 = transft.fit(X, y).transform(X)
    assert_array_almost_equal(X_trans, X_trans2)


def test_feature_union_weights():
    # test feature union with transformer weights
    iris = load_iris()
    X = iris.data
    y = iris.target
    pca = RandomizedPCA(n_components=2, random_state=0)
    select = SelectKBest(k=1)
    # test using fit followed by transform
    fs = FeatureUnion([("pca", pca), ("select", select)],
                      transformer_weights={"pca": 10})
    fs.fit(X, y)
    X_transformed = fs.transform(X)
    # test using fit_transform
    fs = FeatureUnion([("pca", pca), ("select", select)],
                      transformer_weights={"pca": 10})
    X_fit_transformed = fs.fit_transform(X, y)
    # test it works with transformers missing fit_transform
    fs = FeatureUnion([("mock", TransfT()), ("pca", pca), ("select", select)],
                      transformer_weights={"mock": 10})
    X_fit_transformed_wo_method = fs.fit_transform(X, y)
    # check against expected result

    # We use a different pca object to control the random_state stream
    assert_array_almost_equal(X_transformed[:, :-1], 10 * pca.fit_transform(X))
    assert_array_equal(X_transformed[:, -1],
                       select.fit_transform(X, y).ravel())
    assert_array_almost_equal(X_fit_transformed[:, :-1],
                              10 * pca.fit_transform(X))
    assert_array_equal(X_fit_transformed[:, -1],
                       select.fit_transform(X, y).ravel())
    assert_equal(X_fit_transformed_wo_method.shape, (X.shape[0], 7))


def test_feature_union_feature_names():
    JUNK_FOOD_DOCS = (
        "the pizza pizza beer copyright",
        "the pizza burger beer copyright",
        "the the pizza beer beer copyright",
        "the burger beer beer copyright",
        "the coke burger coke copyright",
        "the coke burger burger",
    )
    word_vect = CountVectorizer(analyzer="word")
    char_vect = CountVectorizer(analyzer="char_wb", ngram_range=(3, 3))
    ft = FeatureUnion([("chars", char_vect), ("words", word_vect)])
    ft.fit(JUNK_FOOD_DOCS)
    feature_names = ft.get_feature_names()
    for feat in feature_names:
        assert_true("chars__" in feat or "words__" in feat)
    assert_equal(len(feature_names), 35)
