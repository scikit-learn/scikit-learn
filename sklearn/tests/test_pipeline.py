"""
Test the pipeline module.
"""
import numpy as np
from scipy import sparse

from sklearn.externals.six.moves import zip
from sklearn.utils.testing import assert_raises, assert_raises_regex, assert_raise_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_warns_message

from sklearn.base import clone
from sklearn.pipeline import Pipeline, FeatureUnion, make_pipeline, make_union
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.decomposition import PCA, RandomizedPCA, TruncatedSVD
from sklearn.datasets import load_iris
from sklearn.preprocessing import StandardScaler
from sklearn.feature_extraction.text import CountVectorizer


JUNK_FOOD_DOCS = (
    "the pizza pizza beer copyright",
    "the pizza burger beer copyright",
    "the the pizza beer beer copyright",
    "the burger beer beer copyright",
    "the coke burger coke copyright",
    "the coke burger burger",
)


class IncorrectT(object):
    """Small class to test parameter dispatching.
    """

    def __init__(self, a=None, b=None):
        self.a = a
        self.b = b


class T(IncorrectT):

    def fit(self, X, y):
        return self

    def get_params(self, deep=False):
        return {'a': self.a, 'b': self.b}

    def set_params(self, **params):
        self.a = params['a']
        return self


class TransfT(T):

    def transform(self, X, y=None):
        return X

    def inverse_transform(self, X):
        return X


class FitParamT(object):
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
    # Test the various init parameters of the pipeline.
    assert_raises(TypeError, Pipeline)
    # Check that we can't instantiate pipelines with objects without fit
    # method
    pipe = assert_raises(TypeError, Pipeline, [('svc', IncorrectT)])
    # Smoke test with only an estimator
    clf = T()
    pipe = Pipeline([('svc', clf)])
    assert_equal(pipe.get_params(deep=True),
                 dict(svc__a=None, svc__b=None, svc=clf,
                      **pipe.get_params(deep=False)))

    # Check that params are set
    pipe.set_params(svc__a=0.1)
    assert_equal(clf.a, 0.1)
    assert_equal(clf.b, None)
    # Smoke test the repr:
    repr(pipe)

    # Test with two objects
    clf = SVC()
    filter1 = SelectKBest(f_classif)
    pipe = Pipeline([('anova', filter1), ('svc', clf)])

    # Check that we can't use the same stage name twice
    assert_raises(ValueError, Pipeline, [('svc', SVC()), ('svc', SVC())])

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

    # Check that apart from estimators, the parameters are the same
    params = pipe.get_params(deep=True)
    params2 = pipe2.get_params(deep=True)

    for x in pipe.get_params(deep=False):
        params.pop(x)

    for x in pipe2.get_params(deep=False):
        params2.pop(x)

    # Remove estimators that where copied
    params.pop('svc')
    params.pop('anova')
    params2.pop('svc')
    params2.pop('anova')
    assert_equal(params, params2)


def test_pipeline_methods_anova():
    # Test the various methods of the pipeline (anova).
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
    # Test that the pipeline can take fit parameters
    pipe = Pipeline([('transf', TransfT()), ('clf', FitParamT())])
    pipe.fit(X=None, y=None, clf__should_succeed=True)
    # classifier should return True
    assert_true(pipe.predict(None))
    # and transformer params should not be changed
    assert_true(pipe.named_steps['transf'].a is None)
    assert_true(pipe.named_steps['transf'].b is None)


def test_pipeline_raise_set_params_error():
    # Test pipeline raises set params error message for nested models.
    pipe = Pipeline([('cls', LinearRegression())])

    # expected error message
    error_msg = ('Invalid parameter %s for estimator %s. '
                 'Check the list of available parameters '
                 'with `estimator.get_params().keys()`.')

    assert_raise_message(ValueError,
                         error_msg % ('fake', 'Pipeline'),
                         pipe.set_params,
                         fake='nope')

    # nested model check
    assert_raise_message(ValueError,
                         error_msg % ("fake", pipe),
                         pipe.set_params,
                         fake__estimator='nope')


def test_pipeline_methods_pca_svm():
    # Test the various methods of the pipeline (pca + svm).
    iris = load_iris()
    X = iris.data
    y = iris.target
    # Test with PCA + SVC
    clf = SVC(probability=True, random_state=0)
    pca = PCA(n_components='mle', whiten=True)
    pipe = Pipeline([('pca', pca), ('svc', clf)])
    pipe.fit(X, y)
    pipe.predict(X)
    pipe.predict_proba(X)
    pipe.predict_log_proba(X)
    pipe.score(X, y)


def test_pipeline_methods_preprocessing_svm():
    # Test the various methods of the pipeline (preprocessing + svm).
    iris = load_iris()
    X = iris.data
    y = iris.target
    n_samples = X.shape[0]
    n_classes = len(np.unique(y))
    scaler = StandardScaler()
    pca = RandomizedPCA(n_components=2, whiten=True)
    clf = SVC(probability=True, random_state=0, decision_function_shape='ovr')

    for preprocessing in [scaler, pca]:
        pipe = Pipeline([('preprocess', preprocessing), ('svc', clf)])
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


def test_fit_predict_on_pipeline():
    # test that the fit_predict method is implemented on a pipeline
    # test that the fit_predict on pipeline yields same results as applying
    # transform and clustering steps separately
    iris = load_iris()
    scaler = StandardScaler()
    km = KMeans(random_state=0)

    # first compute the transform and clustering step separately
    scaled = scaler.fit_transform(iris.data)
    separate_pred = km.fit_predict(scaled)

    # use a pipeline to do the transform and clustering in one step
    pipe = Pipeline([('scaler', scaler), ('Kmeans', km)])
    pipeline_pred = pipe.fit_predict(iris.data)

    assert_array_almost_equal(pipeline_pred, separate_pred)


def test_fit_predict_on_pipeline_without_fit_predict():
    # tests that a pipeline does not have fit_predict method when final
    # step of pipeline does not have fit_predict defined
    scaler = StandardScaler()
    pca = PCA()
    pipe = Pipeline([('scaler', scaler), ('pca', pca)])
    assert_raises_regex(AttributeError,
                        "'PCA' object has no attribute 'fit_predict'",
                        getattr, pipe, 'fit_predict')


def test_feature_union():
    # basic sanity check for feature union
    iris = load_iris()
    X = iris.data
    X -= X.mean(axis=0)
    y = iris.target
    svd = TruncatedSVD(n_components=2, random_state=0)
    select = SelectKBest(k=1)
    fs = FeatureUnion([("svd", svd), ("select", select)])
    fs.fit(X, y)
    X_transformed = fs.transform(X)
    assert_equal(X_transformed.shape, (X.shape[0], 3))

    # check if it does the expected thing
    assert_array_almost_equal(X_transformed[:, :-1], svd.fit_transform(X))
    assert_array_equal(X_transformed[:, -1],
                       select.fit_transform(X, y).ravel())

    # test if it also works for sparse input
    # We use a different svd object to control the random_state stream
    fs = FeatureUnion([("svd", svd), ("select", select)])
    X_sp = sparse.csr_matrix(X)
    X_sp_transformed = fs.fit_transform(X_sp, y)
    assert_array_almost_equal(X_transformed, X_sp_transformed.toarray())

    # test setting parameters
    fs.set_params(select__k=2)
    assert_equal(fs.fit_transform(X, y).shape, (X.shape[0], 4))

    # test it works with transformers missing fit_transform
    fs = FeatureUnion([("mock", TransfT()), ("svd", svd), ("select", select)])
    X_transformed = fs.fit_transform(X, y)
    assert_equal(X_transformed.shape, (X.shape[0], 8))


def test_make_union():
    pca = PCA()
    mock = TransfT()
    fu = make_union(pca, mock)
    names, transformers = zip(*fu.transformer_list)
    assert_equal(names, ("pca", "transft"))
    assert_equal(transformers, (pca, mock))


def test_pipeline_transform():
    # Test whether pipeline works with a transformer at the end.
    # Also test pipeline.transform and pipeline.inverse_transform
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


def test_make_pipeline():
    t1 = TransfT()
    t2 = TransfT()

    pipe = make_pipeline(t1, t2)
    assert_true(isinstance(pipe, Pipeline))
    assert_equal(pipe.steps[0][0], "transft-1")
    assert_equal(pipe.steps[1][0], "transft-2")

    pipe = make_pipeline(t1, t2, FitParamT())
    assert_true(isinstance(pipe, Pipeline))
    assert_equal(pipe.steps[0][0], "transft-1")
    assert_equal(pipe.steps[1][0], "transft-2")
    assert_equal(pipe.steps[2][0], "fitparamt")


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


def test_feature_union_parallel():
    # test that n_jobs work for FeatureUnion
    X = JUNK_FOOD_DOCS

    fs = FeatureUnion([
        ("words", CountVectorizer(analyzer='word')),
        ("chars", CountVectorizer(analyzer='char')),
    ])

    fs_parallel = FeatureUnion([
        ("words", CountVectorizer(analyzer='word')),
        ("chars", CountVectorizer(analyzer='char')),
    ], n_jobs=2)

    fs_parallel2 = FeatureUnion([
        ("words", CountVectorizer(analyzer='word')),
        ("chars", CountVectorizer(analyzer='char')),
    ], n_jobs=2)

    fs.fit(X)
    X_transformed = fs.transform(X)
    assert_equal(X_transformed.shape[0], len(X))

    fs_parallel.fit(X)
    X_transformed_parallel = fs_parallel.transform(X)
    assert_equal(X_transformed.shape, X_transformed_parallel.shape)
    assert_array_equal(
        X_transformed.toarray(),
        X_transformed_parallel.toarray()
    )

    # fit_transform should behave the same
    X_transformed_parallel2 = fs_parallel2.fit_transform(X)
    assert_array_equal(
        X_transformed.toarray(),
        X_transformed_parallel2.toarray()
    )

    # transformers should stay fit after fit_transform
    X_transformed_parallel2 = fs_parallel2.transform(X)
    assert_array_equal(
        X_transformed.toarray(),
        X_transformed_parallel2.toarray()
    )


def test_feature_union_feature_names():
    word_vect = CountVectorizer(analyzer="word")
    char_vect = CountVectorizer(analyzer="char_wb", ngram_range=(3, 3))
    ft = FeatureUnion([("chars", char_vect), ("words", word_vect)])
    ft.fit(JUNK_FOOD_DOCS)
    feature_names = ft.get_feature_names()
    for feat in feature_names:
        assert_true("chars__" in feat or "words__" in feat)
    assert_equal(len(feature_names), 35)


def test_classes_property():
    iris = load_iris()
    X = iris.data
    y = iris.target

    reg = make_pipeline(SelectKBest(k=1), LinearRegression())
    reg.fit(X, y)
    assert_raises(AttributeError, getattr, reg, "classes_")

    clf = make_pipeline(SelectKBest(k=1), LogisticRegression(random_state=0))
    assert_raises(AttributeError, getattr, clf, "classes_")
    clf.fit(X, y)
    assert_array_equal(clf.classes_, np.unique(y))


def test_X1d_inverse_transform():
    transformer = TransfT()
    pipeline = make_pipeline(transformer)
    X = np.ones(10)
    msg = "1d X will not be reshaped in pipeline.inverse_transform"
    assert_warns_message(FutureWarning, msg, pipeline.inverse_transform, X)
