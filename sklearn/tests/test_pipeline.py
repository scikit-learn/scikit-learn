"""
Test the pipeline module.
"""
from distutils.version import LooseVersion
from tempfile import mkdtemp
import shutil
import time

import pytest
import numpy as np
from scipy import sparse

from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_dict_equal
from sklearn.utils.testing import assert_no_warnings

from sklearn.base import clone, BaseEstimator
from sklearn.pipeline import Pipeline, FeatureUnion, make_pipeline, make_union
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression, Lasso
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.dummy import DummyRegressor
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.datasets import load_iris
from sklearn.preprocessing import StandardScaler
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.utils._joblib import Memory
from sklearn.utils._joblib import __version__ as joblib_version


JUNK_FOOD_DOCS = (
    "the pizza pizza beer copyright",
    "the pizza burger beer copyright",
    "the the pizza beer beer copyright",
    "the burger beer beer copyright",
    "the coke burger coke copyright",
    "the coke burger burger",
)


class NoFit(object):
    """Small class to test parameter dispatching.
    """

    def __init__(self, a=None, b=None):
        self.a = a
        self.b = b


class NoTrans(NoFit):

    def fit(self, X, y):
        return self

    def get_params(self, deep=False):
        return {'a': self.a, 'b': self.b}

    def set_params(self, **params):
        self.a = params['a']
        return self


class NoInvTransf(NoTrans):
    def transform(self, X):
        return X


class Transf(NoInvTransf):
    def transform(self, X):
        return X

    def inverse_transform(self, X):
        return X


class TransfFitParams(Transf):

    def fit(self, X, y, **fit_params):
        self.fit_params = fit_params
        return self


class Mult(BaseEstimator):
    def __init__(self, mult=1):
        self.mult = mult

    def fit(self, X, y):
        return self

    def transform(self, X):
        return np.asarray(X) * self.mult

    def inverse_transform(self, X):
        return np.asarray(X) / self.mult

    def predict(self, X):
        return (np.asarray(X) * self.mult).sum(axis=1)

    predict_proba = predict_log_proba = decision_function = predict

    def score(self, X, y=None):
        return np.sum(X)


class FitParamT(BaseEstimator):
    """Mock classifier
    """

    def __init__(self):
        self.successful = False

    def fit(self, X, y, should_succeed=False):
        self.successful = should_succeed

    def predict(self, X):
        return self.successful

    def fit_predict(self, X, y, should_succeed=False):
        self.fit(X, y, should_succeed=should_succeed)
        return self.predict(X)

    def score(self, X, y=None, sample_weight=None):
        if sample_weight is not None:
            X = X * sample_weight
        return np.sum(X)


class DummyTransf(Transf):
    """Transformer which store the column means"""

    def fit(self, X, y):
        self.means_ = np.mean(X, axis=0)
        # store timestamp to figure out whether the result of 'fit' has been
        # cached or not
        self.timestamp_ = time.time()
        return self


class DummyEstimatorParams(BaseEstimator):
    """Mock classifier that takes params on predict"""

    def fit(self, X, y):
        return self

    def predict(self, X, got_attribute=False):
        self.got_attribute = got_attribute
        return self


def test_pipeline_init():
    # Test the various init parameters of the pipeline.
    assert_raises(TypeError, Pipeline)
    # Check that we can't instantiate pipelines with objects without fit
    # method
    assert_raises_regex(TypeError,
                        'Last step of Pipeline should implement fit '
                        'or be the string \'passthrough\''
                        '.*NoFit.*',
                        Pipeline, [('clf', NoFit())])
    # Smoke test with only an estimator
    clf = NoTrans()
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

    # Check that we can't instantiate with non-transformers on the way
    # Note that NoTrans implements fit, but not transform
    assert_raises_regex(TypeError,
                        'All intermediate steps should be transformers'
                        '.*\\bNoTrans\\b.*',
                        Pipeline, [('t', NoTrans()), ('svc', clf)])

    # Check that params are set
    pipe.set_params(svc__C=0.1)
    assert_equal(clf.C, 0.1)
    # Smoke test the repr:
    repr(pipe)

    # Check that params are not set when naming them wrong
    assert_raises(ValueError, pipe.set_params, anova__C=0.1)

    # Test clone
    pipe2 = assert_no_warnings(clone, pipe)
    assert not pipe.named_steps['svc'] is pipe2.named_steps['svc']

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


def test_pipeline_init_tuple():
    # Pipeline accepts steps as tuple
    X = np.array([[1, 2]])
    pipe = Pipeline((('transf', Transf()), ('clf', FitParamT())))
    pipe.fit(X, y=None)
    pipe.score(X)

    pipe.set_params(transf='passthrough')
    pipe.fit(X, y=None)
    pipe.score(X)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
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
    pipe = Pipeline([('transf', Transf()), ('clf', FitParamT())])
    pipe.fit(X=None, y=None, clf__should_succeed=True)
    # classifier should return True
    assert pipe.predict(None)
    # and transformer params should not be changed
    assert pipe.named_steps['transf'].a is None
    assert pipe.named_steps['transf'].b is None
    # invalid parameters should raise an error message
    assert_raise_message(
        TypeError,
        "fit() got an unexpected keyword argument 'bad'",
        pipe.fit, None, None, clf__bad=True
    )


def test_pipeline_sample_weight_supported():
    # Pipeline should pass sample_weight
    X = np.array([[1, 2]])
    pipe = Pipeline([('transf', Transf()), ('clf', FitParamT())])
    pipe.fit(X, y=None)
    assert_equal(pipe.score(X), 3)
    assert_equal(pipe.score(X, y=None), 3)
    assert_equal(pipe.score(X, y=None, sample_weight=None), 3)
    assert_equal(pipe.score(X, sample_weight=np.array([2, 3])), 8)


def test_pipeline_sample_weight_unsupported():
    # When sample_weight is None it shouldn't be passed
    X = np.array([[1, 2]])
    pipe = Pipeline([('transf', Transf()), ('clf', Mult())])
    pipe.fit(X, y=None)
    assert_equal(pipe.score(X), 3)
    assert_equal(pipe.score(X, sample_weight=None), 3)
    assert_raise_message(
        TypeError,
        "score() got an unexpected keyword argument 'sample_weight'",
        pipe.score, X, sample_weight=np.array([2, 3])
    )


def test_pipeline_raise_set_params_error():
    # Test pipeline raises set params error message for nested models.
    pipe = Pipeline([('cls', LinearRegression())])

    # expected error message
    error_msg = ('Invalid parameter %s for estimator %s. '
                 'Check the list of available parameters '
                 'with `estimator.get_params().keys()`.')

    assert_raise_message(ValueError,
                         error_msg % ('fake', pipe),
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
    clf = SVC(gamma='scale', probability=True, random_state=0)
    pca = PCA(svd_solver='full', n_components='mle', whiten=True)
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
    pca = PCA(n_components=2, svd_solver='randomized', whiten=True)
    clf = SVC(gamma='scale', probability=True, random_state=0,
              decision_function_shape='ovr')

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
    # As pipeline doesn't clone estimators on construction,
    # it must have its own estimators
    scaler_for_pipeline = StandardScaler()
    km_for_pipeline = KMeans(random_state=0)

    # first compute the transform and clustering step separately
    scaled = scaler.fit_transform(iris.data)
    separate_pred = km.fit_predict(scaled)

    # use a pipeline to do the transform and clustering in one step
    pipe = Pipeline([
        ('scaler', scaler_for_pipeline),
        ('Kmeans', km_for_pipeline)
    ])
    pipeline_pred = pipe.fit_predict(iris.data)

    assert_array_almost_equal(pipeline_pred, separate_pred)


def test_fit_predict_on_pipeline_without_fit_predict():
    # tests that a pipeline does not have fit_predict method when final
    # step of pipeline does not have fit_predict defined
    scaler = StandardScaler()
    pca = PCA(svd_solver='full')
    pipe = Pipeline([('scaler', scaler), ('pca', pca)])
    assert_raises_regex(AttributeError,
                        "'PCA' object has no attribute 'fit_predict'",
                        getattr, pipe, 'fit_predict')


def test_fit_predict_with_intermediate_fit_params():
    # tests that Pipeline passes fit_params to intermediate steps
    # when fit_predict is invoked
    pipe = Pipeline([('transf', TransfFitParams()), ('clf', FitParamT())])
    pipe.fit_predict(X=None,
                     y=None,
                     transf__should_get_this=True,
                     clf__should_succeed=True)
    assert pipe.named_steps['transf'].fit_params['should_get_this']
    assert pipe.named_steps['clf'].successful
    assert 'should_succeed' not in pipe.named_steps['transf'].fit_params


def test_predict_with_predict_params():
    # tests that Pipeline passes predict_params to the final estimator
    # when predict is invoked
    pipe = Pipeline([('transf', Transf()), ('clf', DummyEstimatorParams())])
    pipe.fit(None, None)
    pipe.predict(X=None, got_attribute=True)

    assert pipe.named_steps['clf'].got_attribute


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

    # Test clone
    fs2 = assert_no_warnings(clone, fs)
    assert fs.transformer_list[0][1] is not fs2.transformer_list[0][1]

    # test setting parameters
    fs.set_params(select__k=2)
    assert_equal(fs.fit_transform(X, y).shape, (X.shape[0], 4))

    # test it works with transformers missing fit_transform
    fs = FeatureUnion([("mock", Transf()), ("svd", svd), ("select", select)])
    X_transformed = fs.fit_transform(X, y)
    assert_equal(X_transformed.shape, (X.shape[0], 8))

    # test error if some elements do not support transform
    assert_raises_regex(TypeError,
                        'All estimators should implement fit and '
                        'transform.*\\bNoTrans\\b',
                        FeatureUnion,
                        [("transform", Transf()), ("no_transform", NoTrans())])

    # test that init accepts tuples
    fs = FeatureUnion((("svd", svd), ("select", select)))
    fs.fit(X, y)


def test_make_union():
    pca = PCA(svd_solver='full')
    mock = Transf()
    fu = make_union(pca, mock)
    names, transformers = zip(*fu.transformer_list)
    assert_equal(names, ("pca", "transf"))
    assert_equal(transformers, (pca, mock))


def test_make_union_kwargs():
    pca = PCA(svd_solver='full')
    mock = Transf()
    fu = make_union(pca, mock, n_jobs=3)
    assert_equal(fu.transformer_list, make_union(pca, mock).transformer_list)
    assert_equal(3, fu.n_jobs)
    # invalid keyword parameters should raise an error message
    assert_raise_message(
        TypeError,
        'Unknown keyword arguments: "transformer_weights"',
        make_union, pca, mock, transformer_weights={'pca': 10, 'Transf': 1}
    )


def test_pipeline_transform():
    # Test whether pipeline works with a transformer at the end.
    # Also test pipeline.transform and pipeline.inverse_transform
    iris = load_iris()
    X = iris.data
    pca = PCA(n_components=2, svd_solver='full')
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
    transf = Transf()
    pipeline = Pipeline([('mock', transf)])

    # test fit_transform:
    X_trans = pipeline.fit_transform(X, y)
    X_trans2 = transf.fit(X, y).transform(X)
    assert_array_almost_equal(X_trans, X_trans2)


def test_set_pipeline_steps():
    transf1 = Transf()
    transf2 = Transf()
    pipeline = Pipeline([('mock', transf1)])
    assert pipeline.named_steps['mock'] is transf1

    # Directly setting attr
    pipeline.steps = [('mock2', transf2)]
    assert 'mock' not in pipeline.named_steps
    assert pipeline.named_steps['mock2'] is transf2
    assert_equal([('mock2', transf2)], pipeline.steps)

    # Using set_params
    pipeline.set_params(steps=[('mock', transf1)])
    assert_equal([('mock', transf1)], pipeline.steps)

    # Using set_params to replace single step
    pipeline.set_params(mock=transf2)
    assert_equal([('mock', transf2)], pipeline.steps)

    # With invalid data
    pipeline.set_params(steps=[('junk', ())])
    assert_raises(TypeError, pipeline.fit, [[1]], [1])
    assert_raises(TypeError, pipeline.fit_transform, [[1]], [1])


def test_pipeline_named_steps():
    transf = Transf()
    mult2 = Mult(mult=2)
    pipeline = Pipeline([('mock', transf), ("mult", mult2)])

    # Test access via named_steps bunch object
    assert 'mock' in pipeline.named_steps
    assert 'mock2' not in pipeline.named_steps
    assert pipeline.named_steps.mock is transf
    assert pipeline.named_steps.mult is mult2

    # Test bunch with conflict attribute of dict
    pipeline = Pipeline([('values', transf), ("mult", mult2)])
    assert pipeline.named_steps.values is not transf
    assert pipeline.named_steps.mult is mult2


@pytest.mark.parametrize('passthrough', [None, 'passthrough'])
def test_pipeline_correctly_adjusts_steps(passthrough):
    X = np.array([[1]])
    y = np.array([1])
    mult2 = Mult(mult=2)
    mult3 = Mult(mult=3)
    mult5 = Mult(mult=5)

    pipeline = Pipeline([
        ('m2', mult2),
        ('bad', passthrough),
        ('m3', mult3),
        ('m5', mult5)
    ])

    pipeline.fit(X, y)
    expected_names = ['m2', 'bad', 'm3', 'm5']
    actual_names = [name for name, _ in pipeline.steps]
    assert expected_names == actual_names


@pytest.mark.parametrize('passthrough', [None, 'passthrough'])
def test_set_pipeline_step_passthrough(passthrough):
    X = np.array([[1]])
    y = np.array([1])
    mult2 = Mult(mult=2)
    mult3 = Mult(mult=3)
    mult5 = Mult(mult=5)

    def make():
        return Pipeline([('m2', mult2), ('m3', mult3), ('last', mult5)])

    pipeline = make()

    exp = 2 * 3 * 5
    assert_array_equal([[exp]], pipeline.fit_transform(X, y))
    assert_array_equal([exp], pipeline.fit(X).predict(X))
    assert_array_equal(X, pipeline.inverse_transform([[exp]]))

    pipeline.set_params(m3=passthrough)
    exp = 2 * 5
    assert_array_equal([[exp]], pipeline.fit_transform(X, y))
    assert_array_equal([exp], pipeline.fit(X).predict(X))
    assert_array_equal(X, pipeline.inverse_transform([[exp]]))
    assert_dict_equal(pipeline.get_params(deep=True),
                      {'steps': pipeline.steps,
                       'm2': mult2,
                       'm3': passthrough,
                       'last': mult5,
                       'memory': None,
                       'm2__mult': 2,
                       'last__mult': 5,
                       })

    pipeline.set_params(m2=passthrough)
    exp = 5
    assert_array_equal([[exp]], pipeline.fit_transform(X, y))
    assert_array_equal([exp], pipeline.fit(X).predict(X))
    assert_array_equal(X, pipeline.inverse_transform([[exp]]))

    # for other methods, ensure no AttributeErrors on None:
    other_methods = ['predict_proba', 'predict_log_proba',
                     'decision_function', 'transform', 'score']
    for method in other_methods:
        getattr(pipeline, method)(X)

    pipeline.set_params(m2=mult2)
    exp = 2 * 5
    assert_array_equal([[exp]], pipeline.fit_transform(X, y))
    assert_array_equal([exp], pipeline.fit(X).predict(X))
    assert_array_equal(X, pipeline.inverse_transform([[exp]]))

    pipeline = make()
    pipeline.set_params(last=passthrough)
    # mult2 and mult3 are active
    exp = 6
    assert_array_equal([[exp]], pipeline.fit(X, y).transform(X))
    assert_array_equal([[exp]], pipeline.fit_transform(X, y))
    assert_array_equal(X, pipeline.inverse_transform([[exp]]))
    assert_raise_message(AttributeError,
                         "'str' object has no attribute 'predict'",
                         getattr, pipeline, 'predict')

    # Check 'passthrough' step at construction time
    exp = 2 * 5
    pipeline = Pipeline(
        [('m2', mult2), ('m3', passthrough), ('last', mult5)])
    assert_array_equal([[exp]], pipeline.fit_transform(X, y))
    assert_array_equal([exp], pipeline.fit(X).predict(X))
    assert_array_equal(X, pipeline.inverse_transform([[exp]]))


def test_pipeline_ducktyping():
    pipeline = make_pipeline(Mult(5))
    pipeline.predict
    pipeline.transform
    pipeline.inverse_transform

    pipeline = make_pipeline(Transf())
    assert not hasattr(pipeline, 'predict')
    pipeline.transform
    pipeline.inverse_transform

    pipeline = make_pipeline('passthrough')
    assert pipeline.steps[0] == ('passthrough', 'passthrough')
    assert not hasattr(pipeline, 'predict')
    pipeline.transform
    pipeline.inverse_transform

    pipeline = make_pipeline(Transf(), NoInvTransf())
    assert not hasattr(pipeline, 'predict')
    pipeline.transform
    assert not hasattr(pipeline, 'inverse_transform')

    pipeline = make_pipeline(NoInvTransf(), Transf())
    assert not hasattr(pipeline, 'predict')
    pipeline.transform
    assert not hasattr(pipeline, 'inverse_transform')


def test_make_pipeline():
    t1 = Transf()
    t2 = Transf()
    pipe = make_pipeline(t1, t2)
    assert isinstance(pipe, Pipeline)
    assert_equal(pipe.steps[0][0], "transf-1")
    assert_equal(pipe.steps[1][0], "transf-2")

    pipe = make_pipeline(t1, t2, FitParamT())
    assert isinstance(pipe, Pipeline)
    assert_equal(pipe.steps[0][0], "transf-1")
    assert_equal(pipe.steps[1][0], "transf-2")
    assert_equal(pipe.steps[2][0], "fitparamt")

    assert_raise_message(
        TypeError,
        'Unknown keyword arguments: "random_parameter"',
        make_pipeline, t1, t2, random_parameter='rnd'
    )


def test_feature_union_weights():
    # test feature union with transformer weights
    iris = load_iris()
    X = iris.data
    y = iris.target
    pca = PCA(n_components=2, svd_solver='randomized', random_state=0)
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
    fs = FeatureUnion([("mock", Transf()), ("pca", pca), ("select", select)],
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
        assert "chars__" in feat or "words__" in feat
    assert_equal(len(feature_names), 35)

    ft = FeatureUnion([("tr1", Transf())]).fit([[1]])
    assert_raise_message(AttributeError,
                         'Transformer tr1 (type Transf) does not provide '
                         'get_feature_names', ft.get_feature_names)


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
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


def test_set_feature_union_steps():
    mult2 = Mult(2)
    mult2.get_feature_names = lambda: ['x2']
    mult3 = Mult(3)
    mult3.get_feature_names = lambda: ['x3']
    mult5 = Mult(5)
    mult5.get_feature_names = lambda: ['x5']

    ft = FeatureUnion([('m2', mult2), ('m3', mult3)])
    assert_array_equal([[2, 3]], ft.transform(np.asarray([[1]])))
    assert_equal(['m2__x2', 'm3__x3'], ft.get_feature_names())

    # Directly setting attr
    ft.transformer_list = [('m5', mult5)]
    assert_array_equal([[5]], ft.transform(np.asarray([[1]])))
    assert_equal(['m5__x5'], ft.get_feature_names())

    # Using set_params
    ft.set_params(transformer_list=[('mock', mult3)])
    assert_array_equal([[3]], ft.transform(np.asarray([[1]])))
    assert_equal(['mock__x3'], ft.get_feature_names())

    # Using set_params to replace single step
    ft.set_params(mock=mult5)
    assert_array_equal([[5]], ft.transform(np.asarray([[1]])))
    assert_equal(['mock__x5'], ft.get_feature_names())


@pytest.mark.parametrize('drop', ['drop', None])
def test_set_feature_union_step_drop(drop):
    mult2 = Mult(2)
    mult2.get_feature_names = lambda: ['x2']
    mult3 = Mult(3)
    mult3.get_feature_names = lambda: ['x3']
    X = np.asarray([[1]])

    ft = FeatureUnion([('m2', mult2), ('m3', mult3)])
    assert_array_equal([[2, 3]], ft.fit(X).transform(X))
    assert_array_equal([[2, 3]], ft.fit_transform(X))
    assert_equal(['m2__x2', 'm3__x3'], ft.get_feature_names())

    ft.set_params(m2=drop)
    assert_array_equal([[3]], ft.fit(X).transform(X))
    assert_array_equal([[3]], ft.fit_transform(X))
    assert_equal(['m3__x3'], ft.get_feature_names())

    ft.set_params(m3=drop)
    assert_array_equal([[]], ft.fit(X).transform(X))
    assert_array_equal([[]], ft.fit_transform(X))
    assert_equal([], ft.get_feature_names())

    # check we can change back
    ft.set_params(m3=mult3)
    assert_array_equal([[3]], ft.fit(X).transform(X))

    # Check 'drop' step at construction time
    ft = FeatureUnion([('m2', drop), ('m3', mult3)])
    assert_array_equal([[3]], ft.fit(X).transform(X))
    assert_array_equal([[3]], ft.fit_transform(X))
    assert_equal(['m3__x3'], ft.get_feature_names())


def test_step_name_validation():
    bad_steps1 = [('a__q', Mult(2)), ('b', Mult(3))]
    bad_steps2 = [('a', Mult(2)), ('a', Mult(3))]
    for cls, param in [(Pipeline, 'steps'),
                       (FeatureUnion, 'transformer_list')]:
        # we validate in construction (despite scikit-learn convention)
        bad_steps3 = [('a', Mult(2)), (param, Mult(3))]
        for bad_steps, message in [
            (bad_steps1, "Estimator names must not contain __: got ['a__q']"),
            (bad_steps2, "Names provided are not unique: ['a', 'a']"),
            (bad_steps3, "Estimator names conflict with constructor "
                         "arguments: ['%s']" % param),
        ]:
            # three ways to make invalid:
            # - construction
            assert_raise_message(ValueError, message, cls,
                                 **{param: bad_steps})

            # - setattr
            est = cls(**{param: [('a', Mult(1))]})
            setattr(est, param, bad_steps)
            assert_raise_message(ValueError, message, est.fit, [[1]], [1])
            assert_raise_message(ValueError, message, est.fit_transform,
                                 [[1]], [1])

            # - set_params
            est = cls(**{param: [('a', Mult(1))]})
            est.set_params(**{param: bad_steps})
            assert_raise_message(ValueError, message, est.fit, [[1]], [1])
            assert_raise_message(ValueError, message, est.fit_transform,
                                 [[1]], [1])


@pytest.mark.filterwarnings('ignore: Default solver will be changed')  # 0.22
@pytest.mark.filterwarnings('ignore: Default multi_class will')  # 0.22
def test_set_params_nested_pipeline():
    estimator = Pipeline([
        ('a', Pipeline([
            ('b', DummyRegressor())
        ]))
    ])
    estimator.set_params(a__b__alpha=0.001, a__b=Lasso())
    estimator.set_params(a__steps=[('b', LogisticRegression())], a__b__C=5)


def test_pipeline_wrong_memory():
    # Test that an error is raised when memory is not a string or a Memory
    # instance
    iris = load_iris()
    X = iris.data
    y = iris.target
    # Define memory as an integer
    memory = 1
    cached_pipe = Pipeline([('transf', DummyTransf()),
                            ('svc', SVC())], memory=memory)
    assert_raises_regex(ValueError, "'memory' should be None, a string or"
                        " have the same interface as joblib.Memory."
                        " Got memory='1' instead.", cached_pipe.fit, X, y)


class DummyMemory(object):
    def cache(self, func):
        return func


class WrongDummyMemory(object):
    pass


def test_pipeline_with_cache_attribute():
    X = np.array([[1, 2]])
    pipe = Pipeline([('transf', Transf()), ('clf', Mult())],
                    memory=DummyMemory())
    pipe.fit(X, y=None)
    dummy = WrongDummyMemory()
    pipe = Pipeline([('transf', Transf()), ('clf', Mult())],
                    memory=dummy)
    assert_raises_regex(ValueError, "'memory' should be None, a string or"
                        " have the same interface as joblib.Memory."
                        " Got memory='{}' instead.".format(dummy), pipe.fit, X)


def test_pipeline_memory():
    iris = load_iris()
    X = iris.data
    y = iris.target
    cachedir = mkdtemp()
    try:
        if LooseVersion(joblib_version) < LooseVersion('0.12'):
            # Deal with change of API in joblib
            memory = Memory(cachedir=cachedir, verbose=10)
        else:
            memory = Memory(location=cachedir, verbose=10)
        # Test with Transformer + SVC
        clf = SVC(gamma='scale', probability=True, random_state=0)
        transf = DummyTransf()
        pipe = Pipeline([('transf', clone(transf)), ('svc', clf)])
        cached_pipe = Pipeline([('transf', transf), ('svc', clf)],
                               memory=memory)

        # Memoize the transformer at the first fit
        cached_pipe.fit(X, y)
        pipe.fit(X, y)
        # Get the time stamp of the transformer in the cached pipeline
        ts = cached_pipe.named_steps['transf'].timestamp_
        # Check that cached_pipe and pipe yield identical results
        assert_array_equal(pipe.predict(X), cached_pipe.predict(X))
        assert_array_equal(pipe.predict_proba(X), cached_pipe.predict_proba(X))
        assert_array_equal(pipe.predict_log_proba(X),
                           cached_pipe.predict_log_proba(X))
        assert_array_equal(pipe.score(X, y), cached_pipe.score(X, y))
        assert_array_equal(pipe.named_steps['transf'].means_,
                           cached_pipe.named_steps['transf'].means_)
        assert not hasattr(transf, 'means_')
        # Check that we are reading the cache while fitting
        # a second time
        cached_pipe.fit(X, y)
        # Check that cached_pipe and pipe yield identical results
        assert_array_equal(pipe.predict(X), cached_pipe.predict(X))
        assert_array_equal(pipe.predict_proba(X), cached_pipe.predict_proba(X))
        assert_array_equal(pipe.predict_log_proba(X),
                           cached_pipe.predict_log_proba(X))
        assert_array_equal(pipe.score(X, y), cached_pipe.score(X, y))
        assert_array_equal(pipe.named_steps['transf'].means_,
                           cached_pipe.named_steps['transf'].means_)
        assert_equal(ts, cached_pipe.named_steps['transf'].timestamp_)
        # Create a new pipeline with cloned estimators
        # Check that even changing the name step does not affect the cache hit
        clf_2 = SVC(gamma='scale', probability=True, random_state=0)
        transf_2 = DummyTransf()
        cached_pipe_2 = Pipeline([('transf_2', transf_2), ('svc', clf_2)],
                                 memory=memory)
        cached_pipe_2.fit(X, y)

        # Check that cached_pipe and pipe yield identical results
        assert_array_equal(pipe.predict(X), cached_pipe_2.predict(X))
        assert_array_equal(pipe.predict_proba(X),
                           cached_pipe_2.predict_proba(X))
        assert_array_equal(pipe.predict_log_proba(X),
                           cached_pipe_2.predict_log_proba(X))
        assert_array_equal(pipe.score(X, y), cached_pipe_2.score(X, y))
        assert_array_equal(pipe.named_steps['transf'].means_,
                           cached_pipe_2.named_steps['transf_2'].means_)
        assert_equal(ts, cached_pipe_2.named_steps['transf_2'].timestamp_)
    finally:
        shutil.rmtree(cachedir)


def test_make_pipeline_memory():
    cachedir = mkdtemp()
    if LooseVersion(joblib_version) < LooseVersion('0.12'):
        # Deal with change of API in joblib
        memory = Memory(cachedir=cachedir, verbose=10)
    else:
        memory = Memory(location=cachedir, verbose=10)
    pipeline = make_pipeline(DummyTransf(), SVC(), memory=memory)
    assert pipeline.memory is memory
    pipeline = make_pipeline(DummyTransf(), SVC())
    assert pipeline.memory is None

    shutil.rmtree(cachedir)
