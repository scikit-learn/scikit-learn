# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Raghav RV <rvraghav93@gmail.com>
# License: BSD 3 clause

import inspect
import warnings
import importlib

from pkgutil import walk_packages
from inspect import signature

import numpy as np

import sklearn
from sklearn.utils import IS_PYPY
from sklearn.utils._testing import check_docstring_parameters
from sklearn.utils._testing import _get_func_name
from sklearn.utils._testing import ignore_warnings
from sklearn.utils import all_estimators
from sklearn.utils.estimator_checks import _enforce_estimator_tags_y
from sklearn.utils.estimator_checks import _enforce_estimator_tags_x
from sklearn.utils.estimator_checks import _construct_instance
from sklearn.utils.deprecation import _is_deprecated
from sklearn.externals._pep562 import Pep562
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression

import pytest


# walk_packages() ignores DeprecationWarnings, now we need to ignore
# FutureWarnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', FutureWarning)
    PUBLIC_MODULES = set([
        pckg[1] for pckg in walk_packages(
            prefix='sklearn.',
            # mypy error: Module has no attribute "__path__"
            path=sklearn.__path__)  # type: ignore  # mypy issue #1422
        if not ("._" in pckg[1] or ".tests." in pckg[1])
    ])

# functions to ignore args / docstring of
_DOCSTRING_IGNORES = [
    'sklearn.utils.deprecation.load_mlcomp',
    'sklearn.pipeline.make_pipeline',
    'sklearn.pipeline.make_union',
    'sklearn.utils.extmath.safe_sparse_dot',
    'sklearn.utils._joblib'
]

# Methods where y param should be ignored if y=None by default
_METHODS_IGNORE_NONE_Y = [
    'fit',
    'score',
    'fit_predict',
    'fit_transform',
    'partial_fit',
    'predict'
]


# numpydoc 0.8.0's docscrape tool raises because of collections.abc under
# Python 3.7
@pytest.mark.filterwarnings('ignore::FutureWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.skipif(IS_PYPY, reason='test segfaults on PyPy')
def test_docstring_parameters():
    # Test module docstring formatting

    # Skip test if numpydoc is not found
    pytest.importorskip('numpydoc',
                        reason="numpydoc is required to test the docstrings")

    # XXX unreached code as of v0.22
    from numpydoc import docscrape

    incorrect = []
    for name in PUBLIC_MODULES:
        if name == 'sklearn.utils.fixes':
            # We cannot always control these docstrings
            continue
        with warnings.catch_warnings(record=True):
            module = importlib.import_module(name)
        classes = inspect.getmembers(module, inspect.isclass)
        # Exclude non-scikit-learn classes
        classes = [cls for cls in classes
                   if cls[1].__module__.startswith('sklearn')]
        for cname, cls in classes:
            this_incorrect = []
            if cname in _DOCSTRING_IGNORES or cname.startswith('_'):
                continue
            if inspect.isabstract(cls):
                continue
            with warnings.catch_warnings(record=True) as w:
                cdoc = docscrape.ClassDoc(cls)
            if len(w):
                raise RuntimeError('Error for __init__ of %s in %s:\n%s'
                                   % (cls, name, w[0]))

            cls_init = getattr(cls, '__init__', None)

            if _is_deprecated(cls_init):
                continue
            elif cls_init is not None:
                this_incorrect += check_docstring_parameters(
                    cls.__init__, cdoc)

            for method_name in cdoc.methods:
                method = getattr(cls, method_name)
                if _is_deprecated(method):
                    continue
                param_ignore = None
                # Now skip docstring test for y when y is None
                # by default for API reason
                if method_name in _METHODS_IGNORE_NONE_Y:
                    sig = signature(method)
                    if ('y' in sig.parameters and
                            sig.parameters['y'].default is None):
                        param_ignore = ['y']  # ignore y for fit and score
                result = check_docstring_parameters(
                    method, ignore=param_ignore)
                this_incorrect += result

            incorrect += this_incorrect

        functions = inspect.getmembers(module, inspect.isfunction)
        # Exclude imported functions
        functions = [fn for fn in functions if fn[1].__module__ == name]
        for fname, func in functions:
            # Don't test private methods / functions
            if fname.startswith('_'):
                continue
            if fname == "configuration" and name.endswith("setup"):
                continue
            name_ = _get_func_name(func)
            if (not any(d in name_ for d in _DOCSTRING_IGNORES) and
                    not _is_deprecated(func)):
                incorrect += check_docstring_parameters(func)

    msg = '\n'.join(incorrect)
    if len(incorrect) > 0:
        raise AssertionError("Docstring Error:\n" + msg)


@ignore_warnings(category=FutureWarning)
def test_tabs():
    # Test that there are no tabs in our source files
    for importer, modname, ispkg in walk_packages(sklearn.__path__,
                                                  prefix='sklearn.'):

        if IS_PYPY and ('_svmlight_format_io' in modname or
                        'feature_extraction._hashing_fast' in modname):
            continue

        # because we don't import
        mod = importlib.import_module(modname)

        # TODO: Remove when minimum python version is 3.7
        # unwrap to get module because Pep562 backport wraps the original
        # module
        if isinstance(mod, Pep562):
            mod = mod._module

        try:
            source = inspect.getsource(mod)
        except IOError:  # user probably should have run "make clean"
            continue
        assert '\t' not in source, ('"%s" has tabs, please remove them ',
                                    'or add it to the ignore list'
                                    % modname)


def _construct_searchcv_instance(SearchCV):
    return SearchCV(LogisticRegression(), {"C": [0.1, 1]})


N_FEATURES_MODULES_TO_IGNORE = {
    'calibration',
    'cluster',
    'compose',
    'covariance',
    'decomposition',
    'discriminant_analysis',
    'dummy',
    'ensemble',
    'feature_extraction',
    'feature_selection',
    'gaussian_process',
    'impute',
    'isotonic',
    'kernel_approximation',
    'kernel_ridge',
    'linear_model',
    'manifold',
    'mixture',
    'model_selection',
    'multiclass',
    'multioutput',
    'naive_bayes',
    'neighbors',
    'neural_network',
    'pipeline',
    'preprocessing',
    'random_projection',
    'semi_supervised',
    'svm',
    'tree'
}


@pytest.mark.parametrize('name, Estimator',
                         all_estimators())
def test_fit_docstring_attributes(name, Estimator):
    pytest.importorskip('numpydoc')
    from numpydoc import docscrape

    doc = docscrape.ClassDoc(Estimator)
    attributes = doc['Attributes']

    IGNORED = {'ClassifierChain', 'ColumnTransformer',
               'CountVectorizer', 'DictVectorizer', 'FeatureUnion',
               'GaussianRandomProjection',
               'MultiOutputClassifier', 'MultiOutputRegressor',
               'NoSampleWeightWrapper', 'OneVsOneClassifier',
               'OutputCodeClassifier', 'Pipeline', 'RFE', 'RFECV',
               'RegressorChain', 'SelectFromModel',
               'SparseCoder', 'SparseRandomProjection',
               'SpectralBiclustering', 'StackingClassifier',
               'StackingRegressor', 'TfidfVectorizer', 'VotingClassifier',
               'VotingRegressor', 'SequentialFeatureSelector',
               'HalvingGridSearchCV', 'HalvingRandomSearchCV'}
    if Estimator.__name__ in IGNORED or Estimator.__name__.startswith('_'):
        pytest.skip("Estimator cannot be fit easily to test fit attributes")

    if Estimator.__name__ in ("RandomizedSearchCV", "GridSearchCV"):
        est = _construct_searchcv_instance(Estimator)
    else:
        est = _construct_instance(Estimator)

    if Estimator.__name__ == 'SelectKBest':
        est.k = 2

    if Estimator.__name__ == 'DummyClassifier':
        est.strategy = "stratified"

    if 'PLS' in Estimator.__name__ or 'CCA' in Estimator.__name__:
        est.n_components = 1  # default = 2 is invalid for single target.

    # FIXME: TO BE REMOVED for 1.0 (avoid FutureWarning)
    if Estimator.__name__ == 'AffinityPropagation':
        est.random_state = 63

    # FIXME: TO BE REMOVED for 1.1 (avoid FutureWarning)
    if Estimator.__name__ == 'NMF':
        est.init = 'nndsvda'

    X, y = make_classification(n_samples=20, n_features=3,
                               n_redundant=0, n_classes=2,
                               random_state=2)

    y = _enforce_estimator_tags_y(est, y)
    X = _enforce_estimator_tags_x(est, X)

    if '1dlabels' in est._get_tags()['X_types']:
        est.fit(y)
    elif '2dlabels' in est._get_tags()['X_types']:
        est.fit(np.c_[y, y])
    else:
        est.fit(X, y)

    skipped_attributes = {'x_scores_',  # For PLS, TODO remove in 1.1
                          'y_scores_'}  # For PLS, TODO remove in 1.1

    module = est.__module__.split(".")[1]
    if module in N_FEATURES_MODULES_TO_IGNORE:
        skipped_attributes.add("n_features_in_")

    for attr in attributes:
        if attr.name in skipped_attributes:
            continue
        desc = ' '.join(attr.desc).lower()
        # As certain attributes are present "only" if a certain parameter is
        # provided, this checks if the word "only" is present in the attribute
        # description, and if not the attribute is required to be present.
        if 'only ' in desc:
            continue
        # ignore deprecation warnings
        with ignore_warnings(category=FutureWarning):
            assert hasattr(est, attr.name)

    IGNORED = {'Birch', 'LarsCV', 'Lasso',
               'OrthogonalMatchingPursuit'}

    if Estimator.__name__ in IGNORED:
        pytest.xfail(
            reason="Estimator has too many undocumented attributes.")

    fit_attr = [k for k in est.__dict__.keys() if k.endswith('_')
                and not k.startswith('_')]
    fit_attr_names = [attr.name for attr in attributes]
    undocumented_attrs = set(fit_attr).difference(fit_attr_names)
    undocumented_attrs = set(undocumented_attrs).difference(skipped_attributes)
    assert not undocumented_attrs,\
        "Undocumented attributes: {}".format(undocumented_attrs)
