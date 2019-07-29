import types
import warnings
import sys
import traceback
import pickle
from copy import deepcopy
from functools import partial
from inspect import signature

import numpy as np
from scipy import sparse
from scipy.stats import rankdata

from . import IS_PYPY
from . import _joblib
from .testing import assert_raises, _get_args
from .testing import assert_raises_regex
from .testing import assert_raise_message
from .testing import assert_equal
from .testing import assert_not_equal
from .testing import assert_in
from .testing import assert_array_equal
from .testing import assert_array_almost_equal
from .testing import assert_allclose
from .testing import assert_allclose_dense_sparse
from .testing import assert_warns_message
from .testing import set_random_state
from .testing import assert_greater
from .testing import assert_greater_equal
from .testing import SkipTest
from .testing import ignore_warnings
from .testing import assert_dict_equal
from .testing import create_memmap_backed_data
from . import is_scalar_nan
from ..discriminant_analysis import LinearDiscriminantAnalysis
from ..linear_model import Ridge


from ..base import (clone, ClusterMixin, is_classifier, is_regressor,
                          _DEFAULT_TAGS, RegressorMixin, is_outlier_detector)

from ..metrics import accuracy_score, adjusted_rand_score, f1_score

from ..random_projection import BaseRandomProjection
from ..feature_selection import SelectKBest
from ..pipeline import make_pipeline
from ..exceptions import DataConversionWarning
from ..exceptions import SkipTestWarning
from ..model_selection import train_test_split
from ..model_selection import ShuffleSplit
from ..model_selection._validation import _safe_split
from ..metrics.pairwise import (rbf_kernel, linear_kernel,
                                      pairwise_distances)

from .import shuffle
from .validation import has_fit_parameter, _num_samples
from ..preprocessing import StandardScaler
from ..datasets import load_iris, load_boston, make_blobs


BOSTON = None
CROSS_DECOMPOSITION = ['PLSCanonical', 'PLSRegression', 'CCA', 'PLSSVD']


def _safe_tags(estimator, key=None):
    # if estimator doesn't have _get_tags, use _DEFAULT_TAGS
    # if estimator has tags but not key, use _DEFAULT_TAGS[key]
    if hasattr(estimator, "_get_tags"):
        if key is not None:
            return estimator._get_tags().get(key, _DEFAULT_TAGS[key])
        tags = estimator._get_tags()
        return {key: tags.get(key, _DEFAULT_TAGS[key])
                for key in _DEFAULT_TAGS.keys()}
    if key is not None:
        return _DEFAULT_TAGS[key]
    return _DEFAULT_TAGS


def _yield_checks(name, estimator):
    tags = _safe_tags(estimator)
    yield check_estimators_dtypes
    yield check_fit_score_takes_y
    yield check_sample_weights_pandas_series
    yield check_sample_weights_list
    yield check_sample_weights_invariance
    yield check_estimators_fit_returns_self
    yield partial(check_estimators_fit_returns_self, readonly_memmap=True)

    # Check that all estimator yield informative messages when
    # trained on empty datasets
    if not tags["no_validation"]:
        yield check_complex_data
        yield check_dtype_object
        yield check_estimators_empty_data_messages

    if name not in CROSS_DECOMPOSITION:
        # cross-decomposition's "transform" returns X and Y
        yield check_pipeline_consistency

    if not tags["allow_nan"] and not tags["no_validation"]:
        # Test that all estimators check their input for NaN's and infs
        yield check_estimators_nan_inf

    yield check_estimators_overwrite_params
    if hasattr(estimator, 'sparsify'):
        yield check_sparsify_coefficients

    yield check_estimator_sparse_data

    # Test that estimators can be pickled, and once pickled
    # give the same answer as before.
    yield check_estimators_pickle


def _yield_classifier_checks(name, classifier):
    tags = _safe_tags(classifier)

    # test classifiers can handle non-array data
    yield check_classifier_data_not_an_array
    # test classifiers trained on a single label always return this label
    yield check_classifiers_one_label
    yield check_classifiers_classes
    yield check_estimators_partial_fit_n_features
    # basic consistency testing
    yield check_classifiers_train
    yield partial(check_classifiers_train, readonly_memmap=True)
    yield check_classifiers_regression_target
    if not tags["no_validation"]:
        yield check_supervised_y_no_nan
        yield check_supervised_y_2d
    yield check_estimators_unfitted
    if 'class_weight' in classifier.get_params().keys():
        yield check_class_weight_classifiers

    yield check_non_transformer_estimators_n_iter
    # test if predict_proba is a monotonic transformation of decision_function
    yield check_decision_proba_consistency


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_supervised_y_no_nan(name, estimator_orig):
    # Checks that the Estimator targets are not NaN.
    estimator = clone(estimator_orig)
    rng = np.random.RandomState(888)
    X = rng.randn(10, 5)
    y = np.full(10, np.inf)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    errmsg = "Input contains NaN, infinity or a value too large for " \
             "dtype('float64')."
    try:
        estimator.fit(X, y)
    except ValueError as e:
        if str(e) != errmsg:
            raise ValueError("Estimator {0} raised error as expected, but "
                             "does not match expected error message"
                             .format(name))
    else:
        raise ValueError("Estimator {0} should have raised error on fitting "
                         "array y with NaN value.".format(name))


def _yield_regressor_checks(name, regressor):
    tags = _safe_tags(regressor)
    # TODO: test with intercept
    # TODO: test with multiple responses
    # basic testing
    yield check_regressors_train
    yield partial(check_regressors_train, readonly_memmap=True)
    yield check_regressor_data_not_an_array
    yield check_estimators_partial_fit_n_features
    yield check_regressors_no_decision_function
    if not tags["no_validation"]:
        yield check_supervised_y_2d
    yield check_supervised_y_no_nan
    if name != 'CCA':
        # check that the regressor handles int input
        yield check_regressors_int
    yield check_estimators_unfitted
    yield check_non_transformer_estimators_n_iter


def _yield_transformer_checks(name, transformer):
    # All transformers should either deal with sparse data or raise an
    # exception with type TypeError and an intelligible error message
    yield check_transformer_data_not_an_array
    # these don't actually fit the data, so don't raise errors
    yield check_transformer_general
    yield partial(check_transformer_general, readonly_memmap=True)

    if not _safe_tags(transformer, "stateless"):
        yield check_transformers_unfitted
    # Dependent on external solvers and hence accessing the iter
    # param is non-trivial.
    external_solver = ['Isomap', 'KernelPCA', 'LocallyLinearEmbedding',
                       'RandomizedLasso', 'LogisticRegressionCV']
    if name not in external_solver:
        yield check_transformer_n_iter


def _yield_clustering_checks(name, clusterer):
    yield check_clusterer_compute_labels_predict
    if name not in ('WardAgglomeration', "FeatureAgglomeration"):
        # this is clustering on the features
        # let's not test that here.
        yield check_clustering
        yield partial(check_clustering, readonly_memmap=True)
        yield check_estimators_partial_fit_n_features
    yield check_non_transformer_estimators_n_iter


def _yield_outliers_checks(name, estimator):

    # checks for outlier detectors that have a fit_predict method
    if hasattr(estimator, 'fit_predict'):
        yield check_outliers_fit_predict

    # checks for estimators that can be used on a test set
    if hasattr(estimator, 'predict'):
        yield check_outliers_train
        yield partial(check_outliers_train, readonly_memmap=True)
        # test outlier detectors can handle non-array data
        yield check_classifier_data_not_an_array
        # test if NotFittedError is raised
        yield check_estimators_unfitted


def _yield_all_checks(name, estimator):
    tags = _safe_tags(estimator)
    if "2darray" not in tags["X_types"]:
        warnings.warn("Can't test estimator {} which requires input "
                      " of type {}".format(name, tags["X_types"]),
                      SkipTestWarning)
        return
    if tags["_skip_test"]:
        warnings.warn("Explicit SKIP via _skip_test tag for estimator "
                      "{}.".format(name),
                      SkipTestWarning)
        return

    for check in _yield_checks(name, estimator):
        yield check
    if is_classifier(estimator):
        for check in _yield_classifier_checks(name, estimator):
            yield check
    if is_regressor(estimator):
        for check in _yield_regressor_checks(name, estimator):
            yield check
    if hasattr(estimator, 'transform'):
        for check in _yield_transformer_checks(name, estimator):
            yield check
    if isinstance(estimator, ClusterMixin):
        for check in _yield_clustering_checks(name, estimator):
            yield check
    if is_outlier_detector(estimator):
        for check in _yield_outliers_checks(name, estimator):
            yield check
    yield check_fit2d_predict1d
    yield check_methods_subset_invariance
    yield check_fit2d_1sample
    yield check_fit2d_1feature
    yield check_fit1d
    yield check_get_params_invariance
    yield check_set_params
    yield check_dict_unchanged
    yield check_dont_overwrite_parameters
    yield check_fit_idempotent


def check_estimator(Estimator):
    """Check if estimator adheres to scikit-learn conventions.

    This estimator will run an extensive test-suite for input validation,
    shapes, etc.
    Additional tests for classifiers, regressors, clustering or transformers
    will be run if the Estimator class inherits from the corresponding mixin
    from sklearn.base.

    This test can be applied to classes or instances.
    Classes currently have some additional tests that related to construction,
    while passing instances allows the testing of multiple options.

    Parameters
    ----------
    estimator : estimator object or class
        Estimator to check. Estimator is a class object or instance.

    """
    if isinstance(Estimator, type):
        # got a class
        name = Estimator.__name__
        estimator = Estimator()
        check_parameters_default_constructible(name, Estimator)
        check_no_attributes_set_in_init(name, estimator)
    else:
        # got an instance
        estimator = Estimator
        name = type(estimator).__name__

    for check in _yield_all_checks(name, estimator):
        try:
            check(name, estimator)
        except SkipTest as exception:
            # the only SkipTest thrown currently results from not
            # being able to import pandas.
            warnings.warn(str(exception), SkipTestWarning)


def _boston_subset(n_samples=200):
    global BOSTON
    if BOSTON is None:
        boston = load_boston()
        X, y = boston.data, boston.target
        X, y = shuffle(X, y, random_state=0)
        X, y = X[:n_samples], y[:n_samples]
        X = StandardScaler().fit_transform(X)
        BOSTON = X, y
    return BOSTON


def set_checking_parameters(estimator):
    # set parameters to speed up some estimators and
    # avoid deprecated behaviour
    params = estimator.get_params()
    name = estimator.__class__.__name__
    if ("n_iter" in params and name != "TSNE"):
        estimator.set_params(n_iter=5)
    if "max_iter" in params:
        if estimator.max_iter is not None:
            estimator.set_params(max_iter=min(5, estimator.max_iter))
        # LinearSVR, LinearSVC
        if estimator.__class__.__name__ in ['LinearSVR', 'LinearSVC']:
            estimator.set_params(max_iter=20)
        # NMF
        if estimator.__class__.__name__ == 'NMF':
            estimator.set_params(max_iter=100)
        # MLP
        if estimator.__class__.__name__ in ['MLPClassifier', 'MLPRegressor']:
            estimator.set_params(max_iter=100)
    if "n_resampling" in params:
        # randomized lasso
        estimator.set_params(n_resampling=5)
    if "n_estimators" in params:
        # especially gradient boosting with default 100
        # FIXME: The default number of trees was changed and is set to 'warn'
        # for some of the ensemble methods. We need to catch this case to avoid
        # an error during the comparison. To be reverted in 0.22.
        if estimator.n_estimators == 'warn':
            estimator.set_params(n_estimators=5)
        else:
            estimator.set_params(n_estimators=min(5, estimator.n_estimators))
    if "max_trials" in params:
        # RANSAC
        estimator.set_params(max_trials=10)
    if "n_init" in params:
        # K-Means
        estimator.set_params(n_init=2)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 2

    if name == 'TruncatedSVD':
        # TruncatedSVD doesn't run with n_components = n_features
        # This is ugly :-/
        estimator.n_components = 1

    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = min(estimator.n_clusters, 2)

    if hasattr(estimator, "n_best"):
        estimator.n_best = 1

    if name == "SelectFdr":
        # be tolerant of noisy datasets (not actually speed)
        estimator.set_params(alpha=.5)

    if name == "TheilSenRegressor":
        estimator.max_subpopulation = 100

    if estimator.__class__.__name__ == "IsolationForest":
        # XXX to be removed in 0.22.
        # this is used because the old IsolationForest does not
        # respect the outlier detection API and thus and does not
        # pass the outlier detection common tests.
        estimator.set_params(behaviour='new')

    if isinstance(estimator, BaseRandomProjection):
        # Due to the jl lemma and often very few samples, the number
        # of components of the random matrix projection will be probably
        # greater than the number of features.
        # So we impose a smaller number (avoid "auto" mode)
        estimator.set_params(n_components=2)

    if isinstance(estimator, SelectKBest):
        # SelectKBest has a default of k=10
        # which is more feature than we have in most case.
        estimator.set_params(k=1)

    if name in ('HistGradientBoostingClassifier',
                'HistGradientBoostingRegressor'):
        # The default min_samples_leaf (20) isn't appropriate for small
        # datasets (only very shallow trees are built) that the checks use.
        estimator.set_params(min_samples_leaf=5)


class NotAnArray:
    """An object that is convertible to an array

    Parameters
    ----------
    data : array_like
        The data.
    """

    def __init__(self, data):
        self.data = data

    def __array__(self, dtype=None):
        return self.data


def _is_pairwise(estimator):
    """Returns True if estimator has a _pairwise attribute set to True.

    Parameters
    ----------
    estimator : object
        Estimator object to test.

    Returns
    -------
    out : bool
        True if _pairwise is set to True and False otherwise.
    """
    return bool(getattr(estimator, "_pairwise", False))


def _is_pairwise_metric(estimator):
    """Returns True if estimator accepts pairwise metric.

    Parameters
    ----------
    estimator : object
        Estimator object to test.

    Returns
    -------
    out : bool
        True if _pairwise is set to True and False otherwise.
    """
    metric = getattr(estimator, "metric", None)

    return bool(metric == 'precomputed')


def pairwise_estimator_convert_X(X, estimator, kernel=linear_kernel):

    if _is_pairwise_metric(estimator):
        return pairwise_distances(X, metric='euclidean')
    if _is_pairwise(estimator):
        return kernel(X, X)

    return X


def _generate_sparse_matrix(X_csr):
    """Generate sparse matrices with {32,64}bit indices of diverse format

        Parameters
        ----------
        X_csr: CSR Matrix
            Input matrix in CSR format

        Returns
        -------
        out: iter(Matrices)
            In format['dok', 'lil', 'dia', 'bsr', 'csr', 'csc', 'coo',
             'coo_64', 'csc_64', 'csr_64']
    """

    assert X_csr.format == 'csr'
    yield 'csr', X_csr.copy()
    for sparse_format in ['dok', 'lil', 'dia', 'bsr', 'csc', 'coo']:
        yield sparse_format, X_csr.asformat(sparse_format)

    # Generate large indices matrix only if its supported by scipy
    X_coo = X_csr.asformat('coo')
    X_coo.row = X_coo.row.astype('int64')
    X_coo.col = X_coo.col.astype('int64')
    yield "coo_64", X_coo

    for sparse_format in ['csc', 'csr']:
        X = X_csr.asformat(sparse_format)
        X.indices = X.indices.astype('int64')
        X.indptr = X.indptr.astype('int64')
        yield sparse_format + "_64", X


def check_estimator_sparse_data(name, estimator_orig):

    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X = pairwise_estimator_convert_X(X, estimator_orig)
    X_csr = sparse.csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)
    # catch deprecation warnings
    with ignore_warnings(category=DeprecationWarning):
        estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)
    for matrix_format, X in _generate_sparse_matrix(X_csr):
        # catch deprecation warnings
        with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
            estimator = clone(estimator_orig)
            if name in ['Scaler', 'StandardScaler']:
                estimator.set_params(with_mean=False)
        # fit and predict
        try:
            with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
                estimator.fit(X, y)
            if hasattr(estimator, "predict"):
                pred = estimator.predict(X)
                if _safe_tags(estimator, "multioutput_only"):
                    assert_equal(pred.shape, (X.shape[0], 1))
                else:
                    assert_equal(pred.shape, (X.shape[0],))
            if hasattr(estimator, 'predict_proba'):
                probs = estimator.predict_proba(X)
                assert_equal(probs.shape, (X.shape[0], 4))
        except (TypeError, ValueError) as e:
            if 'sparse' not in repr(e).lower():
                if "64" in matrix_format:
                    msg = ("Estimator %s doesn't seem to support %s matrix, "
                           "and is not failing gracefully, e.g. by using "
                           "check_array(X, accept_large_sparse=False)")
                    raise AssertionError(msg % (name, matrix_format))
                else:
                    print("Estimator %s doesn't seem to fail gracefully on "
                          "sparse data: error message state explicitly that "
                          "sparse input is not supported if this is not"
                          " the case." % name)
                    raise
        except Exception:
            print("Estimator %s doesn't seem to fail gracefully on "
                  "sparse data: it should raise a TypeError if sparse input "
                  "is explicitly not supported." % name)
            raise


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_sample_weights_pandas_series(name, estimator_orig):
    # check that estimators will accept a 'sample_weight' parameter of
    # type pandas.Series in the 'fit' function.
    estimator = clone(estimator_orig)
    if has_fit_parameter(estimator, "sample_weight"):
        try:
            import pandas as pd
            X = np.array([[1, 1], [1, 2], [1, 3], [1, 4],
                          [2, 1], [2, 2], [2, 3], [2, 4]])
            X = pd.DataFrame(pairwise_estimator_convert_X(X, estimator_orig))
            y = pd.Series([1, 1, 1, 1, 2, 2, 2, 2])
            weights = pd.Series([1] * 8)
            if _safe_tags(estimator, "multioutput_only"):
                y = pd.DataFrame(y)
            try:
                estimator.fit(X, y, sample_weight=weights)
            except ValueError:
                raise ValueError("Estimator {0} raises error if "
                                 "'sample_weight' parameter is of "
                                 "type pandas.Series".format(name))
        except ImportError:
            raise SkipTest("pandas is not installed: not testing for "
                           "input of type pandas.Series to class weight.")


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_sample_weights_list(name, estimator_orig):
    # check that estimators will accept a 'sample_weight' parameter of
    # type list in the 'fit' function.
    if has_fit_parameter(estimator_orig, "sample_weight"):
        estimator = clone(estimator_orig)
        rnd = np.random.RandomState(0)
        X = pairwise_estimator_convert_X(rnd.uniform(size=(10, 3)),
                                         estimator_orig)
        y = np.arange(10) % 3
        y = multioutput_estimator_convert_y_2d(estimator, y)
        sample_weight = [3] * 10
        # Test that estimators don't raise any exception
        estimator.fit(X, y, sample_weight=sample_weight)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_sample_weights_invariance(name, estimator_orig):
    # check that the estimators yield same results for
    # unit weights and no weights
    if (has_fit_parameter(estimator_orig, "sample_weight") and
            not (hasattr(estimator_orig, "_pairwise")
                 and estimator_orig._pairwise)):
        # We skip pairwise because the data is not pairwise

        estimator1 = clone(estimator_orig)
        estimator2 = clone(estimator_orig)
        set_random_state(estimator1, random_state=0)
        set_random_state(estimator2, random_state=0)

        X = np.array([[1, 3], [1, 3], [1, 3], [1, 3],
                      [2, 1], [2, 1], [2, 1], [2, 1],
                      [3, 3], [3, 3], [3, 3], [3, 3],
                      [4, 1], [4, 1], [4, 1], [4, 1]], dtype=np.dtype('float'))
        y = np.array([1, 1, 1, 1, 2, 2, 2, 2,
                      1, 1, 1, 1, 2, 2, 2, 2], dtype=np.dtype('int'))
        y = multioutput_estimator_convert_y_2d(estimator1, y)

        estimator1.fit(X, y=y, sample_weight=np.ones(shape=len(y)))
        estimator2.fit(X, y=y, sample_weight=None)

        for method in ["predict", "transform"]:
            if hasattr(estimator_orig, method):
                X_pred1 = getattr(estimator1, method)(X)
                X_pred2 = getattr(estimator2, method)(X)
                if sparse.issparse(X_pred1):
                    X_pred1 = X_pred1.toarray()
                    X_pred2 = X_pred2.toarray()
                assert_allclose(X_pred1, X_pred2,
                                err_msg="For %s sample_weight=None is not"
                                        " equivalent to sample_weight=ones"
                                        % name)


@ignore_warnings(category=(DeprecationWarning, FutureWarning, UserWarning))
def check_dtype_object(name, estimator_orig):
    # check that estimators treat dtype object as numeric if possible
    rng = np.random.RandomState(0)
    X = pairwise_estimator_convert_X(rng.rand(40, 10), estimator_orig)
    X = X.astype(object)
    y = (X[:, 0] * 4).astype(np.int)
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    estimator.fit(X, y)
    if hasattr(estimator, "predict"):
        estimator.predict(X)

    if hasattr(estimator, "transform"):
        estimator.transform(X)

    try:
        estimator.fit(X, y.astype(object))
    except Exception as e:
        if "Unknown label type" not in str(e):
            raise

    tags = _safe_tags(estimator)
    if 'string' not in tags['X_types']:
        X[0, 0] = {'foo': 'bar'}
        msg = "argument must be a string.* number"
        assert_raises_regex(TypeError, msg, estimator.fit, X, y)
    else:
        # Estimators supporting string will not call np.asarray to convert the
        # data to numeric and therefore, the error will not be raised.
        # Checking for each element dtype in the input array will be costly.
        # Refer to #11401 for full discussion.
        estimator.fit(X, y)


def check_complex_data(name, estimator_orig):
    # check that estimators raise an exception on providing complex data
    X = np.random.sample(10) + 1j * np.random.sample(10)
    X = X.reshape(-1, 1)
    y = np.random.sample(10) + 1j * np.random.sample(10)
    estimator = clone(estimator_orig)
    assert_raises_regex(ValueError, "Complex data not supported",
                        estimator.fit, X, y)


@ignore_warnings
def check_dict_unchanged(name, estimator_orig):
    # this estimator raises
    # ValueError: Found array with 0 feature(s) (shape=(23, 0))
    # while a minimum of 1 is required.
    # error
    if name in ['SpectralCoclustering']:
        return
    rnd = np.random.RandomState(0)
    if name in ['RANSACRegressor']:
        X = 3 * rnd.uniform(size=(20, 3))
    else:
        X = 2 * rnd.uniform(size=(20, 3))

    X = pairwise_estimator_convert_X(X, estimator_orig)

    y = X[:, 0].astype(np.int)
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)
    if hasattr(estimator, "n_components"):
        estimator.n_components = 1

    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    if hasattr(estimator, "n_best"):
        estimator.n_best = 1

    set_random_state(estimator, 1)

    estimator.fit(X, y)
    for method in ["predict", "transform", "decision_function",
                   "predict_proba"]:
        if hasattr(estimator, method):
            dict_before = estimator.__dict__.copy()
            getattr(estimator, method)(X)
            assert_dict_equal(estimator.__dict__, dict_before,
                              'Estimator changes __dict__ during %s' % method)


def is_public_parameter(attr):
    return not (attr.startswith('_') or attr.endswith('_'))


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_dont_overwrite_parameters(name, estimator_orig):
    # check that fit method only changes or sets private attributes
    if hasattr(estimator_orig.__init__, "deprecated_original"):
        # to not check deprecated classes
        return
    estimator = clone(estimator_orig)
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    X = pairwise_estimator_convert_X(X, estimator_orig)
    y = X[:, 0].astype(np.int)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    dict_before_fit = estimator.__dict__.copy()
    estimator.fit(X, y)

    dict_after_fit = estimator.__dict__

    public_keys_after_fit = [key for key in dict_after_fit.keys()
                             if is_public_parameter(key)]

    attrs_added_by_fit = [key for key in public_keys_after_fit
                          if key not in dict_before_fit.keys()]

    # check that fit doesn't add any public attribute
    assert not attrs_added_by_fit, (
            'Estimator adds public attribute(s) during'
            ' the fit method.'
            ' Estimators are only allowed to add private attributes'
            ' either started with _ or ended'
            ' with _ but %s added'
            % ', '.join(attrs_added_by_fit))

    # check that fit doesn't change any public attribute
    attrs_changed_by_fit = [key for key in public_keys_after_fit
                            if (dict_before_fit[key]
                                is not dict_after_fit[key])]

    assert not attrs_changed_by_fit, (
            'Estimator changes public attribute(s) during'
            ' the fit method. Estimators are only allowed'
            ' to change attributes started'
            ' or ended with _, but'
            ' %s changed'
            % ', '.join(attrs_changed_by_fit))


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_fit2d_predict1d(name, estimator_orig):
    # check by fitting a 2d array and predicting with a 1d array
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    X = pairwise_estimator_convert_X(X, estimator_orig)
    y = X[:, 0].astype(np.int)
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    estimator.fit(X, y)
    tags = _safe_tags(estimator)
    if tags["no_validation"]:
        # FIXME this is a bit loose
        return

    for method in ["predict", "transform", "decision_function",
                   "predict_proba"]:
        if hasattr(estimator, method):
            assert_raise_message(ValueError, "Reshape your data",
                                 getattr(estimator, method), X[0])


def _apply_on_subsets(func, X):
    # apply function on the whole set and on mini batches
    result_full = func(X)
    n_features = X.shape[1]
    result_by_batch = [func(batch.reshape(1, n_features))
                       for batch in X]
    # func can output tuple (e.g. score_samples)
    if type(result_full) == tuple:
        result_full = result_full[0]
        result_by_batch = list(map(lambda x: x[0], result_by_batch))

    if sparse.issparse(result_full):
        result_full = result_full.A
        result_by_batch = [x.A for x in result_by_batch]
    return np.ravel(result_full), np.ravel(result_by_batch)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_methods_subset_invariance(name, estimator_orig):
    # check that method gives invariant results if applied
    # on mini bathes or the whole set
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20, 3))
    X = pairwise_estimator_convert_X(X, estimator_orig)
    y = X[:, 0].astype(np.int)
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    estimator.fit(X, y)

    for method in ["predict", "transform", "decision_function",
                   "score_samples", "predict_proba"]:

        msg = ("{method} of {name} is not invariant when applied "
               "to a subset.").format(method=method, name=name)
        # TODO remove cases when corrected
        if (name, method) in [('NuSVC', 'decision_function'),
                              ('SparsePCA', 'transform'),
                              ('MiniBatchSparsePCA', 'transform'),
                              ('DummyClassifier', 'predict'),
                              ('BernoulliRBM', 'score_samples')]:
            raise SkipTest(msg)

        if hasattr(estimator, method):
            result_full, result_by_batch = _apply_on_subsets(
                getattr(estimator, method), X)
            assert_allclose(result_full, result_by_batch,
                            atol=1e-7, err_msg=msg)


@ignore_warnings
def check_fit2d_1sample(name, estimator_orig):
    # Check that fitting a 2d array with only one sample either works or
    # returns an informative message. The error message should either mention
    # the number of samples or the number of classes.
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(1, 10))
    y = X[:, 0].astype(np.int)
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)

    # min_cluster_size cannot be less than the data size for OPTICS.
    if name == 'OPTICS':
        estimator.set_params(min_samples=1)

    msgs = ["1 sample", "n_samples = 1", "n_samples=1", "one sample",
            "1 class", "one class"]

    try:
        estimator.fit(X, y)
    except ValueError as e:
        if all(msg not in repr(e) for msg in msgs):
            raise e


@ignore_warnings
def check_fit2d_1feature(name, estimator_orig):
    # check fitting a 2d array with only 1 feature either works or returns
    # informative message
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(10, 1))
    X = pairwise_estimator_convert_X(X, estimator_orig)
    y = X[:, 0].astype(np.int)
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1
    # ensure two labels in subsample for RandomizedLogisticRegression
    if name == 'RandomizedLogisticRegression':
        estimator.sample_fraction = 1
    # ensure non skipped trials for RANSACRegressor
    if name == 'RANSACRegressor':
        estimator.residual_threshold = 0.5

    y = multioutput_estimator_convert_y_2d(estimator, y)
    set_random_state(estimator, 1)

    msgs = ["1 feature(s)", "n_features = 1", "n_features=1"]

    try:
        estimator.fit(X, y)
    except ValueError as e:
        if all(msg not in repr(e) for msg in msgs):
            raise e


@ignore_warnings
def check_fit1d(name, estimator_orig):
    # check fitting 1d X array raises a ValueError
    rnd = np.random.RandomState(0)
    X = 3 * rnd.uniform(size=(20))
    y = X.astype(np.int)
    estimator = clone(estimator_orig)
    tags = _safe_tags(estimator)
    if tags["no_validation"]:
        # FIXME this is a bit loose
        return
    y = multioutput_estimator_convert_y_2d(estimator, y)

    if hasattr(estimator, "n_components"):
        estimator.n_components = 1
    if hasattr(estimator, "n_clusters"):
        estimator.n_clusters = 1

    set_random_state(estimator, 1)
    assert_raises(ValueError, estimator.fit, X, y)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_transformer_general(name, transformer, readonly_memmap=False):
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    X = StandardScaler().fit_transform(X)
    X -= X.min()

    if readonly_memmap:
        X, y = create_memmap_backed_data([X, y])

    _check_transformer(name, transformer, X, y)
    _check_transformer(name, transformer, X.tolist(), y.tolist())


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_transformer_data_not_an_array(name, transformer):
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    X = StandardScaler().fit_transform(X)
    # We need to make sure that we have non negative data, for things
    # like NMF
    X -= X.min() - .1
    this_X = NotAnArray(X)
    this_y = NotAnArray(np.asarray(y))
    _check_transformer(name, transformer, this_X, this_y)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_transformers_unfitted(name, transformer):
    X, y = _boston_subset()

    transformer = clone(transformer)
    with assert_raises((AttributeError, ValueError), msg="The unfitted "
                       "transformer {} does not raise an error when "
                       "transform is called. Perhaps use "
                       "check_is_fitted in transform.".format(name)):
        transformer.transform(X)


def _check_transformer(name, transformer_orig, X, y):
    n_samples, n_features = np.asarray(X).shape
    transformer = clone(transformer_orig)
    set_random_state(transformer)

    # fit

    if name in CROSS_DECOMPOSITION:
        y_ = np.c_[y, y]
        y_[::2, 1] *= 2
    else:
        y_ = y

    transformer.fit(X, y_)
    # fit_transform method should work on non fitted estimator
    transformer_clone = clone(transformer)
    X_pred = transformer_clone.fit_transform(X, y=y_)

    if isinstance(X_pred, tuple):
        for x_pred in X_pred:
            assert_equal(x_pred.shape[0], n_samples)
    else:
        # check for consistent n_samples
        assert_equal(X_pred.shape[0], n_samples)

    if hasattr(transformer, 'transform'):
        if name in CROSS_DECOMPOSITION:
            X_pred2 = transformer.transform(X, y_)
            X_pred3 = transformer.fit_transform(X, y=y_)
        else:
            X_pred2 = transformer.transform(X)
            X_pred3 = transformer.fit_transform(X, y=y_)

        if _safe_tags(transformer_orig, 'non_deterministic'):
            msg = name + ' is non deterministic'
            raise SkipTest(msg)
        if isinstance(X_pred, tuple) and isinstance(X_pred2, tuple):
            for x_pred, x_pred2, x_pred3 in zip(X_pred, X_pred2, X_pred3):
                assert_allclose_dense_sparse(
                    x_pred, x_pred2, atol=1e-2,
                    err_msg="fit_transform and transform outcomes "
                            "not consistent in %s"
                    % transformer)
                assert_allclose_dense_sparse(
                    x_pred, x_pred3, atol=1e-2,
                    err_msg="consecutive fit_transform outcomes "
                            "not consistent in %s"
                    % transformer)
        else:
            assert_allclose_dense_sparse(
                X_pred, X_pred2,
                err_msg="fit_transform and transform outcomes "
                        "not consistent in %s"
                % transformer, atol=1e-2)
            assert_allclose_dense_sparse(
                X_pred, X_pred3, atol=1e-2,
                err_msg="consecutive fit_transform outcomes "
                        "not consistent in %s"
                % transformer)
            assert_equal(_num_samples(X_pred2), n_samples)
            assert_equal(_num_samples(X_pred3), n_samples)

        # raises error on malformed input for transform
        if hasattr(X, 'T') and not _safe_tags(transformer, "stateless"):
            # If it's not an array, it does not have a 'T' property
            with assert_raises(ValueError, msg="The transformer {} does "
                               "not raise an error when the number of "
                               "features in transform is different from"
                               " the number of features in "
                               "fit.".format(name)):
                transformer.transform(X.T)


@ignore_warnings
def check_pipeline_consistency(name, estimator_orig):
    if _safe_tags(estimator_orig, 'non_deterministic'):
        msg = name + ' is non deterministic'
        raise SkipTest(msg)

    # check that make_pipeline(est) gives same score as est
    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)
    X -= X.min()
    X = pairwise_estimator_convert_X(X, estimator_orig, kernel=rbf_kernel)
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)
    set_random_state(estimator)
    pipeline = make_pipeline(estimator)
    estimator.fit(X, y)
    pipeline.fit(X, y)

    funcs = ["score", "fit_transform"]

    for func_name in funcs:
        func = getattr(estimator, func_name, None)
        if func is not None:
            func_pipeline = getattr(pipeline, func_name)
            result = func(X, y)
            result_pipe = func_pipeline(X, y)
            assert_allclose_dense_sparse(result, result_pipe)


@ignore_warnings
def check_fit_score_takes_y(name, estimator_orig):
    # check that all estimators accept an optional y
    # in fit and score so they can be used in pipelines
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 3))
    X = pairwise_estimator_convert_X(X, estimator_orig)
    y = np.arange(10) % 3
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)
    set_random_state(estimator)

    funcs = ["fit", "score", "partial_fit", "fit_predict", "fit_transform"]
    for func_name in funcs:
        func = getattr(estimator, func_name, None)
        if func is not None:
            func(X, y)
            args = [p.name for p in signature(func).parameters.values()]
            if args[0] == "self":
                # if_delegate_has_method makes methods into functions
                # with an explicit "self", so need to shift arguments
                args = args[1:]
            assert args[1] in ["y", "Y"], (
                    "Expected y or Y as second argument for method "
                    "%s of %s. Got arguments: %r."
                    % (func_name, type(estimator).__name__, args))


@ignore_warnings
def check_estimators_dtypes(name, estimator_orig):
    rnd = np.random.RandomState(0)
    X_train_32 = 3 * rnd.uniform(size=(20, 5)).astype(np.float32)
    X_train_32 = pairwise_estimator_convert_X(X_train_32, estimator_orig)
    X_train_64 = X_train_32.astype(np.float64)
    X_train_int_64 = X_train_32.astype(np.int64)
    X_train_int_32 = X_train_32.astype(np.int32)
    y = X_train_int_64[:, 0]
    y = multioutput_estimator_convert_y_2d(estimator_orig, y)

    methods = ["predict", "transform", "decision_function", "predict_proba"]

    for X_train in [X_train_32, X_train_64, X_train_int_64, X_train_int_32]:
        estimator = clone(estimator_orig)
        set_random_state(estimator, 1)
        estimator.fit(X_train, y)

        for method in methods:
            if hasattr(estimator, method):
                getattr(estimator, method)(X_train)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_estimators_empty_data_messages(name, estimator_orig):
    e = clone(estimator_orig)
    set_random_state(e, 1)

    X_zero_samples = np.empty(0).reshape(0, 3)
    # The precise message can change depending on whether X or y is
    # validated first. Let us test the type of exception only:
    with assert_raises(ValueError, msg="The estimator {} does not"
                       " raise an error when an empty data is used "
                       "to train. Perhaps use "
                       "check_array in train.".format(name)):
        e.fit(X_zero_samples, [])

    X_zero_features = np.empty(0).reshape(3, 0)
    # the following y should be accepted by both classifiers and regressors
    # and ignored by unsupervised models
    y = multioutput_estimator_convert_y_2d(e, np.array([1, 0, 1]))
    msg = (r"0 feature\(s\) \(shape=\(3, 0\)\) while a minimum of \d* "
           "is required.")
    assert_raises_regex(ValueError, msg, e.fit, X_zero_features, y)


@ignore_warnings(category=DeprecationWarning)
def check_estimators_nan_inf(name, estimator_orig):
    # Checks that Estimator X's do not contain NaN or inf.
    rnd = np.random.RandomState(0)
    X_train_finite = pairwise_estimator_convert_X(rnd.uniform(size=(10, 3)),
                                                  estimator_orig)
    X_train_nan = rnd.uniform(size=(10, 3))
    X_train_nan[0, 0] = np.nan
    X_train_inf = rnd.uniform(size=(10, 3))
    X_train_inf[0, 0] = np.inf
    y = np.ones(10)
    y[:5] = 0
    y = multioutput_estimator_convert_y_2d(estimator_orig, y)
    error_string_fit = "Estimator doesn't check for NaN and inf in fit."
    error_string_predict = ("Estimator doesn't check for NaN and inf in"
                            " predict.")
    error_string_transform = ("Estimator doesn't check for NaN and inf in"
                              " transform.")
    for X_train in [X_train_nan, X_train_inf]:
        # catch deprecation warnings
        with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
            estimator = clone(estimator_orig)
            set_random_state(estimator, 1)
            # try to fit
            try:
                estimator.fit(X_train, y)
            except ValueError as e:
                if 'inf' not in repr(e) and 'NaN' not in repr(e):
                    print(error_string_fit, estimator, e)
                    traceback.print_exc(file=sys.stdout)
                    raise e
            except Exception as exc:
                print(error_string_fit, estimator, exc)
                traceback.print_exc(file=sys.stdout)
                raise exc
            else:
                raise AssertionError(error_string_fit, estimator)
            # actually fit
            estimator.fit(X_train_finite, y)

            # predict
            if hasattr(estimator, "predict"):
                try:
                    estimator.predict(X_train)
                except ValueError as e:
                    if 'inf' not in repr(e) and 'NaN' not in repr(e):
                        print(error_string_predict, estimator, e)
                        traceback.print_exc(file=sys.stdout)
                        raise e
                except Exception as exc:
                    print(error_string_predict, estimator, exc)
                    traceback.print_exc(file=sys.stdout)
                else:
                    raise AssertionError(error_string_predict, estimator)

            # transform
            if hasattr(estimator, "transform"):
                try:
                    estimator.transform(X_train)
                except ValueError as e:
                    if 'inf' not in repr(e) and 'NaN' not in repr(e):
                        print(error_string_transform, estimator, e)
                        traceback.print_exc(file=sys.stdout)
                        raise e
                except Exception as exc:
                    print(error_string_transform, estimator, exc)
                    traceback.print_exc(file=sys.stdout)
                else:
                    raise AssertionError(error_string_transform, estimator)


@ignore_warnings
def check_estimators_pickle(name, estimator_orig):
    """Test that we can pickle all estimators"""
    check_methods = ["predict", "transform", "decision_function",
                     "predict_proba"]

    X, y = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                      random_state=0, n_features=2, cluster_std=0.1)

    # some estimators can't do features less than 0
    X -= X.min()
    X = pairwise_estimator_convert_X(X, estimator_orig, kernel=rbf_kernel)

    tags = _safe_tags(estimator_orig)
    # include NaN values when the estimator should deal with them
    if tags['allow_nan']:
        # set randomly 10 elements to np.nan
        rng = np.random.RandomState(42)
        mask = rng.choice(X.size, 10, replace=False)
        X.reshape(-1)[mask] = np.nan

    estimator = clone(estimator_orig)

    # some estimators only take multioutputs
    y = multioutput_estimator_convert_y_2d(estimator, y)

    set_random_state(estimator)
    estimator.fit(X, y)

    result = dict()
    for method in check_methods:
        if hasattr(estimator, method):
            result[method] = getattr(estimator, method)(X)

    # pickle and unpickle!
    pickled_estimator = pickle.dumps(estimator)
    if estimator.__module__.startswith('sklearn.'):
        assert b"version" in pickled_estimator
    unpickled_estimator = pickle.loads(pickled_estimator)

    result = dict()
    for method in check_methods:
        if hasattr(estimator, method):
            result[method] = getattr(estimator, method)(X)

    for method in result:
        unpickled_result = getattr(unpickled_estimator, method)(X)
        assert_allclose_dense_sparse(result[method], unpickled_result)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_estimators_partial_fit_n_features(name, estimator_orig):
    # check if number of features changes between calls to partial_fit.
    if not hasattr(estimator_orig, 'partial_fit'):
        return
    estimator = clone(estimator_orig)
    X, y = make_blobs(n_samples=50, random_state=1)
    X -= X.min()

    try:
        if is_classifier(estimator):
            classes = np.unique(y)
            estimator.partial_fit(X, y, classes=classes)
        else:
            estimator.partial_fit(X, y)
    except NotImplementedError:
        return

    with assert_raises(ValueError,
                       msg="The estimator {} does not raise an"
                           " error when the number of features"
                           " changes between calls to "
                           "partial_fit.".format(name)):
        estimator.partial_fit(X[:, :-1], y)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_clustering(name, clusterer_orig, readonly_memmap=False):
    clusterer = clone(clusterer_orig)
    X, y = make_blobs(n_samples=50, random_state=1)
    X, y = shuffle(X, y, random_state=7)
    X = StandardScaler().fit_transform(X)
    rng = np.random.RandomState(7)
    X_noise = np.concatenate([X, rng.uniform(low=-3, high=3, size=(5, 2))])

    if readonly_memmap:
        X, y, X_noise = create_memmap_backed_data([X, y, X_noise])

    n_samples, n_features = X.shape
    # catch deprecation and neighbors warnings
    if hasattr(clusterer, "n_clusters"):
        clusterer.set_params(n_clusters=3)
    set_random_state(clusterer)
    if name == 'AffinityPropagation':
        clusterer.set_params(preference=-100)
        clusterer.set_params(max_iter=100)

    # fit
    clusterer.fit(X)
    # with lists
    clusterer.fit(X.tolist())

    pred = clusterer.labels_
    assert_equal(pred.shape, (n_samples,))
    assert_greater(adjusted_rand_score(pred, y), 0.4)
    if _safe_tags(clusterer, 'non_deterministic'):
        return
    set_random_state(clusterer)
    with warnings.catch_warnings(record=True):
        pred2 = clusterer.fit_predict(X)
    assert_array_equal(pred, pred2)

    # fit_predict(X) and labels_ should be of type int
    assert_in(pred.dtype, [np.dtype('int32'), np.dtype('int64')])
    assert_in(pred2.dtype, [np.dtype('int32'), np.dtype('int64')])

    # Add noise to X to test the possible values of the labels
    labels = clusterer.fit_predict(X_noise)

    # There should be at least one sample in every cluster. Equivalently
    # labels_ should contain all the consecutive values between its
    # min and its max.
    labels_sorted = np.unique(labels)
    assert_array_equal(labels_sorted, np.arange(labels_sorted[0],
                                                labels_sorted[-1] + 1))

    # Labels are expected to start at 0 (no noise) or -1 (if noise)
    assert labels_sorted[0] in [0, -1]
    # Labels should be less than n_clusters - 1
    if hasattr(clusterer, 'n_clusters'):
        n_clusters = getattr(clusterer, 'n_clusters')
        assert_greater_equal(n_clusters - 1, labels_sorted[-1])
    # else labels should be less than max(labels_) which is necessarily true


@ignore_warnings(category=DeprecationWarning)
def check_clusterer_compute_labels_predict(name, clusterer_orig):
    """Check that predict is invariant of compute_labels"""
    X, y = make_blobs(n_samples=20, random_state=0)
    clusterer = clone(clusterer_orig)
    set_random_state(clusterer)

    if hasattr(clusterer, "compute_labels"):
        # MiniBatchKMeans
        X_pred1 = clusterer.fit(X).predict(X)
        clusterer.set_params(compute_labels=False)
        X_pred2 = clusterer.fit(X).predict(X)
        assert_array_equal(X_pred1, X_pred2)


@ignore_warnings(category=DeprecationWarning)
def check_classifiers_one_label(name, classifier_orig):
    error_string_fit = "Classifier can't train when only one class is present."
    error_string_predict = ("Classifier can't predict when only one class is "
                            "present.")
    rnd = np.random.RandomState(0)
    X_train = rnd.uniform(size=(10, 3))
    X_test = rnd.uniform(size=(10, 3))
    y = np.ones(10)
    # catch deprecation warnings
    with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
        classifier = clone(classifier_orig)
        # try to fit
        try:
            classifier.fit(X_train, y)
        except ValueError as e:
            if 'class' not in repr(e):
                print(error_string_fit, classifier, e)
                traceback.print_exc(file=sys.stdout)
                raise e
            else:
                return
        except Exception as exc:
            print(error_string_fit, classifier, exc)
            traceback.print_exc(file=sys.stdout)
            raise exc
        # predict
        try:
            assert_array_equal(classifier.predict(X_test), y)
        except Exception as exc:
            print(error_string_predict, classifier, exc)
            raise exc


@ignore_warnings  # Warnings are raised by decision function
def check_classifiers_train(name, classifier_orig, readonly_memmap=False):
    X_m, y_m = make_blobs(n_samples=300, random_state=0)
    X_m, y_m = shuffle(X_m, y_m, random_state=7)
    X_m = StandardScaler().fit_transform(X_m)
    # generate binary problem from multi-class one
    y_b = y_m[y_m != 2]
    X_b = X_m[y_m != 2]
    tags = _safe_tags(classifier_orig)

    if name in ['BernoulliNB', 'MultinomialNB', 'ComplementNB']:
        X_m -= X_m.min()
        X_b -= X_b.min()

    if readonly_memmap:
        X_m, y_m, X_b, y_b = create_memmap_backed_data([X_m, y_m, X_b, y_b])

    for (X, y) in [(X_m, y_m), (X_b, y_b)]:
        classes = np.unique(y)
        n_classes = len(classes)
        n_samples, n_features = X.shape
        classifier = clone(classifier_orig)
        X = pairwise_estimator_convert_X(X, classifier)
        y = multioutput_estimator_convert_y_2d(classifier, y)

        set_random_state(classifier)
        # raises error on malformed input for fit
        if not tags["no_validation"]:
            with assert_raises(
                ValueError,
                msg="The classifier {} does not "
                    "raise an error when incorrect/malformed input "
                    "data for fit is passed. The number of training "
                    "examples is not the same as the number of labels. "
                    "Perhaps use check_X_y in fit.".format(name)):
                classifier.fit(X, y[:-1])

        # fit
        classifier.fit(X, y)
        # with lists
        classifier.fit(X.tolist(), y.tolist())
        assert hasattr(classifier, "classes_")
        y_pred = classifier.predict(X)

        assert_equal(y_pred.shape, (n_samples,))
        # training set performance
        if not tags['poor_score']:
            assert_greater(accuracy_score(y, y_pred), 0.83)

        # raises error on malformed input for predict
        msg_pairwise = (
            "The classifier {} does not raise an error when shape of X in "
            " {} is not equal to (n_test_samples, n_training_samples)")
        msg = ("The classifier {} does not raise an error when the number of "
               "features in {} is different from the number of features in "
               "fit.")

        if not tags["no_validation"]:
            if _is_pairwise(classifier):
                with assert_raises(ValueError,
                                   msg=msg_pairwise.format(name, "predict")):
                    classifier.predict(X.reshape(-1, 1))
            else:
                with assert_raises(ValueError,
                                   msg=msg.format(name, "predict")):
                    classifier.predict(X.T)
        if hasattr(classifier, "decision_function"):
            try:
                # decision_function agrees with predict
                decision = classifier.decision_function(X)
                if n_classes == 2:
                    if not tags["multioutput_only"]:
                        assert_equal(decision.shape, (n_samples,))
                    else:
                        assert_equal(decision.shape, (n_samples, 1))
                    dec_pred = (decision.ravel() > 0).astype(np.int)
                    assert_array_equal(dec_pred, y_pred)
                else:
                    assert_equal(decision.shape, (n_samples, n_classes))
                    assert_array_equal(np.argmax(decision, axis=1), y_pred)

                # raises error on malformed input for decision_function
                if not tags["no_validation"]:
                    if _is_pairwise(classifier):
                        with assert_raises(ValueError, msg=msg_pairwise.format(
                                name, "decision_function")):
                            classifier.decision_function(X.reshape(-1, 1))
                    else:
                        with assert_raises(ValueError, msg=msg.format(
                                name, "decision_function")):
                            classifier.decision_function(X.T)
            except NotImplementedError:
                pass

        if hasattr(classifier, "predict_proba"):
            # predict_proba agrees with predict
            y_prob = classifier.predict_proba(X)
            assert_equal(y_prob.shape, (n_samples, n_classes))
            assert_array_equal(np.argmax(y_prob, axis=1), y_pred)
            # check that probas for all classes sum to one
            assert_array_almost_equal(np.sum(y_prob, axis=1),
                                      np.ones(n_samples))
            if not tags["no_validation"]:
                # raises error on malformed input for predict_proba
                if _is_pairwise(classifier_orig):
                    with assert_raises(ValueError, msg=msg_pairwise.format(
                            name, "predict_proba")):
                        classifier.predict_proba(X.reshape(-1, 1))
                else:
                    with assert_raises(ValueError, msg=msg.format(
                            name, "predict_proba")):
                        classifier.predict_proba(X.T)
            if hasattr(classifier, "predict_log_proba"):
                # predict_log_proba is a transformation of predict_proba
                y_log_prob = classifier.predict_log_proba(X)
                assert_allclose(y_log_prob, np.log(y_prob), 8, atol=1e-9)
                assert_array_equal(np.argsort(y_log_prob), np.argsort(y_prob))


def check_outlier_corruption(num_outliers, expected_outliers, decision):
    # Check for deviation from the precise given contamination level that may
    # be due to ties in the anomaly scores.
    if num_outliers < expected_outliers:
        start = num_outliers
        end = expected_outliers + 1
    else:
        start = expected_outliers
        end = num_outliers + 1

    # ensure that all values in the 'critical area' are tied,
    # leading to the observed discrepancy between provided
    # and actual contamination levels.
    sorted_decision = np.sort(decision)
    msg = ('The number of predicted outliers is not equal to the expected '
           'number of outliers and this difference is not explained by the '
           'number of ties in the decision_function values')
    assert len(np.unique(sorted_decision[start:end])) == 1, msg


def check_outliers_train(name, estimator_orig, readonly_memmap=True):
    n_samples = 300
    X, _ = make_blobs(n_samples=n_samples, random_state=0)
    X = shuffle(X, random_state=7)

    if readonly_memmap:
        X = create_memmap_backed_data(X)

    n_samples, n_features = X.shape
    estimator = clone(estimator_orig)
    set_random_state(estimator)

    # fit
    estimator.fit(X)
    # with lists
    estimator.fit(X.tolist())

    y_pred = estimator.predict(X)
    assert y_pred.shape == (n_samples,)
    assert y_pred.dtype.kind == 'i'
    assert_array_equal(np.unique(y_pred), np.array([-1, 1]))

    decision = estimator.decision_function(X)
    scores = estimator.score_samples(X)
    for output in [decision, scores]:
        assert output.dtype == np.dtype('float')
        assert output.shape == (n_samples,)

    # raises error on malformed input for predict
    assert_raises(ValueError, estimator.predict, X.T)

    # decision_function agrees with predict
    dec_pred = (decision >= 0).astype(np.int)
    dec_pred[dec_pred == 0] = -1
    assert_array_equal(dec_pred, y_pred)

    # raises error on malformed input for decision_function
    assert_raises(ValueError, estimator.decision_function, X.T)

    # decision_function is a translation of score_samples
    y_dec = scores - estimator.offset_
    assert_allclose(y_dec, decision)

    # raises error on malformed input for score_samples
    assert_raises(ValueError, estimator.score_samples, X.T)

    # contamination parameter (not for OneClassSVM which has the nu parameter)
    if (hasattr(estimator, 'contamination')
            and not hasattr(estimator, 'novelty')):
        # proportion of outliers equal to contamination parameter when not
        # set to 'auto'. This is true for the training set and cannot thus be
        # checked as follows for estimators with a novelty parameter such as
        # LocalOutlierFactor (tested in check_outliers_fit_predict)
        expected_outliers = 30
        contamination = expected_outliers / n_samples
        estimator.set_params(contamination=contamination)
        estimator.fit(X)
        y_pred = estimator.predict(X)

        num_outliers = np.sum(y_pred != 1)
        # num_outliers should be equal to expected_outliers unless
        # there are ties in the decision_function values. this can
        # only be tested for estimators with a decision_function
        # method, i.e. all estimators except LOF which is already
        # excluded from this if branch.
        if num_outliers != expected_outliers:
            decision = estimator.decision_function(X)
            check_outlier_corruption(num_outliers, expected_outliers, decision)

        # raises error when contamination is a scalar and not in [0,1]
        for contamination in [-0.5, 2.3]:
            estimator.set_params(contamination=contamination)
            assert_raises(ValueError, estimator.fit, X)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_estimators_fit_returns_self(name, estimator_orig,
                                      readonly_memmap=False):
    """Check if self is returned when calling fit"""
    X, y = make_blobs(random_state=0, n_samples=9, n_features=4)
    # some want non-negative input
    X -= X.min()
    X = pairwise_estimator_convert_X(X, estimator_orig)

    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    if readonly_memmap:
        X, y = create_memmap_backed_data([X, y])

    set_random_state(estimator)
    assert estimator.fit(X, y) is estimator


@ignore_warnings
def check_estimators_unfitted(name, estimator_orig):
    """Check that predict raises an exception in an unfitted estimator.

    Unfitted estimators should raise either AttributeError or ValueError.
    The specific exception type NotFittedError inherits from both and can
    therefore be adequately raised for that purpose.
    """

    # Common test for Regressors, Classifiers and Outlier detection estimators
    X, y = _boston_subset()

    estimator = clone(estimator_orig)

    msg = "fit"
    if hasattr(estimator, 'predict'):
        can_predict = False
        try:
            # some models can predict without fitting
            # like GaussianProcess regressors
            # in this case, we skip this test
            pred = estimator.predict(X)
            assert pred.shape[0] == X.shape[0]
            can_predict = True
        except ValueError:
            pass
        if can_predict:
            raise SkipTest(
                "{} can predict without fitting, skipping "
                "check_estimator_unfitted.".format(name))

        assert_raise_message((AttributeError, ValueError), msg,
                             estimator.predict, X)

    if hasattr(estimator, 'decision_function'):
        assert_raise_message((AttributeError, ValueError), msg,
                             estimator.decision_function, X)

    if hasattr(estimator, 'predict_proba'):
        assert_raise_message((AttributeError, ValueError), msg,
                             estimator.predict_proba, X)

    if hasattr(estimator, 'predict_log_proba'):
        assert_raise_message((AttributeError, ValueError), msg,
                             estimator.predict_log_proba, X)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_supervised_y_2d(name, estimator_orig):
    if _safe_tags(estimator_orig, "multioutput_only"):
        # These only work on 2d, so this test makes no sense
        return
    rnd = np.random.RandomState(0)
    X = pairwise_estimator_convert_X(rnd.uniform(size=(10, 3)), estimator_orig)
    y = np.arange(10) % 3
    estimator = clone(estimator_orig)
    set_random_state(estimator)
    # fit
    estimator.fit(X, y)
    y_pred = estimator.predict(X)

    set_random_state(estimator)
    # Check that when a 2D y is given, a DataConversionWarning is
    # raised
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", DataConversionWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        estimator.fit(X, y[:, np.newaxis])
    y_pred_2d = estimator.predict(X)
    msg = "expected 1 DataConversionWarning, got: %s" % (
        ", ".join([str(w_x) for w_x in w]))
    if not _safe_tags(estimator, "multioutput"):
        # check that we warned if we don't support multi-output
        assert_greater(len(w), 0, msg)
        assert "DataConversionWarning('A column-vector y" \
               " was passed when a 1d array was expected" in msg
    assert_allclose(y_pred.ravel(), y_pred_2d.ravel())


@ignore_warnings
def check_classifiers_predictions(X, y, name, classifier_orig):
    classes = np.unique(y)
    classifier = clone(classifier_orig)
    if name == 'BernoulliNB':
        X = X > X.mean()
    set_random_state(classifier)

    classifier.fit(X, y)
    y_pred = classifier.predict(X)

    if hasattr(classifier, "decision_function"):
        decision = classifier.decision_function(X)
        assert isinstance(decision, np.ndarray)
        if len(classes) == 2:
            dec_pred = (decision.ravel() > 0).astype(np.int)
            dec_exp = classifier.classes_[dec_pred]
            assert_array_equal(dec_exp, y_pred,
                               err_msg="decision_function does not match "
                               "classifier for %r: expected '%s', got '%s'" %
                               (classifier, ", ".join(map(str, dec_exp)),
                                ", ".join(map(str, y_pred))))
        elif getattr(classifier, 'decision_function_shape', 'ovr') == 'ovr':
            decision_y = np.argmax(decision, axis=1).astype(int)
            y_exp = classifier.classes_[decision_y]
            assert_array_equal(y_exp, y_pred,
                               err_msg="decision_function does not match "
                               "classifier for %r: expected '%s', got '%s'" %
                               (classifier, ", ".join(map(str, y_exp)),
                                ", ".join(map(str, y_pred))))

    # training set performance
    if name != "ComplementNB":
        # This is a pathological data set for ComplementNB.
        # For some specific cases 'ComplementNB' predicts less classes
        # than expected
        assert_array_equal(np.unique(y), np.unique(y_pred))
    assert_array_equal(classes, classifier.classes_,
                       err_msg="Unexpected classes_ attribute for %r: "
                       "expected '%s', got '%s'" %
                       (classifier, ", ".join(map(str, classes)),
                        ", ".join(map(str, classifier.classes_))))


def choose_check_classifiers_labels(name, y, y_names):
    return y if name in ["LabelPropagation", "LabelSpreading"] else y_names


def check_classifiers_classes(name, classifier_orig):
    X_multiclass, y_multiclass = make_blobs(n_samples=30, random_state=0,
                                            cluster_std=0.1)
    X_multiclass, y_multiclass = shuffle(X_multiclass, y_multiclass,
                                         random_state=7)
    X_multiclass = StandardScaler().fit_transform(X_multiclass)
    # We need to make sure that we have non negative data, for things
    # like NMF
    X_multiclass -= X_multiclass.min() - .1

    X_binary = X_multiclass[y_multiclass != 2]
    y_binary = y_multiclass[y_multiclass != 2]

    X_multiclass = pairwise_estimator_convert_X(X_multiclass, classifier_orig)
    X_binary = pairwise_estimator_convert_X(X_binary, classifier_orig)

    labels_multiclass = ["one", "two", "three"]
    labels_binary = ["one", "two"]

    y_names_multiclass = np.take(labels_multiclass, y_multiclass)
    y_names_binary = np.take(labels_binary, y_binary)

    for X, y, y_names in [(X_multiclass, y_multiclass, y_names_multiclass),
                          (X_binary, y_binary, y_names_binary)]:
        for y_names_i in [y_names, y_names.astype('O')]:
            y_ = choose_check_classifiers_labels(name, y, y_names_i)
            check_classifiers_predictions(X, y_, name, classifier_orig)

    labels_binary = [-1, 1]
    y_names_binary = np.take(labels_binary, y_binary)
    y_binary = choose_check_classifiers_labels(name, y_binary, y_names_binary)
    check_classifiers_predictions(X_binary, y_binary, name, classifier_orig)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_regressors_int(name, regressor_orig):
    X, _ = _boston_subset()
    X = pairwise_estimator_convert_X(X[:50], regressor_orig)
    rnd = np.random.RandomState(0)
    y = rnd.randint(3, size=X.shape[0])
    y = multioutput_estimator_convert_y_2d(regressor_orig, y)
    rnd = np.random.RandomState(0)
    # separate estimators to control random seeds
    regressor_1 = clone(regressor_orig)
    regressor_2 = clone(regressor_orig)
    set_random_state(regressor_1)
    set_random_state(regressor_2)

    if name in CROSS_DECOMPOSITION:
        y_ = np.vstack([y, 2 * y + rnd.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y

    # fit
    regressor_1.fit(X, y_)
    pred1 = regressor_1.predict(X)
    regressor_2.fit(X, y_.astype(np.float))
    pred2 = regressor_2.predict(X)
    assert_allclose(pred1, pred2, atol=1e-2, err_msg=name)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_regressors_train(name, regressor_orig, readonly_memmap=False):
    X, y = _boston_subset()
    X = pairwise_estimator_convert_X(X, regressor_orig)
    y = StandardScaler().fit_transform(y.reshape(-1, 1))  # X is already scaled
    y = y.ravel()
    regressor = clone(regressor_orig)
    y = multioutput_estimator_convert_y_2d(regressor, y)
    if name in CROSS_DECOMPOSITION:
        rnd = np.random.RandomState(0)
        y_ = np.vstack([y, 2 * y + rnd.randint(2, size=len(y))])
        y_ = y_.T
    else:
        y_ = y

    if readonly_memmap:
        X, y, y_ = create_memmap_backed_data([X, y, y_])

    if not hasattr(regressor, 'alphas') and hasattr(regressor, 'alpha'):
        # linear regressors need to set alpha, but not generalized CV ones
        regressor.alpha = 0.01
    if name == 'PassiveAggressiveRegressor':
        regressor.C = 0.01

    # raises error on malformed input for fit
    with assert_raises(ValueError, msg="The classifier {} does not"
                       " raise an error when incorrect/malformed input "
                       "data for fit is passed. The number of training "
                       "examples is not the same as the number of "
                       "labels. Perhaps use check_X_y in fit.".format(name)):
        regressor.fit(X, y[:-1])
    # fit
    set_random_state(regressor)
    regressor.fit(X, y_)
    regressor.fit(X.tolist(), y_.tolist())
    y_pred = regressor.predict(X)
    assert_equal(y_pred.shape, y_.shape)

    # TODO: find out why PLS and CCA fail. RANSAC is random
    # and furthermore assumes the presence of outliers, hence
    # skipped
    if not _safe_tags(regressor, "poor_score"):
        assert_greater(regressor.score(X, y_), 0.5)


@ignore_warnings
def check_regressors_no_decision_function(name, regressor_orig):
    # checks whether regressors have decision_function or predict_proba
    rng = np.random.RandomState(0)
    X = rng.normal(size=(10, 4))
    regressor = clone(regressor_orig)
    y = multioutput_estimator_convert_y_2d(regressor, X[:, 0])

    if hasattr(regressor, "n_components"):
        # FIXME CCA, PLS is not robust to rank 1 effects
        regressor.n_components = 1

    regressor.fit(X, y)
    funcs = ["decision_function", "predict_proba", "predict_log_proba"]
    for func_name in funcs:
        func = getattr(regressor, func_name, None)
        if func is None:
            # doesn't have function
            continue
        # has function. Should raise deprecation warning
        msg = func_name
        assert_warns_message(DeprecationWarning, msg, func, X)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_class_weight_classifiers(name, classifier_orig):
    if name == "NuSVC":
        # the sparse version has a parameter that doesn't do anything
        raise SkipTest("Not testing NuSVC class weight as it is ignored.")
    if name.endswith("NB"):
        # NaiveBayes classifiers have a somewhat different interface.
        # FIXME SOON!
        raise SkipTest

    for n_centers in [2, 3]:
        # create a very noisy dataset
        X, y = make_blobs(centers=n_centers, random_state=0, cluster_std=20)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,
                                                            random_state=0)

        # can't use gram_if_pairwise() here, setting up gram matrix manually
        if _is_pairwise(classifier_orig):
            X_test = rbf_kernel(X_test, X_train)
            X_train = rbf_kernel(X_train, X_train)

        n_centers = len(np.unique(y_train))

        if n_centers == 2:
            class_weight = {0: 1000, 1: 0.0001}
        else:
            class_weight = {0: 1000, 1: 0.0001, 2: 0.0001}

        classifier = clone(classifier_orig).set_params(
            class_weight=class_weight)
        if hasattr(classifier, "n_iter"):
            classifier.set_params(n_iter=100)
        if hasattr(classifier, "max_iter"):
            classifier.set_params(max_iter=1000)
        if hasattr(classifier, "min_weight_fraction_leaf"):
            classifier.set_params(min_weight_fraction_leaf=0.01)
        if hasattr(classifier, "n_iter_no_change"):
            classifier.set_params(n_iter_no_change=20)

        set_random_state(classifier)
        classifier.fit(X_train, y_train)
        y_pred = classifier.predict(X_test)
        # XXX: Generally can use 0.89 here. On Windows, LinearSVC gets
        #      0.88 (Issue #9111)
        assert_greater(np.mean(y_pred == 0), 0.87)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_class_weight_balanced_classifiers(name, classifier_orig, X_train,
                                            y_train, X_test, y_test, weights):
    classifier = clone(classifier_orig)
    if hasattr(classifier, "n_iter"):
        classifier.set_params(n_iter=100)
    if hasattr(classifier, "max_iter"):
        classifier.set_params(max_iter=1000)

    set_random_state(classifier)
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)

    classifier.set_params(class_weight='balanced')
    classifier.fit(X_train, y_train)
    y_pred_balanced = classifier.predict(X_test)
    assert_greater(f1_score(y_test, y_pred_balanced, average='weighted'),
                   f1_score(y_test, y_pred, average='weighted'))


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_class_weight_balanced_linear_classifier(name, Classifier):
    """Test class weights with non-contiguous class labels."""
    # this is run on classes, not instances, though this should be changed
    X = np.array([[-1.0, -1.0], [-1.0, 0], [-.8, -1.0],
                  [1.0, 1.0], [1.0, 0.0]])
    y = np.array([1, 1, 1, -1, -1])

    classifier = Classifier()

    if hasattr(classifier, "n_iter"):
        # This is a very small dataset, default n_iter are likely to prevent
        # convergence
        classifier.set_params(n_iter=1000)
    if hasattr(classifier, "max_iter"):
        classifier.set_params(max_iter=1000)
    set_random_state(classifier)

    # Let the model compute the class frequencies
    classifier.set_params(class_weight='balanced')
    coef_balanced = classifier.fit(X, y).coef_.copy()

    # Count each label occurrence to reweight manually
    n_samples = len(y)
    n_classes = float(len(np.unique(y)))

    class_weight = {1: n_samples / (np.sum(y == 1) * n_classes),
                    -1: n_samples / (np.sum(y == -1) * n_classes)}
    classifier.set_params(class_weight=class_weight)
    coef_manual = classifier.fit(X, y).coef_.copy()

    assert_allclose(coef_balanced, coef_manual,
                    err_msg="Classifier %s is not computing"
                    " class_weight=balanced properly."
                    % name)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_estimators_overwrite_params(name, estimator_orig):
    X, y = make_blobs(random_state=0, n_samples=9)
    # some want non-negative input
    X -= X.min()
    X = pairwise_estimator_convert_X(X, estimator_orig, kernel=rbf_kernel)
    estimator = clone(estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    set_random_state(estimator)

    # Make a physical copy of the original estimator parameters before fitting.
    params = estimator.get_params()
    original_params = deepcopy(params)

    # Fit the model
    estimator.fit(X, y)

    # Compare the state of the model parameters with the original parameters
    new_params = estimator.get_params()
    for param_name, original_value in original_params.items():
        new_value = new_params[param_name]

        # We should never change or mutate the internal state of input
        # parameters by default. To check this we use the joblib.hash function
        # that introspects recursively any subobjects to compute a checksum.
        # The only exception to this rule of immutable constructor parameters
        # is possible RandomState instance but in this check we explicitly
        # fixed the random_state params recursively to be integer seeds.
        assert_equal(_joblib.hash(new_value), _joblib.hash(original_value),
                     "Estimator %s should not change or mutate "
                     " the parameter %s from %s to %s during fit."
                     % (name, param_name, original_value, new_value))


def check_no_attributes_set_in_init(name, estimator):
    """Check setting during init. """

    if hasattr(type(estimator).__init__, "deprecated_original"):
        return

    init_params = _get_args(type(estimator).__init__)
    if IS_PYPY:
        # __init__ signature has additional objects in PyPy
        for key in ['obj']:
            if key in init_params:
                init_params.remove(key)
    parents_init_params = [param for params_parent in
                           (_get_args(parent) for parent in
                            type(estimator).__mro__)
                           for param in params_parent]

    # Test for no setting apart from parameters during init
    invalid_attr = (set(vars(estimator)) - set(init_params)
                    - set(parents_init_params))
    assert not invalid_attr, (
            "Estimator %s should not set any attribute apart"
            " from parameters during init. Found attributes %s."
            % (name, sorted(invalid_attr)))
    # Ensure that each parameter is set in init
    invalid_attr = set(init_params) - set(vars(estimator)) - {"self"}
    assert not invalid_attr, (
            "Estimator %s should store all parameters"
            " as an attribute during init. Did not find "
            "attributes %s."
            % (name, sorted(invalid_attr)))


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_sparsify_coefficients(name, estimator_orig):
    X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1],
                  [-1, -2], [2, 2], [-2, -2]])
    y = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    est = clone(estimator_orig)

    est.fit(X, y)
    pred_orig = est.predict(X)

    # test sparsify with dense inputs
    est.sparsify()
    assert sparse.issparse(est.coef_)
    pred = est.predict(X)
    assert_array_equal(pred, pred_orig)

    # pickle and unpickle with sparse coef_
    est = pickle.loads(pickle.dumps(est))
    assert sparse.issparse(est.coef_)
    pred = est.predict(X)
    assert_array_equal(pred, pred_orig)


@ignore_warnings(category=DeprecationWarning)
def check_classifier_data_not_an_array(name, estimator_orig):
    X = np.array([[3, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 1]])
    X = pairwise_estimator_convert_X(X, estimator_orig)
    y = [1, 1, 1, 2, 2, 2]
    y = multioutput_estimator_convert_y_2d(estimator_orig, y)
    check_estimators_data_not_an_array(name, estimator_orig, X, y)


@ignore_warnings(category=DeprecationWarning)
def check_regressor_data_not_an_array(name, estimator_orig):
    X, y = _boston_subset(n_samples=50)
    X = pairwise_estimator_convert_X(X, estimator_orig)
    y = multioutput_estimator_convert_y_2d(estimator_orig, y)
    check_estimators_data_not_an_array(name, estimator_orig, X, y)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_estimators_data_not_an_array(name, estimator_orig, X, y):
    if name in CROSS_DECOMPOSITION:
        raise SkipTest("Skipping check_estimators_data_not_an_array "
                       "for cross decomposition module as estimators "
                       "are not deterministic.")
    # separate estimators to control random seeds
    estimator_1 = clone(estimator_orig)
    estimator_2 = clone(estimator_orig)
    set_random_state(estimator_1)
    set_random_state(estimator_2)

    y_ = NotAnArray(np.asarray(y))
    X_ = NotAnArray(np.asarray(X))

    # fit
    estimator_1.fit(X_, y_)
    pred1 = estimator_1.predict(X_)
    estimator_2.fit(X, y)
    pred2 = estimator_2.predict(X)
    assert_allclose(pred1, pred2, atol=1e-2, err_msg=name)


def check_parameters_default_constructible(name, Estimator):
    # this check works on classes, not instances
    # test default-constructibility
    # get rid of deprecation warnings
    with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
        required_parameters = getattr(Estimator, "_required_parameters", [])
        if required_parameters:
            if required_parameters in (["base_estimator"], ["estimator"]):
                if issubclass(Estimator, RegressorMixin):
                    estimator = Estimator(Ridge())
                else:
                    estimator = Estimator(LinearDiscriminantAnalysis())
            else:
                raise SkipTest("Can't instantiate estimator {} which"
                               " requires parameters {}".format(
                                   name, required_parameters))
        else:
            estimator = Estimator()
        # test cloning
        clone(estimator)
        # test __repr__
        repr(estimator)
        # test that set_params returns self
        assert estimator.set_params() is estimator

        # test if init does nothing but set parameters
        # this is important for grid_search etc.
        # We get the default parameters from init and then
        # compare these against the actual values of the attributes.

        # this comes from getattr. Gets rid of deprecation decorator.
        init = getattr(estimator.__init__, 'deprecated_original',
                       estimator.__init__)

        try:
            def param_filter(p):
                """Identify hyper parameters of an estimator"""
                return (p.name != 'self' and
                        p.kind != p.VAR_KEYWORD and
                        p.kind != p.VAR_POSITIONAL)

            init_params = [p for p in signature(init).parameters.values()
                           if param_filter(p)]

        except (TypeError, ValueError):
            # init is not a python function.
            # true for mixins
            return
        params = estimator.get_params()
        if required_parameters == ["estimator"]:
            # they can need a non-default argument
            init_params = init_params[1:]

        for init_param in init_params:
            assert_not_equal(init_param.default, init_param.empty,
                             "parameter %s for %s has no default value"
                             % (init_param.name, type(estimator).__name__))
            if type(init_param.default) is type:
                assert_in(init_param.default, [np.float64, np.int64])
            else:
                assert_in(type(init_param.default),
                          [str, int, float, bool, tuple, type(None),
                           np.float64, types.FunctionType, _joblib.Memory])
            if init_param.name not in params.keys():
                # deprecated parameter, not in get_params
                assert init_param.default is None
                continue

            param_value = params[init_param.name]
            if isinstance(param_value, np.ndarray):
                assert_array_equal(param_value, init_param.default)
            else:
                if is_scalar_nan(param_value):
                    # Allows to set default parameters to np.nan
                    assert param_value is init_param.default, init_param.name
                else:
                    assert param_value == init_param.default, init_param.name


def multioutput_estimator_convert_y_2d(estimator, y):
    # Estimators in mono_output_task_error raise ValueError if y is of 1-D
    # Convert into a 2-D y for those estimators.
    if _safe_tags(estimator, "multioutput_only"):
        return np.reshape(y, (-1, 1))
    return y


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_non_transformer_estimators_n_iter(name, estimator_orig):
    # Test that estimators that are not transformers with a parameter
    # max_iter, return the attribute of n_iter_ at least 1.

    # These models are dependent on external solvers like
    # libsvm and accessing the iter parameter is non-trivial.
    not_run_check_n_iter = ['Ridge', 'SVR', 'NuSVR', 'NuSVC',
                            'RidgeClassifier', 'SVC', 'RandomizedLasso',
                            'LogisticRegressionCV', 'LinearSVC',
                            'LogisticRegression']

    # Tested in test_transformer_n_iter
    not_run_check_n_iter += CROSS_DECOMPOSITION
    if name in not_run_check_n_iter:
        return

    # LassoLars stops early for the default alpha=1.0 the iris dataset.
    if name == 'LassoLars':
        estimator = clone(estimator_orig).set_params(alpha=0.)
    else:
        estimator = clone(estimator_orig)
    if hasattr(estimator, 'max_iter'):
        iris = load_iris()
        X, y_ = iris.data, iris.target
        y_ = multioutput_estimator_convert_y_2d(estimator, y_)

        set_random_state(estimator, 0)
        if name == 'AffinityPropagation':
            estimator.fit(X)
        else:
            estimator.fit(X, y_)

        assert estimator.n_iter_ >= 1


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_transformer_n_iter(name, estimator_orig):
    # Test that transformers with a parameter max_iter, return the
    # attribute of n_iter_ at least 1.
    estimator = clone(estimator_orig)
    if hasattr(estimator, "max_iter"):
        if name in CROSS_DECOMPOSITION:
            # Check using default data
            X = [[0., 0., 1.], [1., 0., 0.], [2., 2., 2.], [2., 5., 4.]]
            y_ = [[0.1, -0.2], [0.9, 1.1], [0.1, -0.5], [0.3, -0.2]]

        else:
            X, y_ = make_blobs(n_samples=30, centers=[[0, 0, 0], [1, 1, 1]],
                               random_state=0, n_features=2, cluster_std=0.1)
            X -= X.min() - 0.1
        set_random_state(estimator, 0)
        estimator.fit(X, y_)

        # These return a n_iter per component.
        if name in CROSS_DECOMPOSITION:
            for iter_ in estimator.n_iter_:
                assert_greater_equal(iter_, 1)
        else:
            assert_greater_equal(estimator.n_iter_, 1)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_get_params_invariance(name, estimator_orig):
    # Checks if get_params(deep=False) is a subset of get_params(deep=True)
    e = clone(estimator_orig)

    shallow_params = e.get_params(deep=False)
    deep_params = e.get_params(deep=True)

    assert all(item in deep_params.items() for item in
               shallow_params.items())


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_set_params(name, estimator_orig):
    # Check that get_params() returns the same thing
    # before and after set_params() with some fuzz
    estimator = clone(estimator_orig)

    orig_params = estimator.get_params(deep=False)
    msg = ("get_params result does not match what was passed to set_params")

    estimator.set_params(**orig_params)
    curr_params = estimator.get_params(deep=False)
    assert_equal(set(orig_params.keys()), set(curr_params.keys()), msg)
    for k, v in curr_params.items():
        assert orig_params[k] is v, msg

    # some fuzz values
    test_values = [-np.inf, np.inf, None]

    test_params = deepcopy(orig_params)
    for param_name in orig_params.keys():
        default_value = orig_params[param_name]
        for value in test_values:
            test_params[param_name] = value
            try:
                estimator.set_params(**test_params)
            except (TypeError, ValueError) as e:
                e_type = e.__class__.__name__
                # Exception occurred, possibly parameter validation
                warnings.warn("{0} occurred during set_params of param {1} on "
                              "{2}. It is recommended to delay parameter "
                              "validation until fit.".format(e_type,
                                                             param_name,
                                                             name))

                change_warning_msg = "Estimator's parameters changed after " \
                                     "set_params raised {}".format(e_type)
                params_before_exception = curr_params
                curr_params = estimator.get_params(deep=False)
                try:
                    assert_equal(set(params_before_exception.keys()),
                                 set(curr_params.keys()))
                    for k, v in curr_params.items():
                        assert params_before_exception[k] is v
                except AssertionError:
                    warnings.warn(change_warning_msg)
            else:
                curr_params = estimator.get_params(deep=False)
                assert_equal(set(test_params.keys()),
                             set(curr_params.keys()),
                             msg)
                for k, v in curr_params.items():
                    assert test_params[k] is v, msg
        test_params[param_name] = default_value


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_classifiers_regression_target(name, estimator_orig):
    # Check if classifier throws an exception when fed regression targets

    boston = load_boston()
    X, y = boston.data, boston.target
    e = clone(estimator_orig)
    msg = 'Unknown label type: '
    if not _safe_tags(e, "no_validation"):
        assert_raises_regex(ValueError, msg, e.fit, X, y)


@ignore_warnings(category=(DeprecationWarning, FutureWarning))
def check_decision_proba_consistency(name, estimator_orig):
    # Check whether an estimator having both decision_function and
    # predict_proba methods has outputs with perfect rank correlation.

    centers = [(2, 2), (4, 4)]
    X, y = make_blobs(n_samples=100, random_state=0, n_features=4,
                      centers=centers, cluster_std=1.0, shuffle=True)
    X_test = np.random.randn(20, 2) + 4
    estimator = clone(estimator_orig)

    if (hasattr(estimator, "decision_function") and
            hasattr(estimator, "predict_proba")):

        estimator.fit(X, y)
        # Since the link function from decision_function() to predict_proba()
        # is sometimes not precise enough (typically expit), we round to the
        # 10th decimal to avoid numerical issues.
        a = estimator.predict_proba(X_test)[:, 1].round(decimals=10)
        b = estimator.decision_function(X_test).round(decimals=10)
        assert_array_equal(rankdata(a), rankdata(b))


def check_outliers_fit_predict(name, estimator_orig):
    # Check fit_predict for outlier detectors.

    n_samples = 300
    X, _ = make_blobs(n_samples=n_samples, random_state=0)
    X = shuffle(X, random_state=7)
    n_samples, n_features = X.shape
    estimator = clone(estimator_orig)

    set_random_state(estimator)

    y_pred = estimator.fit_predict(X)
    assert y_pred.shape == (n_samples,)
    assert y_pred.dtype.kind == 'i'
    assert_array_equal(np.unique(y_pred), np.array([-1, 1]))

    # check fit_predict = fit.predict when the estimator has both a predict and
    # a fit_predict method. recall that it is already assumed here that the
    # estimator has a fit_predict method
    if hasattr(estimator, 'predict'):
        y_pred_2 = estimator.fit(X).predict(X)
        assert_array_equal(y_pred, y_pred_2)

    if hasattr(estimator, "contamination"):
        # proportion of outliers equal to contamination parameter when not
        # set to 'auto'
        expected_outliers = 30
        contamination = float(expected_outliers)/n_samples
        estimator.set_params(contamination=contamination)
        y_pred = estimator.fit_predict(X)

        num_outliers = np.sum(y_pred != 1)
        # num_outliers should be equal to expected_outliers unless
        # there are ties in the decision_function values. this can
        # only be tested for estimators with a decision_function
        # method
        if (num_outliers != expected_outliers and
                hasattr(estimator, 'decision_function')):
            decision = estimator.decision_function(X)
            check_outlier_corruption(num_outliers, expected_outliers, decision)

        # raises error when contamination is a scalar and not in [0,1]
        for contamination in [-0.5, 2.3]:
            estimator.set_params(contamination=contamination)
            assert_raises(ValueError, estimator.fit_predict, X)


def check_fit_idempotent(name, estimator_orig):
    # Check that est.fit(X) is the same as est.fit(X).fit(X). Ideally we would
    # check that the estimated parameters during training (e.g. coefs_) are
    # the same, but having a universal comparison function for those
    # attributes is difficult and full of edge cases. So instead we check that
    # predict(), predict_proba(), decision_function() and transform() return
    # the same results.

    check_methods = ["predict", "transform", "decision_function",
                     "predict_proba"]
    rng = np.random.RandomState(0)

    estimator = clone(estimator_orig)
    set_random_state(estimator)
    if 'warm_start' in estimator.get_params().keys():
        estimator.set_params(warm_start=False)

    n_samples = 100
    X = rng.normal(loc=100, size=(n_samples, 2))
    X = pairwise_estimator_convert_X(X, estimator)
    if is_regressor(estimator_orig):
        y = rng.normal(size=n_samples)
    else:
        y = rng.randint(low=0, high=2, size=n_samples)
    y = multioutput_estimator_convert_y_2d(estimator, y)

    train, test = next(ShuffleSplit(test_size=.2, random_state=rng).split(X))
    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, y_test = _safe_split(estimator, X, y, test, train)

    # Fit for the first time
    estimator.fit(X_train, y_train)

    result = {method: getattr(estimator, method)(X_test)
              for method in check_methods
              if hasattr(estimator, method)}

    # Fit again
    set_random_state(estimator)
    estimator.fit(X_train, y_train)

    for method in check_methods:
        if hasattr(estimator, method):
            new_result = getattr(estimator, method)(X_test)
            assert_allclose_dense_sparse(result[method], new_result)
