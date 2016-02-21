"""Test the validation module"""
from __future__ import division

import sys
import warnings

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_warns
from sklearn.utils.mocking import CheckingClassifier, MockDataFrame

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import permutation_test_score
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import LeaveOneLabelOut
from sklearn.model_selection import LeavePLabelOut
from sklearn.model_selection import LabelKFold
from sklearn.model_selection import LabelShuffleSplit
from sklearn.model_selection import learning_curve
from sklearn.model_selection import validation_curve
from sklearn.model_selection._validation import _check_is_permutation

from sklearn.datasets import make_regression
from sklearn.datasets import load_boston
from sklearn.datasets import load_iris
from sklearn.metrics import explained_variance_score
from sklearn.metrics import make_scorer
from sklearn.metrics import precision_score

from sklearn.linear_model import Ridge
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.cluster import KMeans

from sklearn.preprocessing import Imputer
from sklearn.pipeline import Pipeline

from sklearn.externals.six.moves import cStringIO as StringIO
from sklearn.base import BaseEstimator
from sklearn.multiclass import OneVsRestClassifier
from sklearn.datasets import make_classification
from sklearn.datasets import make_multilabel_classification

from test_split import MockClassifier


class MockImprovingEstimator(BaseEstimator):
    """Dummy classifier to test the learning curve"""
    def __init__(self, n_max_train_sizes):
        self.n_max_train_sizes = n_max_train_sizes
        self.train_sizes = 0
        self.X_subset = None

    def fit(self, X_subset, y_subset=None):
        self.X_subset = X_subset
        self.train_sizes = X_subset.shape[0]
        return self

    def predict(self, X):
        raise NotImplementedError

    def score(self, X=None, Y=None):
        # training score becomes worse (2 -> 1), test error better (0 -> 1)
        if self._is_training_data(X):
            return 2. - float(self.train_sizes) / self.n_max_train_sizes
        else:
            return float(self.train_sizes) / self.n_max_train_sizes

    def _is_training_data(self, X):
        return X is self.X_subset


class MockIncrementalImprovingEstimator(MockImprovingEstimator):
    """Dummy classifier that provides partial_fit"""
    def __init__(self, n_max_train_sizes):
        super(MockIncrementalImprovingEstimator,
              self).__init__(n_max_train_sizes)
        self.x = None

    def _is_training_data(self, X):
        return self.x in X

    def partial_fit(self, X, y=None, **params):
        self.train_sizes += X.shape[0]
        self.x = X[0]


class MockEstimatorWithParameter(BaseEstimator):
    """Dummy classifier to test the validation curve"""
    def __init__(self, param=0.5):
        self.X_subset = None
        self.param = param

    def fit(self, X_subset, y_subset):
        self.X_subset = X_subset
        self.train_sizes = X_subset.shape[0]
        return self

    def predict(self, X):
        raise NotImplementedError

    def score(self, X=None, y=None):
        return self.param if self._is_training_data(X) else 1 - self.param

    def _is_training_data(self, X):
        return X is self.X_subset


# XXX: use 2D array, since 1D X is being detected as a single sample in
# check_consistent_length
X = np.ones((10, 2))
X_sparse = coo_matrix(X)
y = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
# The number of samples per class needs to be > n_folds, for StratifiedKFold(3)
y2 = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3, 3])


def test_cross_val_score():
    clf = MockClassifier()

    for a in range(-10, 10):
        clf.a = a
        # Smoke test
        scores = cross_val_score(clf, X, y2)
        assert_array_equal(scores, clf.score(X, y2))

        # test with multioutput y
        multioutput_y = np.column_stack([y2, y2[::-1]])
        scores = cross_val_score(clf, X_sparse, multioutput_y)
        assert_array_equal(scores, clf.score(X_sparse, multioutput_y))

        scores = cross_val_score(clf, X_sparse, y2)
        assert_array_equal(scores, clf.score(X_sparse, y2))

        # test with multioutput y
        scores = cross_val_score(clf, X_sparse, multioutput_y)
        assert_array_equal(scores, clf.score(X_sparse, multioutput_y))

    # test with X and y as list
    list_check = lambda x: isinstance(x, list)
    clf = CheckingClassifier(check_X=list_check)
    scores = cross_val_score(clf, X.tolist(), y2.tolist())

    clf = CheckingClassifier(check_y=list_check)
    scores = cross_val_score(clf, X, y2.tolist())

    assert_raises(ValueError, cross_val_score, clf, X, y2, scoring="sklearn")

    # test with 3d X and
    X_3d = X[:, :, np.newaxis]
    clf = MockClassifier(allow_nd=True)
    scores = cross_val_score(clf, X_3d, y2)

    clf = MockClassifier(allow_nd=False)
    assert_raises(ValueError, cross_val_score, clf, X_3d, y2)


def test_cross_val_score_predict_labels():
    # Check if ValueError (when labels is None) propagates to cross_val_score
    # and cross_val_predict
    # And also check if labels is correctly passed to the cv object
    X, y = make_classification(n_samples=20, n_classes=2, random_state=0)

    clf = SVC(kernel="linear")

    label_cvs = [LeaveOneLabelOut(), LeavePLabelOut(2), LabelKFold(),
                 LabelShuffleSplit()]
    for cv in label_cvs:
        assert_raise_message(ValueError,
                             "The labels parameter should not be None",
                             cross_val_score, estimator=clf, X=X, y=y, cv=cv)
        assert_raise_message(ValueError,
                             "The labels parameter should not be None",
                             cross_val_predict, estimator=clf, X=X, y=y, cv=cv)


def test_cross_val_score_pandas():
    # check cross_val_score doesn't destroy pandas dataframe
    types = [(MockDataFrame, MockDataFrame)]
    try:
        from pandas import Series, DataFrame
        types.append((Series, DataFrame))
    except ImportError:
        pass
    for TargetType, InputFeatureType in types:
        # X dataframe, y series
        # 3 fold cross val is used so we need atleast 3 samples per class
        X_df, y_ser = InputFeatureType(X), TargetType(y2)
        check_df = lambda x: isinstance(x, InputFeatureType)
        check_series = lambda x: isinstance(x, TargetType)
        clf = CheckingClassifier(check_X=check_df, check_y=check_series)
        cross_val_score(clf, X_df, y_ser)


def test_cross_val_score_mask():
    # test that cross_val_score works with boolean masks
    svm = SVC(kernel="linear")
    iris = load_iris()
    X, y = iris.data, iris.target
    kfold = KFold(5)
    scores_indices = cross_val_score(svm, X, y, cv=kfold)
    kfold = KFold(5)
    cv_masks = []
    for train, test in kfold.split(X, y):
        mask_train = np.zeros(len(y), dtype=np.bool)
        mask_test = np.zeros(len(y), dtype=np.bool)
        mask_train[train] = 1
        mask_test[test] = 1
        cv_masks.append((train, test))
    scores_masks = cross_val_score(svm, X, y, cv=cv_masks)
    assert_array_equal(scores_indices, scores_masks)


def test_cross_val_score_precomputed():
    # test for svm with precomputed kernel
    svm = SVC(kernel="precomputed")
    iris = load_iris()
    X, y = iris.data, iris.target
    linear_kernel = np.dot(X, X.T)
    score_precomputed = cross_val_score(svm, linear_kernel, y)
    svm = SVC(kernel="linear")
    score_linear = cross_val_score(svm, X, y)
    assert_array_equal(score_precomputed, score_linear)

    # Error raised for non-square X
    svm = SVC(kernel="precomputed")
    assert_raises(ValueError, cross_val_score, svm, X, y)

    # test error is raised when the precomputed kernel is not array-like
    # or sparse
    assert_raises(ValueError, cross_val_score, svm,
                  linear_kernel.tolist(), y)


def test_cross_val_score_fit_params():
    clf = MockClassifier()
    n_samples = X.shape[0]
    n_classes = len(np.unique(y))

    W_sparse = coo_matrix((np.array([1]), (np.array([1]), np.array([0]))),
                          shape=(10, 1))
    P_sparse = coo_matrix(np.eye(5))

    DUMMY_INT = 42
    DUMMY_STR = '42'
    DUMMY_OBJ = object()

    def assert_fit_params(clf):
        # Function to test that the values are passed correctly to the
        # classifier arguments for non-array type

        assert_equal(clf.dummy_int, DUMMY_INT)
        assert_equal(clf.dummy_str, DUMMY_STR)
        assert_equal(clf.dummy_obj, DUMMY_OBJ)

    fit_params = {'sample_weight': np.ones(n_samples),
                  'class_prior': np.ones(n_classes) / n_classes,
                  'sparse_sample_weight': W_sparse,
                  'sparse_param': P_sparse,
                  'dummy_int': DUMMY_INT,
                  'dummy_str': DUMMY_STR,
                  'dummy_obj': DUMMY_OBJ,
                  'callback': assert_fit_params}
    cross_val_score(clf, X, y, fit_params=fit_params)


def test_cross_val_score_score_func():
    clf = MockClassifier()
    _score_func_args = []

    def score_func(y_test, y_predict):
        _score_func_args.append((y_test, y_predict))
        return 1.0

    with warnings.catch_warnings(record=True):
        scoring = make_scorer(score_func)
        score = cross_val_score(clf, X, y, scoring=scoring)
    assert_array_equal(score, [1.0, 1.0, 1.0])
    assert len(_score_func_args) == 3


def test_cross_val_score_errors():
    class BrokenEstimator:
        pass

    assert_raises(TypeError, cross_val_score, BrokenEstimator(), X)


def test_cross_val_score_with_score_func_classification():
    iris = load_iris()
    clf = SVC(kernel='linear')

    # Default score (should be the accuracy score)
    scores = cross_val_score(clf, iris.data, iris.target, cv=5)
    assert_array_almost_equal(scores, [0.97, 1., 0.97, 0.97, 1.], 2)

    # Correct classification score (aka. zero / one score) - should be the
    # same as the default estimator score
    zo_scores = cross_val_score(clf, iris.data, iris.target,
                                scoring="accuracy", cv=5)
    assert_array_almost_equal(zo_scores, [0.97, 1., 0.97, 0.97, 1.], 2)

    # F1 score (class are balanced so f1_score should be equal to zero/one
    # score
    f1_scores = cross_val_score(clf, iris.data, iris.target,
                                scoring="f1_weighted", cv=5)
    assert_array_almost_equal(f1_scores, [0.97, 1., 0.97, 0.97, 1.], 2)


def test_cross_val_score_with_score_func_regression():
    X, y = make_regression(n_samples=30, n_features=20, n_informative=5,
                           random_state=0)
    reg = Ridge()

    # Default score of the Ridge regression estimator
    scores = cross_val_score(reg, X, y, cv=5)
    assert_array_almost_equal(scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # R2 score (aka. determination coefficient) - should be the
    # same as the default estimator score
    r2_scores = cross_val_score(reg, X, y, scoring="r2", cv=5)
    assert_array_almost_equal(r2_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # Mean squared error; this is a loss function, so "scores" are negative
    mse_scores = cross_val_score(reg, X, y, cv=5, scoring="mean_squared_error")
    expected_mse = np.array([-763.07, -553.16, -274.38, -273.26, -1681.99])
    assert_array_almost_equal(mse_scores, expected_mse, 2)

    # Explained variance
    scoring = make_scorer(explained_variance_score)
    ev_scores = cross_val_score(reg, X, y, cv=5, scoring=scoring)
    assert_array_almost_equal(ev_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)


def test_permutation_score():
    iris = load_iris()
    X = iris.data
    X_sparse = coo_matrix(X)
    y = iris.target
    svm = SVC(kernel='linear')
    cv = StratifiedKFold(2)

    score, scores, pvalue = permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy")
    assert_greater(score, 0.9)
    assert_almost_equal(pvalue, 0.0, 1)

    score_label, _, pvalue_label = permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy",
        labels=np.ones(y.size), random_state=0)
    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # check that we obtain the same results with a sparse representation
    svm_sparse = SVC(kernel='linear')
    cv_sparse = StratifiedKFold(2)
    score_label, _, pvalue_label = permutation_test_score(
        svm_sparse, X_sparse, y, n_permutations=30, cv=cv_sparse,
        scoring="accuracy", labels=np.ones(y.size), random_state=0)

    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # test with custom scoring object
    def custom_score(y_true, y_pred):
        return (((y_true == y_pred).sum() - (y_true != y_pred).sum())
                / y_true.shape[0])

    scorer = make_scorer(custom_score)
    score, _, pvalue = permutation_test_score(
        svm, X, y, n_permutations=100, scoring=scorer, cv=cv, random_state=0)
    assert_almost_equal(score, .93, 2)
    assert_almost_equal(pvalue, 0.01, 3)

    # set random y
    y = np.mod(np.arange(len(y)), 3)

    score, scores, pvalue = permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy")

    assert_less(score, 0.5)
    assert_greater(pvalue, 0.2)


def test_permutation_test_score_allow_nans():
    # Check that permutation_test_score allows input data with NaNs
    X = np.arange(200, dtype=np.float64).reshape(10, -1)
    X[2, :] = np.nan
    y = np.repeat([0, 1], X.shape[0] / 2)
    p = Pipeline([
        ('imputer', Imputer(strategy='mean', missing_values='NaN')),
        ('classifier', MockClassifier()),
    ])
    permutation_test_score(p, X, y, cv=5)


def test_cross_val_score_allow_nans():
    # Check that cross_val_score allows input data with NaNs
    X = np.arange(200, dtype=np.float64).reshape(10, -1)
    X[2, :] = np.nan
    y = np.repeat([0, 1], X.shape[0] / 2)
    p = Pipeline([
        ('imputer', Imputer(strategy='mean', missing_values='NaN')),
        ('classifier', MockClassifier()),
    ])
    cross_val_score(p, X, y, cv=5)


def test_cross_val_score_multilabel():
    X = np.array([[-3, 4], [2, 4], [3, 3], [0, 2], [-3, 1],
                  [-2, 1], [0, 0], [-2, -1], [-1, -2], [1, -2]])
    y = np.array([[1, 1], [0, 1], [0, 1], [0, 1], [1, 1],
                  [0, 1], [1, 0], [1, 1], [1, 0], [0, 0]])
    clf = KNeighborsClassifier(n_neighbors=1)
    scoring_micro = make_scorer(precision_score, average='micro')
    scoring_macro = make_scorer(precision_score, average='macro')
    scoring_samples = make_scorer(precision_score, average='samples')
    score_micro = cross_val_score(clf, X, y, scoring=scoring_micro, cv=5)
    score_macro = cross_val_score(clf, X, y, scoring=scoring_macro, cv=5)
    score_samples = cross_val_score(clf, X, y, scoring=scoring_samples, cv=5)
    assert_almost_equal(score_micro, [1, 1 / 2, 3 / 4, 1 / 2, 1 / 3])
    assert_almost_equal(score_macro, [1, 1 / 2, 3 / 4, 1 / 2, 1 / 4])
    assert_almost_equal(score_samples, [1, 1 / 2, 3 / 4, 1 / 2, 1 / 4])


def test_cross_val_predict():
    boston = load_boston()
    X, y = boston.data, boston.target
    cv = KFold()

    est = Ridge()

    # Naive loop (should be same as cross_val_predict):
    preds2 = np.zeros_like(y)
    for train, test in cv.split(X, y):
        est.fit(X[train], y[train])
        preds2[test] = est.predict(X[test])

    preds = cross_val_predict(est, X, y, cv=cv)
    assert_array_almost_equal(preds, preds2)

    preds = cross_val_predict(est, X, y)
    assert_equal(len(preds), len(y))

    cv = LeaveOneOut()
    preds = cross_val_predict(est, X, y, cv=cv)
    assert_equal(len(preds), len(y))

    Xsp = X.copy()
    Xsp *= (Xsp > np.median(Xsp))
    Xsp = coo_matrix(Xsp)
    preds = cross_val_predict(est, Xsp, y)
    assert_array_almost_equal(len(preds), len(y))

    preds = cross_val_predict(KMeans(), X)
    assert_equal(len(preds), len(y))

    class BadCV():
        def split(self, X, y=None, labels=None):
            for i in range(4):
                yield np.array([0, 1, 2, 3]), np.array([4, 5, 6, 7, 8])

    assert_raises(ValueError, cross_val_predict, est, X, y, cv=BadCV())


def test_cross_val_predict_input_types():
    iris = load_iris()
    X, y = iris.data, iris.target
    X_sparse = coo_matrix(X)
    multioutput_y = np.column_stack([y, y[::-1]])

    clf = Ridge(fit_intercept=False, random_state=0)
    # 3 fold cv is used --> atleast 3 samples per class
    # Smoke test
    predictions = cross_val_predict(clf, X, y)
    assert_equal(predictions.shape, (150,))

    # test with multioutput y
    predictions = cross_val_predict(clf, X_sparse, multioutput_y)
    assert_equal(predictions.shape, (150, 2))

    predictions = cross_val_predict(clf, X_sparse, y)
    assert_array_equal(predictions.shape, (150,))

    # test with multioutput y
    predictions = cross_val_predict(clf, X_sparse, multioutput_y)
    assert_array_equal(predictions.shape, (150, 2))

    # test with X and y as list
    list_check = lambda x: isinstance(x, list)
    clf = CheckingClassifier(check_X=list_check)
    predictions = cross_val_predict(clf, X.tolist(), y.tolist())

    clf = CheckingClassifier(check_y=list_check)
    predictions = cross_val_predict(clf, X, y.tolist())

    # test with 3d X and
    X_3d = X[:, :, np.newaxis]
    check_3d = lambda x: x.ndim == 3
    clf = CheckingClassifier(check_X=check_3d)
    predictions = cross_val_predict(clf, X_3d, y)
    assert_array_equal(predictions.shape, (150,))


def test_cross_val_predict_pandas():
    # check cross_val_score doesn't destroy pandas dataframe
    types = [(MockDataFrame, MockDataFrame)]
    try:
        from pandas import Series, DataFrame
        types.append((Series, DataFrame))
    except ImportError:
        pass
    for TargetType, InputFeatureType in types:
        # X dataframe, y series
        X_df, y_ser = InputFeatureType(X), TargetType(y2)
        check_df = lambda x: isinstance(x, InputFeatureType)
        check_series = lambda x: isinstance(x, TargetType)
        clf = CheckingClassifier(check_X=check_df, check_y=check_series)
        cross_val_predict(clf, X_df, y_ser)


def test_cross_val_score_sparse_fit_params():
    iris = load_iris()
    X, y = iris.data, iris.target
    clf = MockClassifier()
    fit_params = {'sparse_sample_weight': coo_matrix(np.eye(X.shape[0]))}
    a = cross_val_score(clf, X, y, fit_params=fit_params)
    assert_array_equal(a, np.ones(3))


def test_learning_curve():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)
    with warnings.catch_warnings(record=True) as w:
        train_sizes, train_scores, test_scores = learning_curve(
            estimator, X, y, cv=3, train_sizes=np.linspace(0.1, 1.0, 10))
    if len(w) > 0:
        raise RuntimeError("Unexpected warning: %r" % w[0].message)
    assert_equal(train_scores.shape, (10, 3))
    assert_equal(test_scores.shape, (10, 3))
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores.mean(axis=1),
                              np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores.mean(axis=1),
                              np.linspace(0.1, 1.0, 10))


def test_learning_curve_unsupervised():
    X, _ = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y=None, cv=3, train_sizes=np.linspace(0.1, 1.0, 10))
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores.mean(axis=1),
                              np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores.mean(axis=1),
                              np.linspace(0.1, 1.0, 10))


def test_learning_curve_verbose():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        train_sizes, train_scores, test_scores = \
            learning_curve(estimator, X, y, cv=3, verbose=1)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    assert("[learning_curve]" in out)


def test_learning_curve_incremental_learning_not_possible():
    X, y = make_classification(n_samples=2, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    # The mockup does not have partial_fit()
    estimator = MockImprovingEstimator(1)
    assert_raises(ValueError, learning_curve, estimator, X, y,
                  exploit_incremental_learning=True)


def test_learning_curve_incremental_learning():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockIncrementalImprovingEstimator(20)
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=3, exploit_incremental_learning=True,
        train_sizes=np.linspace(0.1, 1.0, 10))
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores.mean(axis=1),
                              np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores.mean(axis=1),
                              np.linspace(0.1, 1.0, 10))


def test_learning_curve_incremental_learning_unsupervised():
    X, _ = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockIncrementalImprovingEstimator(20)
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y=None, cv=3, exploit_incremental_learning=True,
        train_sizes=np.linspace(0.1, 1.0, 10))
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores.mean(axis=1),
                              np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores.mean(axis=1),
                              np.linspace(0.1, 1.0, 10))


def test_learning_curve_batch_and_incremental_learning_are_equal():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    train_sizes = np.linspace(0.2, 1.0, 5)
    estimator = PassiveAggressiveClassifier(n_iter=1, shuffle=False)

    train_sizes_inc, train_scores_inc, test_scores_inc = \
        learning_curve(
            estimator, X, y, train_sizes=train_sizes,
            cv=3, exploit_incremental_learning=True)
    train_sizes_batch, train_scores_batch, test_scores_batch = \
        learning_curve(
            estimator, X, y, cv=3, train_sizes=train_sizes,
            exploit_incremental_learning=False)

    assert_array_equal(train_sizes_inc, train_sizes_batch)
    assert_array_almost_equal(train_scores_inc.mean(axis=1),
                              train_scores_batch.mean(axis=1))
    assert_array_almost_equal(test_scores_inc.mean(axis=1),
                              test_scores_batch.mean(axis=1))


def test_learning_curve_n_sample_range_out_of_bounds():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0, 1])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0.0, 1.0])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0.1, 1.1])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[0, 20])
    assert_raises(ValueError, learning_curve, estimator, X, y, cv=3,
                  train_sizes=[1, 21])


def test_learning_curve_remove_duplicate_sample_sizes():
    X, y = make_classification(n_samples=3, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(2)
    train_sizes, _, _ = assert_warns(
        RuntimeWarning, learning_curve, estimator, X, y, cv=3,
        train_sizes=np.linspace(0.33, 1.0, 3))
    assert_array_equal(train_sizes, [1, 2])


def test_learning_curve_with_boolean_indices():
    X, y = make_classification(n_samples=30, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    estimator = MockImprovingEstimator(20)
    cv = KFold(n_folds=3)
    train_sizes, train_scores, test_scores = learning_curve(
        estimator, X, y, cv=cv, train_sizes=np.linspace(0.1, 1.0, 10))
    assert_array_equal(train_sizes, np.linspace(2, 20, 10))
    assert_array_almost_equal(train_scores.mean(axis=1),
                              np.linspace(1.9, 1.0, 10))
    assert_array_almost_equal(test_scores.mean(axis=1),
                              np.linspace(0.1, 1.0, 10))


def test_validation_curve():
    X, y = make_classification(n_samples=2, n_features=1, n_informative=1,
                               n_redundant=0, n_classes=2,
                               n_clusters_per_class=1, random_state=0)
    param_range = np.linspace(0, 1, 10)
    with warnings.catch_warnings(record=True) as w:
        train_scores, test_scores = validation_curve(
            MockEstimatorWithParameter(), X, y, param_name="param",
            param_range=param_range, cv=2
        )
    if len(w) > 0:
        raise RuntimeError("Unexpected warning: %r" % w[0].message)

    assert_array_almost_equal(train_scores.mean(axis=1), param_range)
    assert_array_almost_equal(test_scores.mean(axis=1), 1 - param_range)


def test_check_is_permutation():
    p = np.arange(100)
    assert_true(_check_is_permutation(p, 100))
    assert_false(_check_is_permutation(np.delete(p, 23), 100))

    p[0] = 23
    assert_false(_check_is_permutation(p, 100))


def test_cross_val_predict_sparse_prediction():
    # check that cross_val_predict gives same result for sparse and dense input
    X, y = make_multilabel_classification(n_classes=2, n_labels=1,
                                          allow_unlabeled=False,
                                          return_indicator=True,
                                          random_state=1)
    X_sparse = csr_matrix(X)
    y_sparse = csr_matrix(y)
    classif = OneVsRestClassifier(SVC(kernel='linear'))
    preds = cross_val_predict(classif, X, y, cv=10)
    preds_sparse = cross_val_predict(classif, X_sparse, y_sparse, cv=10)
    preds_sparse = preds_sparse.toarray()
    assert_array_almost_equal(preds_sparse, preds)
