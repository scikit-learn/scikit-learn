"""Test the cross_validation module"""
from __future__ import division
import warnings

import numpy as np
from scipy.sparse import coo_matrix
from scipy import stats

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.mocking import CheckingClassifier, MockDataFrame

from sklearn import cross_validation as cval
from sklearn.base import BaseEstimator
from sklearn.datasets import make_regression
from sklearn.datasets import load_digits
from sklearn.datasets import load_iris
from sklearn.metrics import accuracy_score
from sklearn.metrics import explained_variance_score
from sklearn.metrics import make_scorer
from sklearn.metrics import precision_score

from sklearn.externals import six
from sklearn.externals.six.moves import zip

from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

from sklearn.preprocessing import Imputer, LabelBinarizer
from sklearn.pipeline import Pipeline


class MockClassifier(BaseEstimator):
    """Dummy classifier to test the cross-validation"""

    def __init__(self, a=0, allow_nd=False):
        self.a = a
        self.allow_nd = allow_nd

    def fit(self, X, Y=None, sample_weight=None, class_prior=None):
        if self.allow_nd:
            X = X.reshape(len(X), -1)
        if X.ndim >= 3 and not self.allow_nd:
            raise ValueError('X cannot be d')
        if sample_weight is not None:
            assert_true(sample_weight.shape[0] == X.shape[0],
                        'MockClassifier extra fit_param sample_weight.shape[0]'
                        ' is {0}, should be {1}'.format(sample_weight.shape[0],
                                                        X.shape[0]))
        if class_prior is not None:
            assert_true(class_prior.shape[0] == len(np.unique(y)),
                        'MockClassifier extra fit_param class_prior.shape[0]'
                        ' is {0}, should be {1}'.format(class_prior.shape[0],
                                                        len(np.unique(y))))
        return self

    def predict(self, T):
        if self.allow_nd:
            T = T.reshape(len(T), -1)
        return T.shape[0]

    def score(self, X=None, Y=None):
        return 1. / (1 + np.abs(self.a))


X = np.ones((10, 2))
X_sparse = coo_matrix(X)
y = np.arange(10) // 2

##############################################################################
# Tests


def check_valid_split(train, test, n_samples=None):
    # Use python sets to get more informative assertion failure messages
    train, test = set(train), set(test)

    # Train and test split should not overlap
    assert_equal(train.intersection(test), set())

    if n_samples is not None:
        # Check that the union of train an test split cover all the indices
        assert_equal(train.union(test), set(range(n_samples)))


def check_cv_coverage(cv, expected_n_iter=None, n_samples=None):
    # Check that a all the samples appear at least once in a test fold
    if expected_n_iter is not None:
        assert_equal(len(cv), expected_n_iter)
    else:
        expected_n_iter = len(cv)

    collected_test_samples = set()
    iterations = 0
    for train, test in cv:
        check_valid_split(train, test, n_samples=n_samples)
        iterations += 1
        collected_test_samples.update(test)

    # Check that the accumulated test samples cover the whole dataset
    assert_equal(iterations, expected_n_iter)
    if n_samples is not None:
        assert_equal(collected_test_samples, set(range(n_samples)))


def test_kfold_valueerrors():
    # Check that errors are raised if there is not enough samples
    assert_raises(ValueError, cval.KFold, 3, 4)

    # Check that a warning is raised if the least populated class has too few
    # members.
    y = [3, 3, -1, -1, 2]

    cv = assert_warns_message(Warning, "The least populated class",
                              cval.StratifiedKFold, y, 3)

    # Check that despite the warning the folds are still computed even
    # though all the classes are not necessarily represented at on each
    # side of the split at each split
    check_cv_coverage(cv, expected_n_iter=3, n_samples=len(y))

    # Error when number of folds is <= 1
    assert_raises(ValueError, cval.KFold, 2, 0)
    assert_raises(ValueError, cval.KFold, 2, 1)
    assert_raises(ValueError, cval.StratifiedKFold, y, 0)
    assert_raises(ValueError, cval.StratifiedKFold, y, 1)

    # When n is not integer:
    assert_raises(ValueError, cval.KFold, 2.5, 2)

    # When n_folds is not integer:
    assert_raises(ValueError, cval.KFold, 5, 1.5)
    assert_raises(ValueError, cval.StratifiedKFold, y, 1.5)


def test_kfold_indices():
    # Check all indices are returned in the test folds
    kf = cval.KFold(300, 3)
    check_cv_coverage(kf, expected_n_iter=3, n_samples=300)

    # Check all indices are returned in the test folds even when equal-sized
    # folds are not possible
    kf = cval.KFold(17, 3)
    check_cv_coverage(kf, expected_n_iter=3, n_samples=17)


def test_kfold_no_shuffle():
    # Manually check that KFold preserves the data ordering on toy datasets
    splits = iter(cval.KFold(4, 2))
    train, test = next(splits)
    assert_array_equal(test, [0, 1])
    assert_array_equal(train, [2, 3])

    train, test = next(splits)
    assert_array_equal(test, [2, 3])
    assert_array_equal(train, [0, 1])

    splits = iter(cval.KFold(5, 2))
    train, test = next(splits)
    assert_array_equal(test, [0, 1, 2])
    assert_array_equal(train, [3, 4])

    train, test = next(splits)
    assert_array_equal(test, [3, 4])
    assert_array_equal(train, [0, 1, 2])


def test_stratified_kfold_no_shuffle():
    # Manually check that StratifiedKFold preserves the data ordering as much
    # as possible on toy datasets in order to avoid hiding sample dependencies
    # when possible
    splits = iter(cval.StratifiedKFold([1, 1, 0, 0], 2))
    train, test = next(splits)
    assert_array_equal(test, [0, 2])
    assert_array_equal(train, [1, 3])

    train, test = next(splits)
    assert_array_equal(test, [1, 3])
    assert_array_equal(train, [0, 2])

    splits = iter(cval.StratifiedKFold([1, 1, 1, 0, 0, 0, 0], 2))
    train, test = next(splits)
    assert_array_equal(test, [0, 1, 3, 4])
    assert_array_equal(train, [2, 5, 6])

    train, test = next(splits)
    assert_array_equal(test, [2, 5, 6])
    assert_array_equal(train, [0, 1, 3, 4])


def test_stratified_kfold_ratios():
    # Check that stratified kfold preserves label ratios in individual splits
    # Repeat with shuffling turned off and on
    n_samples = 1000
    labels = np.array([4] * int(0.10 * n_samples) +
                      [0] * int(0.89 * n_samples) +
                      [1] * int(0.01 * n_samples))
    for shuffle in [False, True]:
        for train, test in cval.StratifiedKFold(labels, 5, shuffle=shuffle):
            assert_almost_equal(np.sum(labels[train] == 4) / len(train), 0.10,
                                2)
            assert_almost_equal(np.sum(labels[train] == 0) / len(train), 0.89,
                                2)
            assert_almost_equal(np.sum(labels[train] == 1) / len(train), 0.01,
                                2)
            assert_almost_equal(np.sum(labels[test] == 4) / len(test), 0.10, 2)
            assert_almost_equal(np.sum(labels[test] == 0) / len(test), 0.89, 2)
            assert_almost_equal(np.sum(labels[test] == 1) / len(test), 0.01, 2)


def test_kfold_balance():
    # Check that KFold returns folds with balanced sizes
    for kf in [cval.KFold(i, 5) for i in range(11, 17)]:
        sizes = []
        for _, test in kf:
            sizes.append(len(test))

        assert_true((np.max(sizes) - np.min(sizes)) <= 1)
        assert_equal(np.sum(sizes), kf.n)


def test_stratifiedkfold_balance():
    # Check that KFold returns folds with balanced sizes (only when
    # stratification is possible)
    # Repeat with shuffling turned off and on
    labels = [0] * 3 + [1] * 14
    for shuffle in [False, True]:
        for skf in [cval.StratifiedKFold(labels[:i], 3, shuffle=shuffle)
                    for i in range(11, 17)]:
            sizes = []
            for _, test in skf:
                sizes.append(len(test))

            assert_true((np.max(sizes) - np.min(sizes)) <= 1)
            assert_equal(np.sum(sizes), skf.n)


def test_shuffle_kfold():
    # Check the indices are shuffled properly, and that all indices are
    # returned in the different test folds
    kf = cval.KFold(300, 3, shuffle=True, random_state=0)
    ind = np.arange(300)

    all_folds = None
    for train, test in kf:
        sorted_array = np.arange(100)
        assert_true(np.any(sorted_array != ind[train]))
        sorted_array = np.arange(101, 200)
        assert_true(np.any(sorted_array != ind[train]))
        sorted_array = np.arange(201, 300)
        assert_true(np.any(sorted_array != ind[train]))
        if all_folds is None:
            all_folds = ind[test].copy()
        else:
            all_folds = np.concatenate((all_folds, ind[test]))

    all_folds.sort()
    assert_array_equal(all_folds, ind)


def test_shuffle_stratifiedkfold():
    # Check that shuffling is happening when requested, and for proper
    # sample coverage
    labels = [0] * 20 + [1] * 20
    kf0 = list(cval.StratifiedKFold(labels, 5, shuffle=True, random_state=0))
    kf1 = list(cval.StratifiedKFold(labels, 5, shuffle=True, random_state=1))
    for (_, test0), (_, test1) in zip(kf0, kf1):
        assert_true(set(test0) != set(test1))
    check_cv_coverage(kf0, expected_n_iter=5, n_samples=40)


def test_kfold_can_detect_dependent_samples_on_digits():  # see #2372
    # The digits samples are dependent: they are apparently grouped by authors
    # although we don't have any information on the groups segment locations
    # for this data. We can highlight this fact be computing k-fold cross-
    # validation with and without shuffling: we observe that the shuffling case
    # wrongly makes the IID assumption and is therefore too optimistic: it
    # estimates a much higher accuracy (around 0.96) than than the non
    # shuffling variant (around 0.86).

    digits = load_digits()
    X, y = digits.data[:800], digits.target[:800]
    model = SVC(C=10, gamma=0.005)
    n = len(y)

    cv = cval.KFold(n, 5, shuffle=False)
    mean_score = cval.cross_val_score(model, X, y, cv=cv).mean()
    assert_greater(0.88, mean_score)
    assert_greater(mean_score, 0.85)

    # Shuffling the data artificially breaks the dependency and hides the
    # overfitting of the model with regards to the writing style of the authors
    # by yielding a seriously overestimated score:

    cv = cval.KFold(n, 5, shuffle=True, random_state=0)
    mean_score = cval.cross_val_score(model, X, y, cv=cv).mean()
    assert_greater(mean_score, 0.95)

    cv = cval.KFold(n, 5, shuffle=True, random_state=1)
    mean_score = cval.cross_val_score(model, X, y, cv=cv).mean()
    assert_greater(mean_score, 0.95)

    # Similarly, StratifiedKFold should try to shuffle the data as little
    # as possible (while respecting the balanced class constraints)
    # and thus be able to detect the dependency by not overestimating
    # the CV score either. As the digits dataset is approximately balanced
    # the estimated mean score is close to the score measured with
    # non-shuffled KFold

    cv = cval.StratifiedKFold(y, 5)
    mean_score = cval.cross_val_score(model, X, y, cv=cv).mean()
    assert_greater(0.88, mean_score)
    assert_greater(mean_score, 0.85)


def test_shuffle_split():
    ss1 = cval.ShuffleSplit(10, test_size=0.2, random_state=0)
    ss2 = cval.ShuffleSplit(10, test_size=2, random_state=0)
    ss3 = cval.ShuffleSplit(10, test_size=np.int32(2), random_state=0)
    for typ in six.integer_types:
        ss4 = cval.ShuffleSplit(10, test_size=typ(2), random_state=0)
    for t1, t2, t3, t4 in zip(ss1, ss2, ss3, ss4):
        assert_array_equal(t1[0], t2[0])
        assert_array_equal(t2[0], t3[0])
        assert_array_equal(t3[0], t4[0])
        assert_array_equal(t1[1], t2[1])
        assert_array_equal(t2[1], t3[1])
        assert_array_equal(t3[1], t4[1])


def test_stratified_shuffle_split_init():
    y = np.asarray([0, 1, 1, 1, 2, 2, 2])
    # Check that error is raised if there is a class with only one sample
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 0.2)

    # Check that error is raised if the test set size is smaller than n_classes
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 2)
    # Check that error is raised if the train set size is smaller than
    # n_classes
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 3, 2)

    y = np.asarray([0, 0, 0, 1, 1, 1, 2, 2, 2])
    # Check that errors are raised if there is not enough samples
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 0.5, 0.6)
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 8, 0.6)
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, 3, 0.6, 8)

    # Train size or test size too small
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, train_size=2)
    assert_raises(ValueError, cval.StratifiedShuffleSplit, y, test_size=2)


def test_stratified_shuffle_split_iter():
    ys = [np.array([1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3]),
          np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]),
          np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]),
          np.array([1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4]),
          np.array([-1] * 800 + [1] * 50)
          ]

    for y in ys:
        sss = cval.StratifiedShuffleSplit(y, 6, test_size=0.33,
                                          random_state=0)
        for train, test in sss:
            assert_array_equal(np.unique(y[train]), np.unique(y[test]))
            # Checks if folds keep classes proportions
            p_train = (np.bincount(np.unique(y[train], return_inverse=True)[1])
                       / float(len(y[train])))
            p_test = (np.bincount(np.unique(y[test], return_inverse=True)[1])
                      / float(len(y[test])))
            assert_array_almost_equal(p_train, p_test, 1)
            assert_equal(y[train].size + y[test].size, y.size)
            assert_array_equal(np.lib.arraysetops.intersect1d(train, test), [])


@ignore_warnings
def test_stratified_shuffle_split_iter_no_indices():
    y = np.asarray([0, 1, 2] * 10)

    sss1 = cval.StratifiedShuffleSplit(y, indices=False, random_state=0)
    train_mask, test_mask = next(iter(sss1))

    sss2 = cval.StratifiedShuffleSplit(y, indices=True, random_state=0)
    train_indices, test_indices = next(iter(sss2))

    assert_array_equal(sorted(test_indices), np.where(test_mask)[0])


def test_stratified_shuffle_split_even():
    # Test the StratifiedShuffleSplit, indices are drawn with a
    # equal chance
    n_folds = 5
    n_iter = 1000

    def assert_counts_are_ok(idx_counts, p):
        # Here we test that the distribution of the counts
        # per index is close enough to a binomial
        threshold = 0.05 / n_splits
        bf = stats.binom(n_splits, p)
        for count in idx_counts:
            p = bf.pmf(count)
            assert_true(p > threshold,
                        "An index is not drawn with chance corresponding "
                        "to even draws")

    for n_samples in (6, 22):
        labels = np.array((n_samples // 2) * [0, 1])
        splits = cval.StratifiedShuffleSplit(labels, n_iter=n_iter,
                                             test_size=1./n_folds,
                                             random_state=0)

        train_counts = [0] * n_samples
        test_counts = [0] * n_samples
        n_splits = 0
        for train, test in splits:
            n_splits += 1
            for counter, ids in [(train_counts, train), (test_counts, test)]:
                for id in ids:
                    counter[id] += 1
        assert_equal(n_splits, n_iter)

        assert_equal(len(train), splits.n_train)
        assert_equal(len(test), splits.n_test)
        assert_equal(len(set(train).intersection(test)), 0)

        label_counts = np.unique(labels)
        assert_equal(splits.test_size, 1.0 / n_folds)
        assert_equal(splits.n_train + splits.n_test, len(labels))
        assert_equal(len(label_counts), 2)
        ex_test_p = float(splits.n_test) / n_samples
        ex_train_p = float(splits.n_train) / n_samples

        assert_counts_are_ok(train_counts, ex_train_p)
        assert_counts_are_ok(test_counts, ex_test_p)


def test_leave_label_out_changing_labels():
    """Check that LeaveOneLabelOut and LeavePLabelOut work normally if
    the labels variable is changed before calling __iter__"""
    labels = np.array([0, 1, 2, 1, 1, 2, 0, 0])
    labels_changing = np.array(labels, copy=True)
    lolo = cval.LeaveOneLabelOut(labels)
    lolo_changing = cval.LeaveOneLabelOut(labels_changing)
    lplo = cval.LeavePLabelOut(labels, p=2)
    lplo_changing = cval.LeavePLabelOut(labels_changing, p=2)
    labels_changing[:] = 0
    for llo, llo_changing in [(lolo, lolo_changing), (lplo, lplo_changing)]:
        for (train, test), (train_chan, test_chan) in zip(llo, llo_changing):
            assert_array_equal(train, train_chan)
            assert_array_equal(test, test_chan)


def test_cross_val_score():
    clf = MockClassifier()
    for a in range(-10, 10):
        clf.a = a
        # Smoke test
        scores = cval.cross_val_score(clf, X, y)
        assert_array_equal(scores, clf.score(X, y))

        # test with multioutput y
        scores = cval.cross_val_score(clf, X_sparse, X)
        assert_array_equal(scores, clf.score(X_sparse, X))

        scores = cval.cross_val_score(clf, X_sparse, y)
        assert_array_equal(scores, clf.score(X_sparse, y))

        # test with multioutput y
        scores = cval.cross_val_score(clf, X_sparse, X)
        assert_array_equal(scores, clf.score(X_sparse, X))

    # test with X and y as list
    list_check = lambda x: isinstance(x, list)
    clf = CheckingClassifier(check_X=list_check)
    scores = cval.cross_val_score(clf, X.tolist(), y.tolist())

    clf = CheckingClassifier(check_y=list_check)
    scores = cval.cross_val_score(clf, X, y.tolist())

    assert_raises(ValueError, cval.cross_val_score, clf, X, y,
                  scoring="sklearn")

    # test with 3d X and
    X_3d = X[:, :, np.newaxis]
    clf = MockClassifier(allow_nd=True)
    scores = cval.cross_val_score(clf, X_3d, y)

    clf = MockClassifier(allow_nd=False)
    assert_raises(ValueError, cval.cross_val_score, clf, X_3d, y)


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
        X_df, y_ser = InputFeatureType(X), TargetType(y)
        check_df = lambda x: isinstance(x, InputFeatureType)
        check_series = lambda x: isinstance(x, TargetType)
        clf = CheckingClassifier(check_X=check_df, check_y=check_series)
        cval.cross_val_score(clf, X_df, y_ser)


def test_cross_val_score_mask():
    # test that cross_val_score works with boolean masks
    svm = SVC(kernel="linear")
    iris = load_iris()
    X, y = iris.data, iris.target
    cv_indices = cval.KFold(len(y), 5, indices=True)
    scores_indices = cval.cross_val_score(svm, X, y, cv=cv_indices)
    cv_masks = cval.KFold(len(y), 5, indices=False)
    scores_masks = cval.cross_val_score(svm, X, y, cv=cv_masks)
    assert_array_equal(scores_indices, scores_masks)


def test_cross_val_score_precomputed():
    # test for svm with precomputed kernel
    svm = SVC(kernel="precomputed")
    iris = load_iris()
    X, y = iris.data, iris.target
    linear_kernel = np.dot(X, X.T)
    score_precomputed = cval.cross_val_score(svm, linear_kernel, y)
    svm = SVC(kernel="linear")
    score_linear = cval.cross_val_score(svm, X, y)
    assert_array_equal(score_precomputed, score_linear)

    # Error raised for non-square X
    svm = SVC(kernel="precomputed")
    assert_raises(ValueError, cval.cross_val_score, svm, X, y)

    # test error is raised when the precomputed kernel is not array-like
    # or sparse
    assert_raises(ValueError, cval.cross_val_score, svm,
                  linear_kernel.tolist(), y)


def test_cross_val_score_fit_params():
    clf = MockClassifier()
    n_samples = X.shape[0]
    n_classes = len(np.unique(y))
    fit_params = {'sample_weight': np.ones(n_samples),
                  'class_prior': np.ones(n_classes) / n_classes}
    cval.cross_val_score(clf, X, y, fit_params=fit_params)


def test_cross_val_score_score_func():
    clf = MockClassifier()
    _score_func_args = []

    def score_func(y_test, y_predict):
        _score_func_args.append((y_test, y_predict))
        return 1.0

    with warnings.catch_warnings(record=True):
        scoring = make_scorer(score_func)
        score = cval.cross_val_score(clf, X, y, scoring=scoring)
    assert_array_equal(score, [1.0, 1.0, 1.0])
    assert len(_score_func_args) == 3


def test_cross_val_score_errors():
    class BrokenEstimator:
        pass

    assert_raises(TypeError, cval.cross_val_score, BrokenEstimator(), X)


def test_train_test_split_errors():
    assert_raises(ValueError, cval.train_test_split)
    assert_raises(ValueError, cval.train_test_split, range(3), train_size=1.1)
    assert_raises(ValueError, cval.train_test_split, range(3), test_size=0.6,
                  train_size=0.6)
    assert_raises(ValueError, cval.train_test_split, range(3),
                  test_size=np.float32(0.6), train_size=np.float32(0.6))
    assert_raises(ValueError, cval.train_test_split, range(3),
                  test_size="wrong_type")
    assert_raises(ValueError, cval.train_test_split, range(3), test_size=2,
                  train_size=4)
    assert_raises(TypeError, cval.train_test_split, range(3),
                  some_argument=1.1)
    assert_raises(ValueError, cval.train_test_split, range(3), range(42))


def test_train_test_split():
    X = np.arange(100).reshape((10, 10))
    X_s = coo_matrix(X)
    y = range(10)
    split = cval.train_test_split(X, X_s, y, force_arrays=True)
    X_train, X_test, X_s_train, X_s_test, y_train, y_test = split
    assert_array_equal(X_train, X_s_train.toarray())
    assert_array_equal(X_test, X_s_test.toarray())
    assert_array_equal(X_train[:, 0], y_train * 10)
    assert_array_equal(X_test[:, 0], y_test * 10)
    split = cval.train_test_split(X, y, test_size=None, train_size=.5)
    X_train, X_test, y_train, y_test = split
    assert_equal(len(y_test), len(y_train))

    split = cval.train_test_split(X, X_s, y)
    X_train, X_test, X_s_train, X_s_test, y_train, y_test = split
    assert_true(isinstance(y_train, list))
    assert_true(isinstance(y_test, list))


def train_test_split_pandas():
    # check cross_val_score doesn't destroy pandas dataframe
    types = [MockDataFrame]
    try:
        from pandas import DataFrame
        types.append(DataFrame)
    except ImportError:
        pass
    for InputFeatureType in types:
        # X dataframe
        X_df = InputFeatureType(X)
        X_train, X_test = cval.train_test_split(X_df)
        assert_true(isinstance(X_train, InputFeatureType))
        assert_true(isinstance(X_test, InputFeatureType))


def train_test_split_mock_pandas():
    # X mock dataframe
    X_df = MockDataFrame(X)
    X_train, X_test = cval.train_test_split(X_df)
    assert_true(isinstance(X_train, MockDataFrame))
    assert_true(isinstance(X_test, MockDataFrame))
    X_train_array, X_test_array = cval.train_test_split(X_df, force_arrays=True)
    assert_true(isinstance(X_train_array, np.ndarray))
    assert_true(isinstance(X_test_array, np.ndarray))


def test_cross_val_score_with_score_func_classification():
    iris = load_iris()
    clf = SVC(kernel='linear')

    # Default score (should be the accuracy score)
    scores = cval.cross_val_score(clf, iris.data, iris.target, cv=5)
    assert_array_almost_equal(scores, [0.97, 1., 0.97, 0.97, 1.], 2)

    # Correct classification score (aka. zero / one score) - should be the
    # same as the default estimator score
    zo_scores = cval.cross_val_score(clf, iris.data, iris.target,
                                     scoring="accuracy", cv=5)
    assert_array_almost_equal(zo_scores, [0.97, 1., 0.97, 0.97, 1.], 2)

    # F1 score (class are balanced so f1_score should be equal to zero/one
    # score
    f1_scores = cval.cross_val_score(clf, iris.data, iris.target,
                                     scoring="f1", cv=5)
    assert_array_almost_equal(f1_scores, [0.97, 1., 0.97, 0.97, 1.], 2)


def test_cross_val_score_with_score_func_regression():
    X, y = make_regression(n_samples=30, n_features=20, n_informative=5,
                           random_state=0)
    reg = Ridge()

    # Default score of the Ridge regression estimator
    scores = cval.cross_val_score(reg, X, y, cv=5)
    assert_array_almost_equal(scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # R2 score (aka. determination coefficient) - should be the
    # same as the default estimator score
    r2_scores = cval.cross_val_score(reg, X, y, scoring="r2", cv=5)
    assert_array_almost_equal(r2_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)

    # Mean squared error; this is a loss function, so "scores" are negative
    mse_scores = cval.cross_val_score(reg, X, y, cv=5,
                                      scoring="mean_squared_error")
    expected_mse = np.array([-763.07, -553.16, -274.38, -273.26, -1681.99])
    assert_array_almost_equal(mse_scores, expected_mse, 2)

    # Explained variance
    scoring = make_scorer(explained_variance_score)
    ev_scores = cval.cross_val_score(reg, X, y, cv=5, scoring=scoring)
    assert_array_almost_equal(ev_scores, [0.94, 0.97, 0.97, 0.99, 0.92], 2)


def test_permutation_score():
    iris = load_iris()
    X = iris.data
    X_sparse = coo_matrix(X)
    y = iris.target
    svm = SVC(kernel='linear')
    cv = cval.StratifiedKFold(y, 2)

    score, scores, pvalue = cval.permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy")
    assert_greater(score, 0.9)
    assert_almost_equal(pvalue, 0.0, 1)

    score_label, _, pvalue_label = cval.permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy",
        labels=np.ones(y.size), random_state=0)
    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # check that we obtain the same results with a sparse representation
    svm_sparse = SVC(kernel='linear')
    cv_sparse = cval.StratifiedKFold(y, 2)
    score_label, _, pvalue_label = cval.permutation_test_score(
        svm_sparse, X_sparse, y, n_permutations=30, cv=cv_sparse,
        scoring="accuracy", labels=np.ones(y.size), random_state=0)

    assert_true(score_label == score)
    assert_true(pvalue_label == pvalue)

    # test with custom scoring object
    def custom_score(y_true, y_pred):
        return (((y_true == y_pred).sum() - (y_true != y_pred).sum())
                / y_true.shape[0])

    scorer = make_scorer(custom_score)
    score, _, pvalue = cval.permutation_test_score(
        svm, X, y, n_permutations=100, scoring=scorer, cv=cv, random_state=0)
    assert_almost_equal(score, .93, 2)
    assert_almost_equal(pvalue, 0.01, 3)

    # set random y
    y = np.mod(np.arange(len(y)), 3)

    score, scores, pvalue = cval.permutation_test_score(
        svm, X, y, n_permutations=30, cv=cv, scoring="accuracy")

    assert_less(score, 0.5)
    assert_greater(pvalue, 0.2)


def test_cross_val_generator_with_mask():
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = assert_warns(DeprecationWarning, cval.LeaveOneOut,
                       4, indices=False)
    lpo = assert_warns(DeprecationWarning, cval.LeavePOut,
                       4, 2, indices=False)
    kf = assert_warns(DeprecationWarning, cval.KFold,
                      4, 2, indices=False)
    skf = assert_warns(DeprecationWarning, cval.StratifiedKFold,
                       y, 2, indices=False)
    lolo = assert_warns(DeprecationWarning, cval.LeaveOneLabelOut,
                        labels, indices=False)
    lopo = assert_warns(DeprecationWarning, cval.LeavePLabelOut,
                        labels, 2, indices=False)
    ss = assert_warns(DeprecationWarning, cval.ShuffleSplit,
                      4, indices=False)
    for cv in [loo, lpo, kf, skf, lolo, lopo, ss]:
        for train, test in cv:
            assert_equal(np.asarray(train).dtype.kind, 'b')
            assert_equal(np.asarray(train).dtype.kind, 'b')
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]


def test_cross_val_generator_with_indices():
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    # explicitly passing indices value is deprecated
    loo = assert_warns(DeprecationWarning, cval.LeaveOneOut,
                       4, indices=True)
    lpo = assert_warns(DeprecationWarning, cval.LeavePOut,
                       4, 2, indices=True)
    kf = assert_warns(DeprecationWarning, cval.KFold,
                      4, 2, indices=True)
    skf = assert_warns(DeprecationWarning, cval.StratifiedKFold,
                       y, 2, indices=True)
    lolo = assert_warns(DeprecationWarning, cval.LeaveOneLabelOut,
                        labels, indices=True)
    lopo = assert_warns(DeprecationWarning, cval.LeavePLabelOut,
                        labels, 2, indices=True)
    # Bootstrap as a cross-validation is deprecated
    b = assert_warns(DeprecationWarning, cval.Bootstrap, 2)
    ss = assert_warns(DeprecationWarning, cval.ShuffleSplit,
                      2, indices=True)
    for cv in [loo, lpo, kf, skf, lolo, lopo, b, ss]:
        for train, test in cv:
            assert_not_equal(np.asarray(train).dtype.kind, 'b')
            assert_not_equal(np.asarray(train).dtype.kind, 'b')
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]


@ignore_warnings
def test_cross_val_generator_with_default_indices():
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cval.LeaveOneOut(4)
    lpo = cval.LeavePOut(4, 2)
    kf = cval.KFold(4, 2)
    skf = cval.StratifiedKFold(y, 2)
    lolo = cval.LeaveOneLabelOut(labels)
    lopo = cval.LeavePLabelOut(labels, 2)
    b = cval.Bootstrap(2)  # only in index mode
    ss = cval.ShuffleSplit(2)
    for cv in [loo, lpo, kf, skf, lolo, lopo, b, ss]:
        for train, test in cv:
            assert_not_equal(np.asarray(train).dtype.kind, 'b')
            assert_not_equal(np.asarray(train).dtype.kind, 'b')
            X_train, X_test = X[train], X[test]
            y_train, y_test = y[train], y[test]


@ignore_warnings
def test_cross_val_generator_mask_indices_same():
    # Test that the cross validation generators return the same results when
    # indices=True and when indices=False
    y = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2])
    labels = np.array([1, 1, 2, 3, 3, 3, 4])

    loo_mask = cval.LeaveOneOut(5, indices=False)
    loo_ind = cval.LeaveOneOut(5, indices=True)
    lpo_mask = cval.LeavePOut(10, 2, indices=False)
    lpo_ind = cval.LeavePOut(10, 2, indices=True)
    kf_mask = cval.KFold(10, 5, indices=False, shuffle=True, random_state=1)
    kf_ind = cval.KFold(10, 5, indices=True, shuffle=True, random_state=1)
    skf_mask = cval.StratifiedKFold(y, 3, indices=False)
    skf_ind = cval.StratifiedKFold(y, 3, indices=True)
    lolo_mask = cval.LeaveOneLabelOut(labels, indices=False)
    lolo_ind = cval.LeaveOneLabelOut(labels, indices=True)
    lopo_mask = cval.LeavePLabelOut(labels, 2, indices=False)
    lopo_ind = cval.LeavePLabelOut(labels, 2, indices=True)

    for cv_mask, cv_ind in [(loo_mask, loo_ind), (lpo_mask, lpo_ind),
                            (kf_mask, kf_ind), (skf_mask, skf_ind),
                            (lolo_mask, lolo_ind), (lopo_mask, lopo_ind)]:
        for (train_mask, test_mask), (train_ind, test_ind) in \
                zip(cv_mask, cv_ind):
            assert_array_equal(np.where(train_mask)[0], train_ind)
            assert_array_equal(np.where(test_mask)[0], test_ind)


@ignore_warnings
def test_bootstrap_errors():
    assert_raises(ValueError, cval.Bootstrap, 10, train_size=100)
    assert_raises(ValueError, cval.Bootstrap, 10, test_size=100)
    assert_raises(ValueError, cval.Bootstrap, 10, train_size=1.1)
    assert_raises(ValueError, cval.Bootstrap, 10, test_size=1.1)
    assert_raises(ValueError, cval.Bootstrap, 10, train_size=0.6,
                  test_size=0.5)


@ignore_warnings
def test_bootstrap_test_sizes():
    assert_equal(cval.Bootstrap(10, test_size=0.2).test_size, 2)
    assert_equal(cval.Bootstrap(10, test_size=2).test_size, 2)
    assert_equal(cval.Bootstrap(10, test_size=None).test_size, 5)


def test_shufflesplit_errors():
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=2.0)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=1.0)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=0.1,
                  train_size=0.95)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=11)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=10)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=8, train_size=3)
    assert_raises(ValueError, cval.ShuffleSplit, 10, train_size=1j)
    assert_raises(ValueError, cval.ShuffleSplit, 10, test_size=None,
                  train_size=None)


def test_shufflesplit_reproducible():
    # Check that iterating twice on the ShuffleSplit gives the same
    # sequence of train-test when the random_state is given
    ss = cval.ShuffleSplit(10, random_state=21)
    assert_array_equal(list(a for a, b in ss), list(a for a, b in ss))


@ignore_warnings
def test_cross_indices_exception():
    X = coo_matrix(np.array([[1, 2], [3, 4], [5, 6], [7, 8]]))
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = cval.LeaveOneOut(4, indices=False)
    lpo = cval.LeavePOut(4, 2, indices=False)
    kf = cval.KFold(4, 2, indices=False)
    skf = cval.StratifiedKFold(y, 2, indices=False)
    lolo = cval.LeaveOneLabelOut(labels, indices=False)
    lopo = cval.LeavePLabelOut(labels, 2, indices=False)

    assert_raises(ValueError, cval.check_cv, loo, X, y)
    assert_raises(ValueError, cval.check_cv, lpo, X, y)
    assert_raises(ValueError, cval.check_cv, kf, X, y)
    assert_raises(ValueError, cval.check_cv, skf, X, y)
    assert_raises(ValueError, cval.check_cv, lolo, X, y)
    assert_raises(ValueError, cval.check_cv, lopo, X, y)


def test_safe_split_with_precomputed_kernel():
    clf = SVC()
    clfp = SVC(kernel="precomputed")

    iris = load_iris()
    X, y = iris.data, iris.target
    K = np.dot(X, X.T)

    cv = cval.ShuffleSplit(X.shape[0], test_size=0.25, random_state=0)
    tr, te = list(cv)[0]

    X_tr, y_tr = cval._safe_split(clf, X, y, tr)
    K_tr, y_tr2 = cval._safe_split(clfp, K, y, tr)
    assert_array_almost_equal(K_tr, np.dot(X_tr, X_tr.T))

    X_te, y_te = cval._safe_split(clf, X, y, te, tr)
    K_te, y_te2 = cval._safe_split(clfp, K, y, te, tr)
    assert_array_almost_equal(K_te, np.dot(X_te, X_tr.T))


def test_cross_val_score_allow_nans():
    # Check that cross_val_score allows input data with NaNs
    X = np.arange(200, dtype=np.float64).reshape(10, -1)
    X[2, :] = np.nan
    y = np.repeat([0, 1], X.shape[0]/2)
    p = Pipeline([
        ('imputer', Imputer(strategy='mean', missing_values='NaN')),
        ('classifier', MockClassifier()),
    ])
    cval.cross_val_score(p, X, y, cv=5)


def test_train_test_split_allow_nans():
    # Check that train_test_split allows input data with NaNs
    X = np.arange(200, dtype=np.float64).reshape(10, -1)
    X[2, :] = np.nan
    y = np.repeat([0, 1], X.shape[0]/2)
    cval.train_test_split(X, y, test_size=0.2, random_state=42)


def test_permutation_test_score_allow_nans():
    # Check that permutation_test_score allows input data with NaNs
    X = np.arange(200, dtype=np.float64).reshape(10, -1)
    X[2, :] = np.nan
    y = np.repeat([0, 1], X.shape[0]/2)
    p = Pipeline([
        ('imputer', Imputer(strategy='mean', missing_values='NaN')),
        ('classifier', MockClassifier()),
    ])
    cval.permutation_test_score(p, X, y, cv=5)


def test_check_cv_return_types():
    X = np.ones((9, 2))
    cv = cval._check_cv(3, X, classifier=False)
    assert_true(isinstance(cv, cval.KFold))

    y_binary = np.array([0, 1, 0, 1, 0, 0, 1, 1, 1])
    cv = cval._check_cv(3, X, y_binary, classifier=True)
    assert_true(isinstance(cv, cval.StratifiedKFold))

    y_multiclass = np.array([0, 1, 0, 1, 2, 1, 2, 0, 2])
    cv = cval._check_cv(3, X, y_multiclass, classifier=True)
    assert_true(isinstance(cv, cval.StratifiedKFold))

    X = np.ones((5, 2))
    y_seq_of_seqs = [[], [1, 2], [3], [0, 1, 3], [2]]
    cv = cval._check_cv(3, X, y_seq_of_seqs, classifier=True)
    assert_true(isinstance(cv, cval.KFold))

    y_indicator_matrix = LabelBinarizer().fit_transform(y_seq_of_seqs)
    cv = cval._check_cv(3, X, y_indicator_matrix, classifier=True)
    assert_true(isinstance(cv, cval.KFold))

    y_multioutput = np.array([[1, 2], [0, 3], [0, 0], [3, 1], [2, 0]])
    cv = cval._check_cv(3, X, y_multioutput, classifier=True)
    assert_true(isinstance(cv, cval.KFold))


def test_cross_val_score_multilabel():
    X = np.array([[-3, 4], [2, 4], [3, 3], [0, 2], [-3, 1],
                  [-2, 1], [0, 0], [-2, -1], [-1, -2], [1, -2]])
    y = np.array([[1, 1], [0, 1], [0, 1], [0, 1], [1, 1],
                  [0, 1], [1, 0], [1, 1], [1, 0], [0, 0]])
    clf = KNeighborsClassifier(n_neighbors=1)
    scoring_micro = make_scorer(precision_score, average='micro')
    scoring_macro = make_scorer(precision_score, average='macro')
    scoring_samples = make_scorer(precision_score, average='samples')
    score_micro = cval.cross_val_score(clf, X, y, scoring=scoring_micro, cv=5)
    score_macro = cval.cross_val_score(clf, X, y, scoring=scoring_macro, cv=5)
    score_samples = cval.cross_val_score(clf, X, y,
                                         scoring=scoring_samples, cv=5)
    assert_almost_equal(score_micro, [1, 1/2, 3/4, 1/2, 1/3])
    assert_almost_equal(score_macro, [1, 1/2, 3/4, 1/2, 1/4])
    assert_almost_equal(score_samples, [1, 1/2, 3/4, 1/2, 1/4])
