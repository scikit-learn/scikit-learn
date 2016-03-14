"""Test the split module"""
from __future__ import division
import warnings

import numpy as np
from scipy.sparse import coo_matrix
from scipy import stats
from scipy.misc import comb
from itertools import combinations
from sklearn.utils.fixes import combinations_with_replacement

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.validation import _num_samples
from sklearn.utils.mocking import MockDataFrame

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import LabelKFold
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import LeaveOneLabelOut
from sklearn.model_selection import LeavePOut
from sklearn.model_selection import LeavePLabelOut
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import LabelShuffleSplit
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import PredefinedSplit
from sklearn.model_selection import check_cv
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

from sklearn.linear_model import Ridge

from sklearn.model_selection._split import _safe_split
from sklearn.model_selection._split import _validate_shuffle_split
from sklearn.model_selection._split import _CVIterableWrapper
from sklearn.model_selection._split import _build_repr

from sklearn.datasets import load_digits
from sklearn.datasets import load_iris
from sklearn.datasets import make_classification

from sklearn.externals import six
from sklearn.externals.six.moves import zip

from sklearn.svm import SVC

X = np.ones(10)
y = np.arange(10) // 2
P_sparse = coo_matrix(np.eye(5))
iris = load_iris()
digits = load_digits()


class MockClassifier(object):
    """Dummy classifier to test the cross-validation"""

    def __init__(self, a=0, allow_nd=False):
        self.a = a
        self.allow_nd = allow_nd

    def fit(self, X, Y=None, sample_weight=None, class_prior=None,
            sparse_sample_weight=None, sparse_param=None, dummy_int=None,
            dummy_str=None, dummy_obj=None, callback=None):
        """The dummy arguments are to test that this fit function can
        accept non-array arguments through cross-validation, such as:
            - int
            - str (this is actually array-like)
            - object
            - function
        """
        self.dummy_int = dummy_int
        self.dummy_str = dummy_str
        self.dummy_obj = dummy_obj
        if callback is not None:
            callback(self)

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
        if sparse_sample_weight is not None:
            fmt = ('MockClassifier extra fit_param sparse_sample_weight'
                   '.shape[0] is {0}, should be {1}')
            assert_true(sparse_sample_weight.shape[0] == X.shape[0],
                        fmt.format(sparse_sample_weight.shape[0], X.shape[0]))
        if sparse_param is not None:
            fmt = ('MockClassifier extra fit_param sparse_param.shape '
                   'is ({0}, {1}), should be ({2}, {3})')
            assert_true(sparse_param.shape == P_sparse.shape,
                        fmt.format(sparse_param.shape[0],
                                   sparse_param.shape[1],
                                   P_sparse.shape[0], P_sparse.shape[1]))
        return self

    def predict(self, T):
        if self.allow_nd:
            T = T.reshape(len(T), -1)
        return T[:, 0]

    def score(self, X=None, Y=None):
        return 1. / (1 + np.abs(self.a))

    def get_params(self, deep=False):
        return {'a': self.a, 'allow_nd': self.allow_nd}


@ignore_warnings
def test_cross_validator_with_default_params():
    n_samples = 4
    n_unique_labels = 4
    n_folds = 2
    p = 2
    n_iter = 10  # (the default value)

    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    X_1d = np.array([1, 2, 3, 4])
    y = np.array([1, 1, 2, 2])
    labels = np.array([1, 2, 3, 4])
    loo = LeaveOneOut()
    lpo = LeavePOut(p)
    kf = KFold(n_folds)
    skf = StratifiedKFold(n_folds)
    lolo = LeaveOneLabelOut()
    lopo = LeavePLabelOut(p)
    ss = ShuffleSplit(random_state=0)
    ps = PredefinedSplit([1, 1, 2, 2])  # n_splits = np of unique folds = 2

    loo_repr = "LeaveOneOut()"
    lpo_repr = "LeavePOut(p=2)"
    kf_repr = "KFold(n_folds=2, random_state=None, shuffle=False)"
    skf_repr = "StratifiedKFold(n_folds=2, random_state=None, shuffle=False)"
    lolo_repr = "LeaveOneLabelOut()"
    lopo_repr = "LeavePLabelOut(n_labels=2)"
    ss_repr = ("ShuffleSplit(n_iter=10, random_state=0, test_size=0.1, "
               "train_size=None)")
    ps_repr = "PredefinedSplit(test_fold=array([1, 1, 2, 2]))"

    n_splits = [n_samples, comb(n_samples, p), n_folds, n_folds,
                n_unique_labels, comb(n_unique_labels, p), n_iter, 2]

    for i, (cv, cv_repr) in enumerate(zip(
            [loo, lpo, kf, skf, lolo, lopo, ss, ps],
            [loo_repr, lpo_repr, kf_repr, skf_repr, lolo_repr, lopo_repr,
             ss_repr, ps_repr])):
        # Test if get_n_splits works correctly
        assert_equal(n_splits[i], cv.get_n_splits(X, y, labels))

        # Test if the cross-validator works as expected even if
        # the data is 1d
        np.testing.assert_equal(list(cv.split(X, y, labels)),
                                list(cv.split(X_1d, y, labels)))
        # Test that train, test indices returned are integers
        for train, test in cv.split(X, y, labels):
            assert_equal(np.asarray(train).dtype.kind, 'i')
            assert_equal(np.asarray(train).dtype.kind, 'i')

        # Test if the repr works without any errors
        assert_equal(cv_repr, repr(cv))


def check_valid_split(train, test, n_samples=None):
    # Use python sets to get more informative assertion failure messages
    train, test = set(train), set(test)

    # Train and test split should not overlap
    assert_equal(train.intersection(test), set())

    if n_samples is not None:
        # Check that the union of train an test split cover all the indices
        assert_equal(train.union(test), set(range(n_samples)))


def check_cv_coverage(cv, X, y, labels, expected_n_iter=None):
    n_samples = _num_samples(X)
    # Check that a all the samples appear at least once in a test fold
    if expected_n_iter is not None:
        assert_equal(cv.get_n_splits(X, y, labels), expected_n_iter)
    else:
        expected_n_iter = cv.get_n_splits(X, y, labels)

    collected_test_samples = set()
    iterations = 0
    for train, test in cv.split(X, y, labels):
        check_valid_split(train, test, n_samples=n_samples)
        iterations += 1
        collected_test_samples.update(test)

    # Check that the accumulated test samples cover the whole dataset
    assert_equal(iterations, expected_n_iter)
    if n_samples is not None:
        assert_equal(collected_test_samples, set(range(n_samples)))


def test_kfold_valueerrors():
    X1 = np.array([[1, 2], [3, 4], [5, 6]])
    X2 = np.array([[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]])
    # Check that errors are raised if there is not enough samples
    assert_raises(ValueError, next, KFold(4).split(X1))

    # Check that a warning is raised if the least populated class has too few
    # members.
    y = np.array([3, 3, -1, -1, 3])

    skf_3 = StratifiedKFold(3)
    assert_warns_message(Warning, "The least populated class",
                         next, skf_3.split(X2, y))

    # Check that despite the warning the folds are still computed even
    # though all the classes are not necessarily represented at on each
    # side of the split at each split
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        check_cv_coverage(skf_3, X2, y, labels=None, expected_n_iter=3)

    # Check that errors are raised if all n_labels for individual
    # classes are less than n_folds.
    y = np.array([3, 3, -1, -1, 2])

    assert_raises(ValueError, next, skf_3.split(X2, y))

    # Error when number of folds is <= 1
    assert_raises(ValueError, KFold, 0)
    assert_raises(ValueError, KFold, 1)
    error_string = ("k-fold cross-validation requires at least one"
                    " train/test split")
    assert_raise_message(ValueError, error_string,
                         StratifiedKFold, 0)
    assert_raise_message(ValueError, error_string,
                         StratifiedKFold, 1)

    # When n_folds is not integer:
    assert_raises(ValueError, KFold, 1.5)
    assert_raises(ValueError, KFold, 2.0)
    assert_raises(ValueError, StratifiedKFold, 1.5)
    assert_raises(ValueError, StratifiedKFold, 2.0)

    # When shuffle is not  a bool:
    assert_raises(TypeError, KFold, n_folds=4, shuffle=None)


def test_kfold_indices():
    # Check all indices are returned in the test folds
    X1 = np.ones(18)
    kf = KFold(3)
    check_cv_coverage(kf, X1, y=None, labels=None, expected_n_iter=3)

    # Check all indices are returned in the test folds even when equal-sized
    # folds are not possible
    X2 = np.ones(17)
    kf = KFold(3)
    check_cv_coverage(kf, X2, y=None, labels=None, expected_n_iter=3)

    # Check if get_n_splits returns the number of folds
    assert_equal(5, KFold(5).get_n_splits(X2))


def test_kfold_no_shuffle():
    # Manually check that KFold preserves the data ordering on toy datasets
    X2 = [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]

    splits = KFold(2).split(X2[:-1])
    train, test = next(splits)
    assert_array_equal(test, [0, 1])
    assert_array_equal(train, [2, 3])

    train, test = next(splits)
    assert_array_equal(test, [2, 3])
    assert_array_equal(train, [0, 1])

    splits = KFold(2).split(X2)
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
    X, y = np.ones(4), [1, 1, 0, 0]
    splits = StratifiedKFold(2).split(X, y)
    train, test = next(splits)
    assert_array_equal(test, [0, 2])
    assert_array_equal(train, [1, 3])

    train, test = next(splits)
    assert_array_equal(test, [1, 3])
    assert_array_equal(train, [0, 2])

    X, y = np.ones(7), [1, 1, 1, 0, 0, 0, 0]
    splits = StratifiedKFold(2).split(X, y)
    train, test = next(splits)
    assert_array_equal(test, [0, 1, 3, 4])
    assert_array_equal(train, [2, 5, 6])

    train, test = next(splits)
    assert_array_equal(test, [2, 5, 6])
    assert_array_equal(train, [0, 1, 3, 4])

    # Check if get_n_splits returns the number of folds
    assert_equal(5, StratifiedKFold(5).get_n_splits(X, y))


def test_stratified_kfold_ratios():
    # Check that stratified kfold preserves class ratios in individual splits
    # Repeat with shuffling turned off and on
    n_samples = 1000
    X = np.ones(n_samples)
    y = np.array([4] * int(0.10 * n_samples) +
                 [0] * int(0.89 * n_samples) +
                 [1] * int(0.01 * n_samples))

    for shuffle in (False, True):
        for train, test in StratifiedKFold(5, shuffle=shuffle).split(X, y):
            assert_almost_equal(np.sum(y[train] == 4) / len(train), 0.10, 2)
            assert_almost_equal(np.sum(y[train] == 0) / len(train), 0.89, 2)
            assert_almost_equal(np.sum(y[train] == 1) / len(train), 0.01, 2)
            assert_almost_equal(np.sum(y[test] == 4) / len(test), 0.10, 2)
            assert_almost_equal(np.sum(y[test] == 0) / len(test), 0.89, 2)
            assert_almost_equal(np.sum(y[test] == 1) / len(test), 0.01, 2)


def test_kfold_balance():
    # Check that KFold returns folds with balanced sizes
    for i in range(11, 17):
        kf = KFold(5).split(X=np.ones(i))
        sizes = []
        for _, test in kf:
            sizes.append(len(test))

        assert_true((np.max(sizes) - np.min(sizes)) <= 1)
        assert_equal(np.sum(sizes), i)


def test_stratifiedkfold_balance():
    # Check that KFold returns folds with balanced sizes (only when
    # stratification is possible)
    # Repeat with shuffling turned off and on
    X = np.ones(17)
    y = [0] * 3 + [1] * 14

    for shuffle in (True, False):
        cv = StratifiedKFold(3, shuffle=shuffle)
        for i in range(11, 17):
            skf = cv.split(X[:i], y[:i])
            sizes = []
            for _, test in skf:
                sizes.append(len(test))

            assert_true((np.max(sizes) - np.min(sizes)) <= 1)
            assert_equal(np.sum(sizes), i)


def test_shuffle_kfold():
    # Check the indices are shuffled properly
    kf = KFold(3)
    kf2 = KFold(3, shuffle=True, random_state=0)
    kf3 = KFold(3, shuffle=True, random_state=1)

    X = np.ones(300)

    all_folds = np.zeros(300)
    for (tr1, te1), (tr2, te2), (tr3, te3) in zip(
            kf.split(X), kf2.split(X), kf3.split(X)):
        for tr_a, tr_b in combinations((tr1, tr2, tr3), 2):
            # Assert that there is no complete overlap
            assert_not_equal(len(np.intersect1d(tr_a, tr_b)), len(tr1))

        # Set all test indices in successive iterations of kf2 to 1
        all_folds[te2] = 1

    # Check that all indices are returned in the different test folds
    assert_equal(sum(all_folds), 300)


def test_shuffle_kfold_stratifiedkfold_reproducibility():
    # Check that when the shuffle is True multiple split calls produce the
    # same split when random_state is set
    X = np.ones(15)  # Divisible by 3
    y = [0] * 7 + [1] * 8
    X2 = np.ones(16)  # Not divisible by 3
    y2 = [0] * 8 + [1] * 8

    kf = KFold(3, shuffle=True, random_state=0)
    skf = StratifiedKFold(3, shuffle=True, random_state=0)

    for cv in (kf, skf):
        np.testing.assert_equal(list(cv.split(X, y)), list(cv.split(X, y)))
        np.testing.assert_equal(list(cv.split(X2, y2)), list(cv.split(X2, y2)))

    kf = KFold(3, shuffle=True)
    skf = StratifiedKFold(3, shuffle=True)

    for cv in (kf, skf):
        for data in zip((X, X2), (y, y2)):
            try:
                np.testing.assert_equal(list(cv.split(*data)),
                                        list(cv.split(*data)))
            except AssertionError:
                pass
            else:
                raise AssertionError("The splits for data, %s, are same even "
                                     "when random state is not set" % data)


def test_shuffle_stratifiedkfold():
    # Check that shuffling is happening when requested, and for proper
    # sample coverage
    X_40 = np.ones(40)
    y = [0] * 20 + [1] * 20
    kf0 = StratifiedKFold(5, shuffle=True, random_state=0)
    kf1 = StratifiedKFold(5, shuffle=True, random_state=1)
    for (_, test0), (_, test1) in zip(kf0.split(X_40, y),
                                      kf1.split(X_40, y)):
        assert_not_equal(set(test0), set(test1))
    check_cv_coverage(kf0, X_40, y, labels=None, expected_n_iter=5)


def test_kfold_can_detect_dependent_samples_on_digits():  # see #2372
    # The digits samples are dependent: they are apparently grouped by authors
    # although we don't have any information on the groups segment locations
    # for this data. We can highlight this fact by computing k-fold cross-
    # validation with and without shuffling: we observe that the shuffling case
    # wrongly makes the IID assumption and is therefore too optimistic: it
    # estimates a much higher accuracy (around 0.93) than that the non
    # shuffling variant (around 0.81).

    X, y = digits.data[:600], digits.target[:600]
    model = SVC(C=10, gamma=0.005)

    n_folds = 3

    cv = KFold(n_folds=n_folds, shuffle=False)
    mean_score = cross_val_score(model, X, y, cv=cv).mean()
    assert_greater(0.92, mean_score)
    assert_greater(mean_score, 0.80)

    # Shuffling the data artificially breaks the dependency and hides the
    # overfitting of the model with regards to the writing style of the authors
    # by yielding a seriously overestimated score:

    cv = KFold(n_folds, shuffle=True, random_state=0)
    mean_score = cross_val_score(model, X, y, cv=cv).mean()
    assert_greater(mean_score, 0.92)

    cv = KFold(n_folds, shuffle=True, random_state=1)
    mean_score = cross_val_score(model, X, y, cv=cv).mean()
    assert_greater(mean_score, 0.92)

    # Similarly, StratifiedKFold should try to shuffle the data as little
    # as possible (while respecting the balanced class constraints)
    # and thus be able to detect the dependency by not overestimating
    # the CV score either. As the digits dataset is approximately balanced
    # the estimated mean score is close to the score measured with
    # non-shuffled KFold

    cv = StratifiedKFold(n_folds)
    mean_score = cross_val_score(model, X, y, cv=cv).mean()
    assert_greater(0.93, mean_score)
    assert_greater(mean_score, 0.80)


def test_shuffle_split():
    ss1 = ShuffleSplit(test_size=0.2, random_state=0).split(X)
    ss2 = ShuffleSplit(test_size=2, random_state=0).split(X)
    ss3 = ShuffleSplit(test_size=np.int32(2), random_state=0).split(X)
    for typ in six.integer_types:
        ss4 = ShuffleSplit(test_size=typ(2), random_state=0).split(X)
    for t1, t2, t3, t4 in zip(ss1, ss2, ss3, ss4):
        assert_array_equal(t1[0], t2[0])
        assert_array_equal(t2[0], t3[0])
        assert_array_equal(t3[0], t4[0])
        assert_array_equal(t1[1], t2[1])
        assert_array_equal(t2[1], t3[1])
        assert_array_equal(t3[1], t4[1])


def test_stratified_shuffle_split_init():
    X = np.arange(7)
    y = np.asarray([0, 1, 1, 1, 2, 2, 2])
    # Check that error is raised if there is a class with only one sample
    assert_raises(ValueError, next,
                  StratifiedShuffleSplit(3, 0.2).split(X, y))

    # Check that error is raised if the test set size is smaller than n_classes
    assert_raises(ValueError, next, StratifiedShuffleSplit(3, 2).split(X, y))
    # Check that error is raised if the train set size is smaller than
    # n_classes
    assert_raises(ValueError, next,
                  StratifiedShuffleSplit(3, 3, 2).split(X, y))

    X = np.arange(9)
    y = np.asarray([0, 0, 0, 1, 1, 1, 2, 2, 2])
    # Check that errors are raised if there is not enough samples
    assert_raises(ValueError, StratifiedShuffleSplit, 3, 0.5, 0.6)
    assert_raises(ValueError, next,
                  StratifiedShuffleSplit(3, 8, 0.6).split(X, y))
    assert_raises(ValueError, next,
                  StratifiedShuffleSplit(3, 0.6, 8).split(X, y))

    # Train size or test size too small
    assert_raises(ValueError, next,
                  StratifiedShuffleSplit(train_size=2).split(X, y))
    assert_raises(ValueError, next,
                  StratifiedShuffleSplit(test_size=2).split(X, y))


def test_stratified_shuffle_split_iter():
    ys = [np.array([1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3]),
          np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]),
          np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]),
          np.array([1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4]),
          np.array([-1] * 800 + [1] * 50)
          ]

    for y in ys:
        sss = StratifiedShuffleSplit(6, test_size=0.33,
                                     random_state=0).split(np.ones(len(y)), y)
        for train, test in sss:
            assert_array_equal(np.unique(y[train]), np.unique(y[test]))
            # Checks if folds keep classes proportions
            p_train = (np.bincount(np.unique(y[train],
                                   return_inverse=True)[1]) /
                       float(len(y[train])))
            p_test = (np.bincount(np.unique(y[test],
                                  return_inverse=True)[1]) /
                      float(len(y[test])))
            assert_array_almost_equal(p_train, p_test, 1)
            assert_equal(y[train].size + y[test].size, y.size)
            assert_array_equal(np.lib.arraysetops.intersect1d(train, test), [])


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
        splits = StratifiedShuffleSplit(n_iter=n_iter,
                                        test_size=1. / n_folds,
                                        random_state=0)

        train_counts = [0] * n_samples
        test_counts = [0] * n_samples
        n_splits = 0
        for train, test in splits.split(X=np.ones(n_samples), y=labels):
            n_splits += 1
            for counter, ids in [(train_counts, train), (test_counts, test)]:
                for id in ids:
                    counter[id] += 1
        assert_equal(n_splits, n_iter)

        n_train, n_test = _validate_shuffle_split(n_samples,
                                                  test_size=1./n_folds,
                                                  train_size=1.-(1./n_folds))

        assert_equal(len(train), n_train)
        assert_equal(len(test), n_test)
        assert_equal(len(set(train).intersection(test)), 0)

        label_counts = np.unique(labels)
        assert_equal(splits.test_size, 1.0 / n_folds)
        assert_equal(n_train + n_test, len(labels))
        assert_equal(len(label_counts), 2)
        ex_test_p = float(n_test) / n_samples
        ex_train_p = float(n_train) / n_samples

        assert_counts_are_ok(train_counts, ex_train_p)
        assert_counts_are_ok(test_counts, ex_test_p)


def test_stratified_shuffle_split_overlap_train_test_bug():
    # See https://github.com/scikit-learn/scikit-learn/issues/6121 for
    # the original bug report
    y = [0, 1, 2, 3] * 3 + [4, 5] * 5
    X = np.ones_like(y)

    splits = StratifiedShuffleSplit(n_iter=1,
                                    test_size=0.5, random_state=0)

    train, test = next(iter(splits.split(X=X, y=y)))

    assert_array_equal(np.intersect1d(train, test), [])


def test_predefinedsplit_with_kfold_split():
    # Check that PredefinedSplit can reproduce a split generated by Kfold.
    folds = -1 * np.ones(10)
    kf_train = []
    kf_test = []
    for i, (train_ind, test_ind) in enumerate(KFold(5, shuffle=True).split(X)):
        kf_train.append(train_ind)
        kf_test.append(test_ind)
        folds[test_ind] = i
    ps_train = []
    ps_test = []
    ps = PredefinedSplit(folds)
    # n_splits is simply the no of unique folds
    assert_equal(len(np.unique(folds)), ps.get_n_splits())
    for train_ind, test_ind in ps.split():
        ps_train.append(train_ind)
        ps_test.append(test_ind)
    assert_array_equal(ps_train, kf_train)
    assert_array_equal(ps_test, kf_test)


def test_label_shuffle_split():
    labels = [np.array([1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3]),
              np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]),
              np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]),
              np.array([1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4])]

    for l in labels:
        X = y = np.ones(len(l))
        n_iter = 6
        test_size = 1./3
        slo = LabelShuffleSplit(n_iter, test_size=test_size, random_state=0)

        # Make sure the repr works
        repr(slo)

        # Test that the length is correct
        assert_equal(slo.get_n_splits(X, y, labels=l), n_iter)

        l_unique = np.unique(l)

        for train, test in slo.split(X, y, labels=l):
            # First test: no train label is in the test set and vice versa
            l_train_unique = np.unique(l[train])
            l_test_unique = np.unique(l[test])
            assert_false(np.any(np.in1d(l[train], l_test_unique)))
            assert_false(np.any(np.in1d(l[test], l_train_unique)))

            # Second test: train and test add up to all the data
            assert_equal(l[train].size + l[test].size, l.size)

            # Third test: train and test are disjoint
            assert_array_equal(np.intersect1d(train, test), [])

            # Fourth test:
            # unique train and test labels are correct, +- 1 for rounding error
            assert_true(abs(len(l_test_unique) -
                            round(test_size * len(l_unique))) <= 1)
            assert_true(abs(len(l_train_unique) -
                            round((1.0 - test_size) * len(l_unique))) <= 1)


def test_leave_label_out_changing_labels():
    # Check that LeaveOneLabelOut and LeavePLabelOut work normally if
    # the labels variable is changed before calling split
    labels = np.array([0, 1, 2, 1, 1, 2, 0, 0])
    X = np.ones(len(labels))
    labels_changing = np.array(labels, copy=True)
    lolo = LeaveOneLabelOut().split(X, labels=labels)
    lolo_changing = LeaveOneLabelOut().split(X, labels=labels)
    lplo = LeavePLabelOut(n_labels=2).split(X, labels=labels)
    lplo_changing = LeavePLabelOut(n_labels=2).split(X, labels=labels)
    labels_changing[:] = 0
    for llo, llo_changing in [(lolo, lolo_changing), (lplo, lplo_changing)]:
        for (train, test), (train_chan, test_chan) in zip(llo, llo_changing):
            assert_array_equal(train, train_chan)
            assert_array_equal(test, test_chan)

    # n_splits = no of 2 (p) label combinations of the unique labels = 3C2 = 3
    assert_equal(3, LeavePLabelOut(n_labels=2).get_n_splits(X, y, labels))
    # n_splits = no of unique labels (C(uniq_lbls, 1) = n_unique_labels)
    assert_equal(3, LeaveOneLabelOut().get_n_splits(X, y, labels))


def test_train_test_split_errors():
    assert_raises(ValueError, train_test_split)
    assert_raises(ValueError, train_test_split, range(3), train_size=1.1)
    assert_raises(ValueError, train_test_split, range(3), test_size=0.6,
                  train_size=0.6)
    assert_raises(ValueError, train_test_split, range(3),
                  test_size=np.float32(0.6), train_size=np.float32(0.6))
    assert_raises(ValueError, train_test_split, range(3),
                  test_size="wrong_type")
    assert_raises(ValueError, train_test_split, range(3), test_size=2,
                  train_size=4)
    assert_raises(TypeError, train_test_split, range(3),
                  some_argument=1.1)
    assert_raises(ValueError, train_test_split, range(3), range(42))


def test_train_test_split():
    X = np.arange(100).reshape((10, 10))
    X_s = coo_matrix(X)
    y = np.arange(10)

    # simple test
    split = train_test_split(X, y, test_size=None, train_size=.5)
    X_train, X_test, y_train, y_test = split
    assert_equal(len(y_test), len(y_train))
    # test correspondence of X and y
    assert_array_equal(X_train[:, 0], y_train * 10)
    assert_array_equal(X_test[:, 0], y_test * 10)

    # don't convert lists to anything else by default
    split = train_test_split(X, X_s, y.tolist())
    X_train, X_test, X_s_train, X_s_test, y_train, y_test = split
    assert_true(isinstance(y_train, list))
    assert_true(isinstance(y_test, list))

    # allow nd-arrays
    X_4d = np.arange(10 * 5 * 3 * 2).reshape(10, 5, 3, 2)
    y_3d = np.arange(10 * 7 * 11).reshape(10, 7, 11)
    split = train_test_split(X_4d, y_3d)
    assert_equal(split[0].shape, (7, 5, 3, 2))
    assert_equal(split[1].shape, (3, 5, 3, 2))
    assert_equal(split[2].shape, (7, 7, 11))
    assert_equal(split[3].shape, (3, 7, 11))

    # test stratification option
    y = np.array([1, 1, 1, 1, 2, 2, 2, 2])
    for test_size, exp_test_size in zip([2, 4, 0.25, 0.5, 0.75],
                                        [2, 4, 2, 4, 6]):
        train, test = train_test_split(y, test_size=test_size,
                                       stratify=y,
                                       random_state=0)
        assert_equal(len(test), exp_test_size)
        assert_equal(len(test) + len(train), len(y))
        # check the 1:1 ratio of ones and twos in the data is preserved
        assert_equal(np.sum(train == 1), np.sum(train == 2))


@ignore_warnings
def train_test_split_pandas():
    # check train_test_split doesn't destroy pandas dataframe
    types = [MockDataFrame]
    try:
        from pandas import DataFrame
        types.append(DataFrame)
    except ImportError:
        pass
    for InputFeatureType in types:
        # X dataframe
        X_df = InputFeatureType(X)
        X_train, X_test = train_test_split(X_df)
        assert_true(isinstance(X_train, InputFeatureType))
        assert_true(isinstance(X_test, InputFeatureType))


def train_test_split_mock_pandas():
    # X mock dataframe
    X_df = MockDataFrame(X)
    X_train, X_test = train_test_split(X_df)
    assert_true(isinstance(X_train, MockDataFrame))
    assert_true(isinstance(X_test, MockDataFrame))
    X_train_arr, X_test_arr = train_test_split(X_df)


def test_shufflesplit_errors():
    # When the {test|train}_size is a float/invalid, error is raised at init
    assert_raises(ValueError, ShuffleSplit, test_size=None, train_size=None)
    assert_raises(ValueError, ShuffleSplit, test_size=2.0)
    assert_raises(ValueError, ShuffleSplit, test_size=1.0)
    assert_raises(ValueError, ShuffleSplit, test_size=0.1, train_size=0.95)
    assert_raises(ValueError, ShuffleSplit, train_size=1j)

    # When the {test|train}_size is an int, validation is based on the input X
    # and happens at split(...)
    assert_raises(ValueError, next, ShuffleSplit(test_size=11).split(X))
    assert_raises(ValueError, next, ShuffleSplit(test_size=10).split(X))
    assert_raises(ValueError, next, ShuffleSplit(test_size=8,
                                                 train_size=3).split(X))


def test_shufflesplit_reproducible():
    # Check that iterating twice on the ShuffleSplit gives the same
    # sequence of train-test when the random_state is given
    ss = ShuffleSplit(random_state=21)
    assert_array_equal(list(a for a, b in ss.split(X)),
                       list(a for a, b in ss.split(X)))


def test_safe_split_with_precomputed_kernel():
    clf = SVC()
    clfp = SVC(kernel="precomputed")

    X, y = iris.data, iris.target
    K = np.dot(X, X.T)

    cv = ShuffleSplit(test_size=0.25, random_state=0)
    tr, te = list(cv.split(X))[0]

    X_tr, y_tr = _safe_split(clf, X, y, tr)
    K_tr, y_tr2 = _safe_split(clfp, K, y, tr)
    assert_array_almost_equal(K_tr, np.dot(X_tr, X_tr.T))

    X_te, y_te = _safe_split(clf, X, y, te, tr)
    K_te, y_te2 = _safe_split(clfp, K, y, te, tr)
    assert_array_almost_equal(K_te, np.dot(X_te, X_tr.T))


def test_train_test_split_allow_nans():
    # Check that train_test_split allows input data with NaNs
    X = np.arange(200, dtype=np.float64).reshape(10, -1)
    X[2, :] = np.nan
    y = np.repeat([0, 1], X.shape[0] / 2)
    train_test_split(X, y, test_size=0.2, random_state=42)


def test_check_cv():
    X = np.ones(9)
    cv = check_cv(3, classifier=False)
    # Use numpy.testing.assert_equal which recursively compares
    # lists of lists
    np.testing.assert_equal(list(KFold(3).split(X)), list(cv.split(X)))

    y_binary = np.array([0, 1, 0, 1, 0, 0, 1, 1, 1])
    cv = check_cv(3, y_binary, classifier=True)
    np.testing.assert_equal(list(StratifiedKFold(3).split(X, y_binary)),
                            list(cv.split(X, y_binary)))

    y_multiclass = np.array([0, 1, 0, 1, 2, 1, 2, 0, 2])
    cv = check_cv(3, y_multiclass, classifier=True)
    np.testing.assert_equal(list(StratifiedKFold(3).split(X, y_multiclass)),
                            list(cv.split(X, y_multiclass)))

    X = np.ones(5)
    y_multilabel = np.array([[0, 0, 0, 0], [0, 1, 1, 0], [0, 0, 0, 1],
                            [1, 1, 0, 1], [0, 0, 1, 0]])
    cv = check_cv(3, y_multilabel, classifier=True)
    np.testing.assert_equal(list(KFold(3).split(X)), list(cv.split(X)))

    y_multioutput = np.array([[1, 2], [0, 3], [0, 0], [3, 1], [2, 0]])
    cv = check_cv(3, y_multioutput, classifier=True)
    np.testing.assert_equal(list(KFold(3).split(X)), list(cv.split(X)))

    # Check if the old style classes are wrapped to have a split method
    X = np.ones(9)
    y_multiclass = np.array([0, 1, 0, 1, 2, 1, 2, 0, 2])
    cv1 = check_cv(3, y_multiclass, classifier=True)

    with warnings.catch_warnings(record=True):
        from sklearn.cross_validation import StratifiedKFold as OldSKF

    cv2 = check_cv(OldSKF(y_multiclass, n_folds=3))
    np.testing.assert_equal(list(cv1.split(X, y_multiclass)),
                            list(cv2.split()))

    assert_raises(ValueError, check_cv, cv="lolo")


def test_cv_iterable_wrapper():
    y_multiclass = np.array([0, 1, 0, 1, 2, 1, 2, 0, 2])

    with warnings.catch_warnings(record=True):
        from sklearn.cross_validation import StratifiedKFold as OldSKF

    cv = OldSKF(y_multiclass, n_folds=3)
    wrapped_old_skf = _CVIterableWrapper(cv)

    # Check if split works correctly
    np.testing.assert_equal(list(cv), list(wrapped_old_skf.split()))

    # Check if get_n_splits works correctly
    assert_equal(len(cv), wrapped_old_skf.get_n_splits())


def test_label_kfold():
    rng = np.random.RandomState(0)

    # Parameters of the test
    n_labels = 15
    n_samples = 1000
    n_folds = 5

    X = y = np.ones(n_samples)

    # Construct the test data
    tolerance = 0.05 * n_samples  # 5 percent error allowed
    labels = rng.randint(0, n_labels, n_samples)

    ideal_n_labels_per_fold = n_samples // n_folds

    len(np.unique(labels))
    # Get the test fold indices from the test set indices of each fold
    folds = np.zeros(n_samples)
    lkf = LabelKFold(n_folds=n_folds)
    for i, (_, test) in enumerate(lkf.split(X, y, labels)):
        folds[test] = i

    # Check that folds have approximately the same size
    assert_equal(len(folds), len(labels))
    for i in np.unique(folds):
        assert_greater_equal(tolerance,
                             abs(sum(folds == i) - ideal_n_labels_per_fold))

    # Check that each label appears only in 1 fold
    for label in np.unique(labels):
        assert_equal(len(np.unique(folds[labels == label])), 1)

    # Check that no label is on both sides of the split
    labels = np.asarray(labels, dtype=object)
    for train, test in lkf.split(X, y, labels):
        assert_equal(len(np.intersect1d(labels[train], labels[test])), 0)

    # Construct the test data
    labels = np.array(['Albert', 'Jean', 'Bertrand', 'Michel', 'Jean',
                       'Francis', 'Robert', 'Michel', 'Rachel', 'Lois',
                       'Michelle', 'Bernard', 'Marion', 'Laura', 'Jean',
                       'Rachel', 'Franck', 'John', 'Gael', 'Anna', 'Alix',
                       'Robert', 'Marion', 'David', 'Tony', 'Abel', 'Becky',
                       'Madmood', 'Cary', 'Mary', 'Alexandre', 'David',
                       'Francis', 'Barack', 'Abdoul', 'Rasha', 'Xi', 'Silvia'])

    n_labels = len(np.unique(labels))
    n_samples = len(labels)
    n_folds = 5
    tolerance = 0.05 * n_samples  # 5 percent error allowed
    ideal_n_labels_per_fold = n_samples // n_folds

    X = y = np.ones(n_samples)

    # Get the test fold indices from the test set indices of each fold
    folds = np.zeros(n_samples)
    for i, (_, test) in enumerate(lkf.split(X, y, labels)):
        folds[test] = i

    # Check that folds have approximately the same size
    assert_equal(len(folds), len(labels))
    for i in np.unique(folds):
        assert_greater_equal(tolerance,
                             abs(sum(folds == i) - ideal_n_labels_per_fold))

    # Check that each label appears only in 1 fold
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        for label in np.unique(labels):
            assert_equal(len(np.unique(folds[labels == label])), 1)

    # Check that no label is on both sides of the split
    labels = np.asarray(labels, dtype=object)
    for train, test in lkf.split(X, y, labels):
        assert_equal(len(np.intersect1d(labels[train], labels[test])), 0)

    # Should fail if there are more folds than labels
    labels = np.array([1, 1, 1, 2, 2])
    X = y = np.ones(len(labels))
    assert_raises_regexp(ValueError, "Cannot have number of folds.*greater",
                         next, LabelKFold(n_folds=3).split(X, y, labels))


def test_nested_cv():
    # Test if nested cross validation works with different combinations of cv
    rng = np.random.RandomState(0)

    X, y = make_classification(n_samples=15, n_classes=2, random_state=0)
    labels = rng.randint(0, 5, 15)

    cvs = [LeaveOneLabelOut(), LeaveOneOut(), LabelKFold(), StratifiedKFold(),
           StratifiedShuffleSplit(n_iter=3, random_state=0)]

    for inner_cv, outer_cv in combinations_with_replacement(cvs, 2):
        gs = GridSearchCV(Ridge(), param_grid={'alpha': [1, .1]},
                          cv=inner_cv)
        cross_val_score(gs, X=X, y=y, labels=labels, cv=outer_cv,
                        fit_params={'labels': labels})


def test_build_repr():
    class MockSplitter:
        def __init__(self, a, b=0, c=None):
            self.a = a
            self.b = b
            self.c = c

        def __repr__(self):
            return _build_repr(self)

    assert_equal(repr(MockSplitter(5, 6)), "MockSplitter(a=5, b=6, c=None)")
