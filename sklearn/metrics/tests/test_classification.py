from __future__ import division, print_function

import numpy as np
from scipy import linalg
from functools import partial
from itertools import product
import warnings

from sklearn import datasets
from sklearn import svm

from sklearn.datasets import make_multilabel_classification
from sklearn.preprocessing import LabelBinarizer, MultiLabelBinarizer
from sklearn.utils.fixes import np_version
from sklearn.utils.validation import check_random_state

from sklearn.utils.testing import assert_raises, clean_warning_registry
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_no_warnings
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_not_equal
from sklearn.utils.testing import ignore_warnings

from sklearn.metrics import accuracy_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import fbeta_score
from sklearn.metrics import hamming_loss
from sklearn.metrics import hinge_loss
from sklearn.metrics import jaccard_similarity_score
from sklearn.metrics import log_loss
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import zero_one_loss
from sklearn.metrics import brier_score_loss


from sklearn.metrics.classification import _check_targets
from sklearn.metrics.base import UndefinedMetricWarning


###############################################################################
# Utilities for testing

def make_prediction(dataset=None, binary=False):
    """Make some classification predictions on a toy dataset using a SVC

    If binary is True restrict to a binary classification problem instead of a
    multiclass classification problem
    """

    if dataset is None:
        # import some data to play with
        dataset = datasets.load_iris()

    X = dataset.data
    y = dataset.target

    if binary:
        # restrict to a binary classification task
        X, y = X[y < 2], y[y < 2]

    n_samples, n_features = X.shape
    p = np.arange(n_samples)

    rng = check_random_state(37)
    rng.shuffle(p)
    X, y = X[p], y[p]
    half = int(n_samples / 2)

    # add noisy features to make the problem harder and avoid perfect results
    rng = np.random.RandomState(0)
    X = np.c_[X, rng.randn(n_samples, 200 * n_features)]

    # run classifier, get class probabilities and label predictions
    clf = svm.SVC(kernel='linear', probability=True, random_state=0)
    probas_pred = clf.fit(X[:half], y[:half]).predict_proba(X[half:])

    if binary:
        # only interested in probabilities of the positive case
        # XXX: do we really want a special API for the binary case?
        probas_pred = probas_pred[:, 1]

    y_pred = clf.predict(X[half:])
    y_true = y[half:]
    return y_true, y_pred, probas_pred


###############################################################################
# Tests


def test_multilabel_accuracy_score_subset_accuracy():
    # Dense label indicator matrix format
    y1 = np.array([[0, 1, 1], [1, 0, 1]])
    y2 = np.array([[0, 0, 1], [1, 0, 1]])

    assert_equal(accuracy_score(y1, y2), 0.5)
    assert_equal(accuracy_score(y1, y1), 1)
    assert_equal(accuracy_score(y2, y2), 1)
    assert_equal(accuracy_score(y2, np.logical_not(y2)), 0)
    assert_equal(accuracy_score(y1, np.logical_not(y1)), 0)
    assert_equal(accuracy_score(y1, np.zeros(y1.shape)), 0)
    assert_equal(accuracy_score(y2, np.zeros(y1.shape)), 0)

    with ignore_warnings():  # sequence of sequences is deprecated
        # List of tuple of label
        y1 = [(1, 2,), (0, 2,)]
        y2 = [(2,), (0, 2,)]

        assert_equal(accuracy_score(y1, y2), 0.5)
        assert_equal(accuracy_score(y1, y1), 1)
        assert_equal(accuracy_score(y2, y2), 1)
        assert_equal(accuracy_score(y2, [(), ()]), 0)
        assert_equal(accuracy_score(y1, y2, normalize=False), 1)
        assert_equal(accuracy_score(y1, y1, normalize=False), 2)
        assert_equal(accuracy_score(y2, y2, normalize=False), 2)
        assert_equal(accuracy_score(y2, [(), ()], normalize=False), 0)


def test_precision_recall_f1_score_binary():
    # Test Precision Recall and F1 Score for binary classification task
    y_true, y_pred, _ = make_prediction(binary=True)

    # detailed measures for each class
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred, average=None)
    assert_array_almost_equal(p, [0.73, 0.85], 2)
    assert_array_almost_equal(r, [0.88, 0.68], 2)
    assert_array_almost_equal(f, [0.80, 0.76], 2)
    assert_array_equal(s, [25, 25])

    # individual scoring function that can be used for grid search: in the
    # binary class case the score is the value of the measure for the positive
    # class (e.g. label == 1). This is deprecated for average != 'binary'.
    assert_dep_warning = partial(assert_warns, DeprecationWarning)
    for kwargs, my_assert in [({}, assert_no_warnings),
                              ({'average': 'binary'}, assert_no_warnings),
                              ({'average': 'micro'}, assert_dep_warning)]:
        ps = my_assert(precision_score, y_true, y_pred, **kwargs)
        assert_array_almost_equal(ps, 0.85, 2)

        rs = my_assert(recall_score, y_true, y_pred, **kwargs)
        assert_array_almost_equal(rs, 0.68, 2)

        fs = my_assert(f1_score, y_true, y_pred, **kwargs)
        assert_array_almost_equal(fs, 0.76, 2)

        assert_almost_equal(my_assert(fbeta_score, y_true, y_pred, beta=2,
                                      **kwargs),
                            (1 + 2 ** 2) * ps * rs / (2 ** 2 * ps + rs), 2)


@ignore_warnings
def test_precision_recall_f_binary_single_class():
    # Test precision, recall and F1 score behave with a single positive or
    # negative class
    # Such a case may occur with non-stratified cross-validation
    assert_equal(1., precision_score([1, 1], [1, 1]))
    assert_equal(1., recall_score([1, 1], [1, 1]))
    assert_equal(1., f1_score([1, 1], [1, 1]))

    assert_equal(0., precision_score([-1, -1], [-1, -1]))
    assert_equal(0., recall_score([-1, -1], [-1, -1]))
    assert_equal(0., f1_score([-1, -1], [-1, -1]))


def test_average_precision_score_score_non_binary_class():
    # Test that average_precision_score function returns an error when trying
    # to compute average_precision_score for multiclass task.
    rng = check_random_state(404)
    y_pred = rng.rand(10)

    # y_true contains three different class values
    y_true = rng.randint(0, 3, size=10)
    assert_raise_message(ValueError, "multiclass format is not supported",
                         average_precision_score, y_true, y_pred)


def test_average_precision_score_duplicate_values():
    # Duplicate values with precision-recall require a different
    # processing than when computing the AUC of a ROC, because the
    # precision-recall curve is a decreasing curve
    # The following situtation corresponds to a perfect
    # test statistic, the average_precision_score should be 1
    y_true = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]
    y_score = [0, .1, .1, .4, .5, .6, .6, .9, .9, 1, 1]
    assert_equal(average_precision_score(y_true, y_score), 1)


def test_average_precision_score_tied_values():
    # Here if we go from left to right in y_true, the 0 values are
    # are separated from the 1 values, so it appears that we've
    # Correctly sorted our classifications. But in fact the first two
    # values have the same score (0.5) and so the first two values
    # could be swapped around, creating an imperfect sorting. This
    # imperfection should come through in the end score, making it less
    # than one.
    y_true = [0, 1, 1]
    y_score = [.5, .5, .6]
    assert_not_equal(average_precision_score(y_true, y_score), 1.)


@ignore_warnings
def test_precision_recall_fscore_support_errors():
    y_true, y_pred, _ = make_prediction(binary=True)

    # Bad beta
    assert_raises(ValueError, precision_recall_fscore_support,
                  y_true, y_pred, beta=0.0)

    # Bad pos_label
    assert_raises(ValueError, precision_recall_fscore_support,
                  y_true, y_pred, pos_label=2, average='macro')

    # Bad average option
    assert_raises(ValueError, precision_recall_fscore_support,
                  [0, 1, 2], [1, 2, 0], average='mega')


def test_confusion_matrix_binary():
    # Test confusion matrix - binary classification case
    y_true, y_pred, _ = make_prediction(binary=True)

    def test(y_true, y_pred):
        cm = confusion_matrix(y_true, y_pred)
        assert_array_equal(cm, [[22, 3], [8, 17]])

        tp, fp, fn, tn = cm.flatten()
        num = (tp * tn - fp * fn)
        den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

        true_mcc = 0 if den == 0 else num / den
        mcc = matthews_corrcoef(y_true, y_pred)
        assert_array_almost_equal(mcc, true_mcc, decimal=2)
        assert_array_almost_equal(mcc, 0.57, decimal=2)

    test(y_true, y_pred)
    test([str(y) for y in y_true],
         [str(y) for y in y_pred])


@ignore_warnings
def test_matthews_corrcoef_nan():
    assert_equal(matthews_corrcoef([0], [1]), 0.0)
    assert_equal(matthews_corrcoef([0, 0], [0, 1]), 0.0)


def test_precision_recall_f1_score_multiclass():
    # Test Precision Recall and F1 Score for multiclass classification task
    y_true, y_pred, _ = make_prediction(binary=False)

    # compute scores with default labels introspection
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred, average=None)
    assert_array_almost_equal(p, [0.83, 0.33, 0.42], 2)
    assert_array_almost_equal(r, [0.79, 0.09, 0.90], 2)
    assert_array_almost_equal(f, [0.81, 0.15, 0.57], 2)
    assert_array_equal(s, [24, 31, 20])

    # averaging tests
    ps = precision_score(y_true, y_pred, pos_label=1, average='micro')
    assert_array_almost_equal(ps, 0.53, 2)

    rs = recall_score(y_true, y_pred, average='micro')
    assert_array_almost_equal(rs, 0.53, 2)

    fs = f1_score(y_true, y_pred, average='micro')
    assert_array_almost_equal(fs, 0.53, 2)

    ps = precision_score(y_true, y_pred, average='macro')
    assert_array_almost_equal(ps, 0.53, 2)

    rs = recall_score(y_true, y_pred, average='macro')
    assert_array_almost_equal(rs, 0.60, 2)

    fs = f1_score(y_true, y_pred, average='macro')
    assert_array_almost_equal(fs, 0.51, 2)

    ps = precision_score(y_true, y_pred, average='weighted')
    assert_array_almost_equal(ps, 0.51, 2)

    rs = recall_score(y_true, y_pred, average='weighted')
    assert_array_almost_equal(rs, 0.53, 2)

    fs = f1_score(y_true, y_pred, average='weighted')
    assert_array_almost_equal(fs, 0.47, 2)

    assert_raises(ValueError, precision_score, y_true, y_pred,
                  average="samples")
    assert_raises(ValueError, recall_score, y_true, y_pred, average="samples")
    assert_raises(ValueError, f1_score, y_true, y_pred, average="samples")
    assert_raises(ValueError, fbeta_score, y_true, y_pred, average="samples",
                  beta=0.5)

    # same prediction but with and explicit label ordering
    p, r, f, s = precision_recall_fscore_support(
        y_true, y_pred, labels=[0, 2, 1], average=None)
    assert_array_almost_equal(p, [0.83, 0.41, 0.33], 2)
    assert_array_almost_equal(r, [0.79, 0.90, 0.10], 2)
    assert_array_almost_equal(f, [0.81, 0.57, 0.15], 2)
    assert_array_equal(s, [24, 20, 31])


def test_precision_refcall_f1_score_multilabel_unordered_labels():
    # test that labels need not be sorted in the multilabel case
    y_true = np.array([[1, 1, 0, 0]])
    y_pred = np.array([[0, 0, 1, 1]])
    for average in ['samples', 'micro', 'macro', 'weighted', None]:
        p, r, f, s = precision_recall_fscore_support(
            y_true, y_pred, labels=[4, 1, 2, 3], warn_for=[], average=average)
        assert_array_equal(p, 0)
        assert_array_equal(r, 0)
        assert_array_equal(f, 0)
        if average is None:
            assert_array_equal(s, [0, 1, 1, 0])


def test_precision_recall_f1_score_multiclass_pos_label_none():
    # Test Precision Recall and F1 Score for multiclass classification task
    # GH Issue #1296
    # initialize data
    y_true = np.array([0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1])
    y_pred = np.array([1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1])

    # compute scores with default labels introspection
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                 pos_label=None,
                                                 average='weighted')


def test_zero_precision_recall():
    # Check that pathological cases do not bring NaNs

    old_error_settings = np.seterr(all='raise')

    try:
        y_true = np.array([0, 1, 2, 0, 1, 2])
        y_pred = np.array([2, 0, 1, 1, 2, 0])

        assert_almost_equal(precision_score(y_true, y_pred,
                                            average='weighted'), 0.0, 2)
        assert_almost_equal(recall_score(y_true, y_pred, average='weighted'),
                            0.0, 2)
        assert_almost_equal(f1_score(y_true, y_pred, average='weighted'),
                            0.0, 2)

    finally:
        np.seterr(**old_error_settings)


def test_confusion_matrix_multiclass():
    # Test confusion matrix - multi-class case
    y_true, y_pred, _ = make_prediction(binary=False)

    def test(y_true, y_pred, string_type=False):
        # compute confusion matrix with default labels introspection
        cm = confusion_matrix(y_true, y_pred)
        assert_array_equal(cm, [[19, 4, 1],
                                [4, 3, 24],
                                [0, 2, 18]])

        # compute confusion matrix with explicit label ordering
        labels = ['0', '2', '1'] if string_type else [0, 2, 1]
        cm = confusion_matrix(y_true,
                              y_pred,
                              labels=labels)
        assert_array_equal(cm, [[19, 1, 4],
                                [0, 18, 2],
                                [4, 24, 3]])

    test(y_true, y_pred)
    test(list(str(y) for y in y_true),
         list(str(y) for y in y_pred),
         string_type=True)


def test_confusion_matrix_multiclass_subset_labels():
    # Test confusion matrix - multi-class case with subset of labels
    y_true, y_pred, _ = make_prediction(binary=False)

    # compute confusion matrix with only first two labels considered
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    assert_array_equal(cm, [[19, 4],
                            [4, 3]])

    # compute confusion matrix with explicit label ordering for only subset
    # of labels
    cm = confusion_matrix(y_true, y_pred, labels=[2, 1])
    assert_array_equal(cm, [[18, 2],
                            [24, 3]])


def test_classification_report_multiclass():
    # Test performance report
    iris = datasets.load_iris()
    y_true, y_pred, _ = make_prediction(dataset=iris, binary=False)

    # print classification report with class names
    expected_report = """\
             precision    recall  f1-score   support

     setosa       0.83      0.79      0.81        24
 versicolor       0.33      0.10      0.15        31
  virginica       0.42      0.90      0.57        20

avg / total       0.51      0.53      0.47        75
"""
    report = classification_report(
        y_true, y_pred, labels=np.arange(len(iris.target_names)),
        target_names=iris.target_names)
    assert_equal(report, expected_report)
    # print classification report with label detection
    expected_report = """\
             precision    recall  f1-score   support

          0       0.83      0.79      0.81        24
          1       0.33      0.10      0.15        31
          2       0.42      0.90      0.57        20

avg / total       0.51      0.53      0.47        75
"""
    report = classification_report(y_true, y_pred)
    assert_equal(report, expected_report)


def test_classification_report_multiclass_with_digits():
    # Test performance report with added digits in floating point values
    iris = datasets.load_iris()
    y_true, y_pred, _ = make_prediction(dataset=iris, binary=False)

    # print classification report with class names
    expected_report = """\
             precision    recall  f1-score   support

     setosa    0.82609   0.79167   0.80851        24
 versicolor    0.33333   0.09677   0.15000        31
  virginica    0.41860   0.90000   0.57143        20

avg / total    0.51375   0.53333   0.47310        75
"""
    report = classification_report(
        y_true, y_pred, labels=np.arange(len(iris.target_names)),
        target_names=iris.target_names, digits=5)
    assert_equal(report, expected_report)
    # print classification report with label detection
    expected_report = """\
             precision    recall  f1-score   support

          0       0.83      0.79      0.81        24
          1       0.33      0.10      0.15        31
          2       0.42      0.90      0.57        20

avg / total       0.51      0.53      0.47        75
"""
    report = classification_report(y_true, y_pred)
    assert_equal(report, expected_report)


def test_classification_report_multiclass_with_string_label():
    y_true, y_pred, _ = make_prediction(binary=False)

    y_true = np.array(["blue", "green", "red"])[y_true]
    y_pred = np.array(["blue", "green", "red"])[y_pred]

    expected_report = """\
             precision    recall  f1-score   support

       blue       0.83      0.79      0.81        24
      green       0.33      0.10      0.15        31
        red       0.42      0.90      0.57        20

avg / total       0.51      0.53      0.47        75
"""
    report = classification_report(y_true, y_pred)
    assert_equal(report, expected_report)

    expected_report = """\
             precision    recall  f1-score   support

          a       0.83      0.79      0.81        24
          b       0.33      0.10      0.15        31
          c       0.42      0.90      0.57        20

avg / total       0.51      0.53      0.47        75
"""
    report = classification_report(y_true, y_pred,
                                   target_names=["a", "b", "c"])
    assert_equal(report, expected_report)


def test_classification_report_multiclass_with_unicode_label():
    y_true, y_pred, _ = make_prediction(binary=False)

    labels = np.array([u"blue\xa2", u"green\xa2", u"red\xa2"])
    y_true = labels[y_true]
    y_pred = labels[y_pred]

    expected_report = u"""\
             precision    recall  f1-score   support

      blue\xa2       0.83      0.79      0.81        24
     green\xa2       0.33      0.10      0.15        31
       red\xa2       0.42      0.90      0.57        20

avg / total       0.51      0.53      0.47        75
"""
    if np_version[:3] < (1, 7, 0):
        expected_message = ("NumPy < 1.7.0 does not implement"
                            " searchsorted on unicode data correctly.")
        assert_raise_message(RuntimeError, expected_message,
                             classification_report, y_true, y_pred)
    else:
        report = classification_report(y_true, y_pred)
        assert_equal(report, expected_report)
@ignore_warnings  # sequence of sequences is deprecated
def test_multilabel_classification_report():
    n_classes = 4
    n_samples = 50
    make_ml = make_multilabel_classification
    _, y_true_ll = make_ml(n_features=1, n_classes=n_classes, random_state=0,
                           n_samples=n_samples)
    _, y_pred_ll = make_ml(n_features=1, n_classes=n_classes, random_state=1,
                           n_samples=n_samples)
    expected_report = """\
             precision    recall  f1-score   support

          0       0.50      0.67      0.57        24
          1       0.51      0.74      0.61        27
          2       0.29      0.08      0.12        26
          3       0.52      0.56      0.54        27

avg / total       0.45      0.51      0.46       104
"""

    lb = MultiLabelBinarizer()
    lb.fit([range(4)])
    y_true_bi = lb.transform(y_true_ll)
    y_pred_bi = lb.transform(y_pred_ll)

    for y_true, y_pred in [(y_true_ll, y_pred_ll), (y_true_bi, y_pred_bi)]:
        report = classification_report(y_true, y_pred)
        assert_equal(report, expected_report)


def test_multilabel_zero_one_loss_subset():
    # Dense label indicator matrix format
    y1 = np.array([[0, 1, 1], [1, 0, 1]])
    y2 = np.array([[0, 0, 1], [1, 0, 1]])

    assert_equal(zero_one_loss(y1, y2), 0.5)
    assert_equal(zero_one_loss(y1, y1), 0)
    assert_equal(zero_one_loss(y2, y2), 0)
    assert_equal(zero_one_loss(y2, np.logical_not(y2)), 1)
    assert_equal(zero_one_loss(y1, np.logical_not(y1)), 1)
    assert_equal(zero_one_loss(y1, np.zeros(y1.shape)), 1)
    assert_equal(zero_one_loss(y2, np.zeros(y1.shape)), 1)

    with ignore_warnings():  # sequence of sequences is deprecated
        # List of tuple of label
        y1 = [(1, 2,), (0, 2,)]
        y2 = [(2,), (0, 2,)]

        assert_equal(zero_one_loss(y1, y2), 0.5)
        assert_equal(zero_one_loss(y1, y1), 0)
        assert_equal(zero_one_loss(y2, y2), 0)
        assert_equal(zero_one_loss(y2, [(), ()]), 1)
        assert_equal(zero_one_loss(y2, [tuple(), (10, )]), 1)


def test_multilabel_hamming_loss():
    # Dense label indicator matrix format
    y1 = np.array([[0, 1, 1], [1, 0, 1]])
    y2 = np.array([[0, 0, 1], [1, 0, 1]])

    assert_equal(hamming_loss(y1, y2), 1 / 6)
    assert_equal(hamming_loss(y1, y1), 0)
    assert_equal(hamming_loss(y2, y2), 0)
    assert_equal(hamming_loss(y2, np.logical_not(y2)), 1)
    assert_equal(hamming_loss(y1, np.logical_not(y1)), 1)
    assert_equal(hamming_loss(y1, np.zeros(y1.shape)), 4 / 6)
    assert_equal(hamming_loss(y2, np.zeros(y1.shape)), 0.5)

    with ignore_warnings():  # sequence of sequences is deprecated
        # List of tuple of label
        y1 = [(1, 2,), (0, 2,)]
        y2 = [(2,), (0, 2,)]

        assert_equal(hamming_loss(y1, y2), 1 / 6)
        assert_equal(hamming_loss(y1, y1), 0)
        assert_equal(hamming_loss(y2, y2), 0)
        assert_equal(hamming_loss(y2, [(), ()]), 0.75)
        assert_equal(hamming_loss(y1, [tuple(), (10, )]), 0.625)
        assert_almost_equal(hamming_loss(y2, [tuple(), (10, )],
                                         classes=np.arange(11)), 0.1818, 2)


def test_multilabel_jaccard_similarity_score():
    # Dense label indicator matrix format
    y1 = np.array([[0, 1, 1], [1, 0, 1]])
    y2 = np.array([[0, 0, 1], [1, 0, 1]])

    # size(y1 \inter y2) = [1, 2]
    # size(y1 \union y2) = [2, 2]

    assert_equal(jaccard_similarity_score(y1, y2), 0.75)
    assert_equal(jaccard_similarity_score(y1, y1), 1)
    assert_equal(jaccard_similarity_score(y2, y2), 1)
    assert_equal(jaccard_similarity_score(y2, np.logical_not(y2)), 0)
    assert_equal(jaccard_similarity_score(y1, np.logical_not(y1)), 0)
    assert_equal(jaccard_similarity_score(y1, np.zeros(y1.shape)), 0)
    assert_equal(jaccard_similarity_score(y2, np.zeros(y1.shape)), 0)

    with ignore_warnings():  # sequence of sequences is deprecated
        # List of tuple of label
        y1 = [(1, 2,), (0, 2,)]
        y2 = [(2,), (0, 2,)]

        assert_equal(jaccard_similarity_score(y1, y2), 0.75)
        assert_equal(jaccard_similarity_score(y1, y1), 1)
        assert_equal(jaccard_similarity_score(y2, y2), 1)
        assert_equal(jaccard_similarity_score(y2, [(), ()]), 0)

        # |y3 inter y4 | = [0, 1, 1]
        # |y3 union y4 | = [2, 1, 3]
        y3 = [(0,), (1,), (3,)]
        y4 = [(4,), (4,), (5, 6)]
        assert_almost_equal(jaccard_similarity_score(y3, y4), 0)

        # |y5 inter y6 | = [0, 1, 1]
        # |y5 union y6 | = [2, 1, 3]
        y5 = [(0,), (1,), (2, 3)]
        y6 = [(1,), (1,), (2, 0)]

        assert_almost_equal(jaccard_similarity_score(y5, y6), (1 + 1 / 3) / 3)


@ignore_warnings
def test_precision_recall_f1_score_multilabel_1():
    # Test precision_recall_f1_score on a crafted multilabel example
    # First crafted example
    y_true_ll = [(0,), (1,), (2, 3)]
    y_pred_ll = [(1,), (1,), (2, 0)]
    lb = LabelBinarizer()
    lb.fit([range(4)])
    y_true_bi = lb.transform(y_true_ll)
    y_pred_bi = lb.transform(y_pred_ll)

    for y_true, y_pred in [(y_true_ll, y_pred_ll), (y_true_bi, y_pred_bi)]:

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average=None)
        #tp = [0, 1, 1, 0]
        #fn = [1, 0, 0, 1]
        #fp = [1, 1, 0, 0]
        # Check per class

        assert_array_almost_equal(p, [0.0, 0.5, 1.0, 0.0], 2)
        assert_array_almost_equal(r, [0.0, 1.0, 1.0, 0.0], 2)
        assert_array_almost_equal(f, [0.0, 1 / 1.5, 1, 0.0], 2)
        assert_array_almost_equal(s, [1, 1, 1, 1], 2)

        f2 = fbeta_score(y_true, y_pred, beta=2, average=None)
        support = s
        assert_array_almost_equal(f2, [0, 0.83, 1, 0], 2)

        # Check macro
        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="macro")
        assert_almost_equal(p, 1.5 / 4)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 2.5 / 1.5 * 0.25)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="macro"),
                            np.mean(f2))

        # Check micro
        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="micro")
        assert_almost_equal(p, 0.5)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 0.5)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="micro"),
                            (1 + 4) * p * r / (4 * p + r))

        # Check weigted
        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="weighted")
        assert_almost_equal(p, 1.5 / 4)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 2.5 / 1.5 * 0.25)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="weighted"),
                            np.average(f2, weights=support))
        # Check weigted
        # |h(x_i) inter y_i | = [0, 1, 1]
        # |y_i| = [1, 1, 2]
        # |h(x_i)| = [1, 1, 2]
        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="samples")
        assert_almost_equal(p, 0.5)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 0.5)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="samples"),
                            0.5)


@ignore_warnings
def test_precision_recall_f1_score_multilabel_2():
    # Test precision_recall_f1_score on a crafted multilabel example 2
    # Second crafted example
    y_true_ll = [(1,), (2,), (2, 3)]
    y_pred_ll = [(4,), (4,), (2, 1)]
    lb = LabelBinarizer()
    lb.fit([range(1, 5)])
    y_true_bi = lb.transform(y_true_ll)
    y_pred_bi = lb.transform(y_pred_ll)

    for y_true, y_pred in [(y_true_ll, y_pred_ll), (y_true_bi, y_pred_bi)]:
        # tp = [ 0.  1.  0.  0.]
        # fp = [ 1.  0.  0.  2.]
        # fn = [ 1.  1.  1.  0.]

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average=None)
        assert_array_almost_equal(p, [0.0, 1.0, 0.0, 0.0], 2)
        assert_array_almost_equal(r, [0.0, 0.5, 0.0, 0.0], 2)
        assert_array_almost_equal(f, [0.0, 0.66, 0.0, 0.0], 2)
        assert_array_almost_equal(s, [1, 2, 1, 0], 2)

        f2 = fbeta_score(y_true, y_pred, beta=2, average=None)
        support = s
        assert_array_almost_equal(f2, [0, 0.55, 0, 0], 2)

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="micro")
        assert_almost_equal(p, 0.25)
        assert_almost_equal(r, 0.25)
        assert_almost_equal(f, 2 * 0.25 * 0.25 / 0.5)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="micro"),
                            (1 + 4) * p * r / (4 * p + r))

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="macro")
        assert_almost_equal(p, 0.25)
        assert_almost_equal(r, 0.125)
        assert_almost_equal(f, 2 / 12)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="macro"),
                            np.mean(f2))

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="weighted")
        assert_almost_equal(p, 2 / 4)
        assert_almost_equal(r, 1 / 4)
        assert_almost_equal(f, 2 / 3 * 2 / 4)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="weighted"),
                            np.average(f2, weights=support))

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="samples")
        # Check weigted
        # |h(x_i) inter y_i | = [0, 0, 1]
        # |y_i| = [1, 1, 2]
        # |h(x_i)| = [1, 1, 2]

        assert_almost_equal(p, 1 / 6)
        assert_almost_equal(r, 1 / 6)
        assert_almost_equal(f, 2 / 4 * 1 / 3)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="samples"),
                            0.1666, 2)


@ignore_warnings
def test_precision_recall_f1_score_with_an_empty_prediction():
    y_true_ll = [(1,), (0,), (2, 1,)]
    y_pred_ll = [tuple(), (3,), (2, 1)]

    lb = LabelBinarizer()
    lb.fit([range(4)])
    y_true_bi = lb.transform(y_true_ll)
    y_pred_bi = lb.transform(y_pred_ll)

    for y_true, y_pred in [(y_true_ll, y_pred_ll), (y_true_bi, y_pred_bi)]:
        # true_pos = [ 0.  1.  1.  0.]
        # false_pos = [ 0.  0.  0.  1.]
        # false_neg = [ 1.  1.  0.  0.]
        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average=None)
        assert_array_almost_equal(p, [0.0, 1.0, 1.0, 0.0], 2)
        assert_array_almost_equal(r, [0.0, 0.5, 1.0, 0.0], 2)
        assert_array_almost_equal(f, [0.0, 1 / 1.5, 1, 0.0], 2)
        assert_array_almost_equal(s, [1, 2, 1, 0], 2)

        f2 = fbeta_score(y_true, y_pred, beta=2, average=None)
        support = s
        assert_array_almost_equal(f2, [0, 0.55, 1, 0], 2)

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="macro")
        assert_almost_equal(p, 0.5)
        assert_almost_equal(r, 1.5 / 4)
        assert_almost_equal(f, 2.5 / (4 * 1.5))
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="macro"),
                            np.mean(f2))

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="micro")
        assert_almost_equal(p, 2 / 3)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 2 / 3 / (2 / 3 + 0.5))
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="micro"),
                            (1 + 4) * p * r / (4 * p + r))

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="weighted")
        assert_almost_equal(p, 3 / 4)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, (2 / 1.5 + 1) / 4)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="weighted"),
                            np.average(f2, weights=support))

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="samples")
        # |h(x_i) inter y_i | = [0, 0, 2]
        # |y_i| = [1, 1, 2]
        # |h(x_i)| = [0, 1, 2]
        assert_almost_equal(p, 1 / 3)
        assert_almost_equal(r, 1 / 3)
        assert_almost_equal(f, 1 / 3)
        assert_equal(s, None)
        assert_almost_equal(fbeta_score(y_true, y_pred, beta=2,
                                        average="samples"),
                            0.333, 2)


def test_precision_recall_f1_no_labels():
    y_true = np.zeros((20, 3))
    y_pred = np.zeros_like(y_true)

    # tp = [0, 0, 0]
    # fn = [0, 0, 0]
    # fp = [0, 0, 0]
    # support = [0, 0, 0]
    # |y_hat_i inter y_i | = [0, 0, 0]
    # |y_i| = [0, 0, 0]
    # |y_hat_i| = [0, 0, 0]

    for beta in [1]:
        p, r, f, s = assert_warns(UndefinedMetricWarning,
                                  precision_recall_fscore_support,
                                  y_true, y_pred, average=None, beta=beta)
        assert_array_almost_equal(p, [0, 0, 0], 2)
        assert_array_almost_equal(r, [0, 0, 0], 2)
        assert_array_almost_equal(f, [0, 0, 0], 2)
        assert_array_almost_equal(s, [0, 0, 0], 2)

        fbeta = assert_warns(UndefinedMetricWarning, fbeta_score,
                             y_true, y_pred, beta=beta, average=None)
        assert_array_almost_equal(fbeta, [0, 0, 0], 2)

        for average in ["macro", "micro", "weighted", "samples"]:
            p, r, f, s = assert_warns(UndefinedMetricWarning,
                                      precision_recall_fscore_support,
                                      y_true, y_pred, average=average,
                                      beta=beta)
            assert_almost_equal(p, 0)
            assert_almost_equal(r, 0)
            assert_almost_equal(f, 0)
            assert_equal(s, None)

            fbeta = assert_warns(UndefinedMetricWarning, fbeta_score,
                                 y_true, y_pred,
                                 beta=beta, average=average)
            assert_almost_equal(fbeta, 0)


def test_prf_warnings():

    # average of per-label scores
    f, w = precision_recall_fscore_support, UndefinedMetricWarning
    my_assert = assert_warns_message
    for average in [None, 'weighted', 'macro']:
        msg = ('Precision and F-score are ill-defined and '
               'being set to 0.0 in labels with no predicted samples.')
        my_assert(w, msg, f, [0, 1, 2], [1, 1, 2], average=average)

        msg = ('Recall and F-score are ill-defined and '
               'being set to 0.0 in labels with no true samples.')
        my_assert(w, msg, f, [1, 1, 2], [0, 1, 2], average=average)

        # average of per-sample scores
        msg = ('Precision and F-score are ill-defined and '
               'being set to 0.0 in samples with no predicted labels.')
        my_assert(w, msg, f, np.array([[1, 0], [1, 0]]),
                  np.array([[1, 0], [0, 0]]), average='samples')

        msg = ('Recall and F-score are ill-defined and '
               'being set to 0.0 in samples with no true labels.')
        my_assert(w, msg, f, np.array([[1, 0], [0, 0]]),
                  np.array([[1, 0], [1, 0]]),
                  average='samples')

        # single score: micro-average
        msg = ('Precision and F-score are ill-defined and '
               'being set to 0.0 due to no predicted samples.')
        my_assert(w, msg, f, np.array([[1, 1], [1, 1]]),
                  np.array([[0, 0], [0, 0]]), average='micro')

        msg = ('Recall and F-score are ill-defined and '
               'being set to 0.0 due to no true samples.')
        my_assert(w, msg, f, np.array([[0, 0], [0, 0]]),
                  np.array([[1, 1], [1, 1]]), average='micro')

        # single postive label
        msg = ('Precision and F-score are ill-defined and '
               'being set to 0.0 due to no predicted samples.')
        my_assert(w, msg, f, [1, 1], [-1, -1], average='macro')

        msg = ('Recall and F-score are ill-defined and '
               'being set to 0.0 due to no true samples.')
        my_assert(w, msg, f, [-1, -1], [1, 1], average='macro')


def test_recall_warnings():
    assert_no_warnings(recall_score,
                       np.array([[1, 1], [1, 1]]),
                       np.array([[0, 0], [0, 0]]),
                       average='micro')
    clean_warning_registry()
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter('always')
        recall_score(np.array([[0, 0], [0, 0]]),
                     np.array([[1, 1], [1, 1]]),
                     average='micro')
        assert_equal(str(record.pop().message),
                     'Recall is ill-defined and '
                     'being set to 0.0 due to no true samples.')


def test_precision_warnings():
    clean_warning_registry()
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter('always')

        precision_score(np.array([[1, 1], [1, 1]]),
                        np.array([[0, 0], [0, 0]]),
                        average='micro')
        assert_equal(str(record.pop().message),
                     'Precision is ill-defined and '
                     'being set to 0.0 due to no predicted samples.')

    assert_no_warnings(precision_score,
                       np.array([[0, 0], [0, 0]]),
                       np.array([[1, 1], [1, 1]]),
                       average='micro')


def test_fscore_warnings():
    clean_warning_registry()
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter('always')

        for score in [f1_score, partial(fbeta_score, beta=2)]:
            score(np.array([[1, 1], [1, 1]]),
                  np.array([[0, 0], [0, 0]]),
                  average='micro')
            assert_equal(str(record.pop().message),
                         'F-score is ill-defined and '
                         'being set to 0.0 due to no predicted samples.')
            score(np.array([[0, 0], [0, 0]]),
                  np.array([[1, 1], [1, 1]]),
                  average='micro')
            assert_equal(str(record.pop().message),
                         'F-score is ill-defined and '
                         'being set to 0.0 due to no true samples.')


def test_prf_average_compat():
    # Ensure warning if f1_score et al.'s average is implicit for multiclass
    y_true = [1, 2, 3, 3]
    y_pred = [1, 2, 3, 1]
    y_true_bin = [0, 1, 1]
    y_pred_bin = [0, 1, 0]

    for metric in [precision_score, recall_score, f1_score,
                   partial(fbeta_score, beta=2)]:
        score = assert_warns(DeprecationWarning, metric, y_true, y_pred)
        score_weighted = assert_no_warnings(metric, y_true, y_pred,
                                            average='weighted')
        assert_equal(score, score_weighted,
                     'average does not act like "weighted" by default')

        # check binary passes without warning
        assert_no_warnings(metric, y_true_bin, y_pred_bin)

        # but binary with pos_label=None should behave like multiclass
        score = assert_warns(DeprecationWarning, metric,
                             y_true_bin, y_pred_bin, pos_label=None)
        score_weighted = assert_no_warnings(metric, y_true_bin, y_pred_bin,
                                            pos_label=None, average='weighted')
        assert_equal(score, score_weighted,
                     'average does not act like "weighted" by default with '
                     'binary data and pos_label=None')


@ignore_warnings  # sequence of sequences is deprecated
def test__check_targets():
    # Check that _check_targets correctly merges target types, squeezes
    # output and fails if input lengths differ.
    IND = 'multilabel-indicator'
    SEQ = 'multilabel-sequences'
    MC = 'multiclass'
    BIN = 'binary'
    CNT = 'continuous'
    MMC = 'multiclass-multioutput'
    MCN = 'continuous-multioutput'
    # all of length 3
    EXAMPLES = [
        (IND, np.array([[0, 1, 1], [1, 0, 0], [0, 0, 1]])),
        # must not be considered binary
        (IND, np.array([[0, 1], [1, 0], [1, 1]])),
        (SEQ, [[2, 3], [1], [3]]),
        (MC, [2, 3, 1]),
        (BIN, [0, 1, 1]),
        (CNT, [0., 1.5, 1.]),
        (MC, np.array([[2], [3], [1]])),
        (BIN, np.array([[0], [1], [1]])),
        (CNT, np.array([[0.], [1.5], [1.]])),
        (MMC, np.array([[0, 2], [1, 3], [2, 3]])),
        (MCN, np.array([[0.5, 2.], [1.1, 3.], [2., 3.]])),
    ]
    # expected type given input types, or None for error
    # (types will be tried in either order)
    EXPECTED = {
        (IND, IND): IND,
        (SEQ, SEQ): IND,
        (MC, MC): MC,
        (BIN, BIN): BIN,

        (IND, SEQ): None,
        (MC, SEQ): None,
        (BIN, SEQ): None,
        (MC, IND): None,
        (BIN, IND): None,
        (BIN, MC): MC,

        # Disallowed types
        (CNT, CNT): None,
        (MMC, MMC): None,
        (MCN, MCN): None,
        (IND, CNT): None,
        (SEQ, CNT): None,
        (MC, CNT): None,
        (BIN, CNT): None,
        (MMC, CNT): None,
        (MCN, CNT): None,
        (IND, MMC): None,
        (SEQ, MMC): None,
        (MC, MMC): None,
        (BIN, MMC): None,
        (MCN, MMC): None,
        (IND, MCN): None,
        (SEQ, MCN): None,
        (MC, MCN): None,
        (BIN, MCN): None,
    }

    for (type1, y1), (type2, y2) in product(EXAMPLES, repeat=2):
        try:
            expected = EXPECTED[type1, type2]
        except KeyError:
            expected = EXPECTED[type2, type1]
        if expected is None:
            assert_raises(ValueError, _check_targets, y1, y2)

            if type1 != type2:
                assert_raise_message(
                    ValueError,
                    "Can't handle mix of {0} and {1}".format(type1, type2),
                    _check_targets, y1, y2)

            else:
                if type1 not in (BIN, MC, SEQ, IND):
                    assert_raise_message(ValueError,
                                         "{0} is not supported".format(type1),
                                         _check_targets, y1, y2)

        else:
            merged_type, y1out, y2out = _check_targets(y1, y2)
            assert_equal(merged_type, expected)
            if merged_type.startswith('multilabel'):
                assert_equal(y1out.format, 'csr')
                assert_equal(y2out.format, 'csr')
            else:
                assert_array_equal(y1out, np.squeeze(y1))
                assert_array_equal(y2out, np.squeeze(y2))
            assert_raises(ValueError, _check_targets, y1[:-1], y2)


def test_hinge_loss_binary():
    y_true = np.array([-1, 1, 1, -1])
    pred_decision = np.array([-8.5, 0.5, 1.5, -0.3])
    assert_equal(hinge_loss(y_true, pred_decision), 1.2 / 4)

    y_true = np.array([0, 2, 2, 0])
    pred_decision = np.array([-8.5, 0.5, 1.5, -0.3])
    assert_equal(hinge_loss(y_true, pred_decision), 1.2 / 4)


def test_hinge_loss_multiclass():
    pred_decision = np.array([
        [0.36, -0.17, -0.58, -0.99],
        [-0.54, -0.37, -0.48, -0.58],
        [-1.45, -0.58, -0.38, -0.17],
        [-0.54, -0.38, -0.48, -0.58],
        [-2.36, -0.79, -0.27,  0.24],
        [-1.45, -0.58, -0.38, -0.17]
    ])
    y_true = np.array([0, 1, 2, 1, 3, 2])
    dummy_losses = np.array([
        1 - pred_decision[0][0] + pred_decision[0][1],
        1 - pred_decision[1][1] + pred_decision[1][2],
        1 - pred_decision[2][2] + pred_decision[2][3],
        1 - pred_decision[3][1] + pred_decision[3][2],
        1 - pred_decision[4][3] + pred_decision[4][2],
        1 - pred_decision[5][2] + pred_decision[5][3]
    ])
    dummy_losses[dummy_losses <= 0] = 0
    dummy_hinge_loss = np.mean(dummy_losses)
    assert_equal(hinge_loss(y_true, pred_decision),
                 dummy_hinge_loss)


def test_hinge_loss_multiclass_missing_labels_with_labels_none():
    y_true = np.array([0, 1, 2, 2])
    pred_decision = np.array([
        [1.27, 0.034, -0.68, -1.40],
        [-1.45, -0.58, -0.38, -0.17],
        [-2.36, -0.79, -0.27,  0.24],
        [-2.36, -0.79, -0.27,  0.24]
    ])
    error_message = ("Please include all labels in y_true "
                     "or pass labels as third argument")
    assert_raise_message(ValueError,
                         error_message,
                         hinge_loss, y_true, pred_decision)


def test_hinge_loss_multiclass_with_missing_labels():
    pred_decision = np.array([
        [0.36, -0.17, -0.58, -0.99],
        [-0.55, -0.38, -0.48, -0.58],
        [-1.45, -0.58, -0.38, -0.17],
        [-0.55, -0.38, -0.48, -0.58],
        [-1.45, -0.58, -0.38, -0.17]
    ])
    y_true = np.array([0, 1, 2, 1, 2])
    labels = np.array([0, 1, 2, 3])
    dummy_losses = np.array([
        1 - pred_decision[0][0] + pred_decision[0][1],
        1 - pred_decision[1][1] + pred_decision[1][2],
        1 - pred_decision[2][2] + pred_decision[2][3],
        1 - pred_decision[3][1] + pred_decision[3][2],
        1 - pred_decision[4][2] + pred_decision[4][3]
    ])
    dummy_losses[dummy_losses <= 0] = 0
    dummy_hinge_loss = np.mean(dummy_losses)
    assert_equal(hinge_loss(y_true, pred_decision, labels=labels),
                 dummy_hinge_loss)


def test_hinge_loss_multiclass_invariance_lists():
    # Currently, invariance of string and integer labels cannot be tested
    # in common invariance tests because invariance tests for multiclass
    # decision functions is not implemented yet.
    y_true = ['blue', 'green', 'red',
              'green', 'white', 'red']
    pred_decision = [
        [0.36, -0.17, -0.58, -0.99],
        [-0.55, -0.38, -0.48, -0.58],
        [-1.45, -0.58, -0.38,  -0.17],
        [-0.55, -0.38, -0.48, -0.58],
        [-2.36, -0.79, -0.27,  0.24],
        [-1.45, -0.58, -0.38,  -0.17]]
    dummy_losses = np.array([
        1 - pred_decision[0][0] + pred_decision[0][1],
        1 - pred_decision[1][1] + pred_decision[1][2],
        1 - pred_decision[2][2] + pred_decision[2][3],
        1 - pred_decision[3][1] + pred_decision[3][2],
        1 - pred_decision[4][3] + pred_decision[4][2],
        1 - pred_decision[5][2] + pred_decision[5][3]
    ])
    dummy_losses[dummy_losses <= 0] = 0
    dummy_hinge_loss = np.mean(dummy_losses)
    assert_equal(hinge_loss(y_true, pred_decision),
                 dummy_hinge_loss)


def test_log_loss():
    # binary case with symbolic labels ("no" < "yes")
    y_true = ["no", "no", "no", "yes", "yes", "yes"]
    y_pred = np.array([[0.5, 0.5], [0.1, 0.9], [0.01, 0.99],
                       [0.9, 0.1], [0.75, 0.25], [0.001, 0.999]])
    loss = log_loss(y_true, y_pred)
    assert_almost_equal(loss, 1.8817971)

    # multiclass case; adapted from http://bit.ly/RJJHWA
    y_true = [1, 0, 2]
    y_pred = [[0.2, 0.7, 0.1], [0.6, 0.2, 0.2], [0.6, 0.1, 0.3]]
    loss = log_loss(y_true, y_pred, normalize=True)
    assert_almost_equal(loss, 0.6904911)

    # check that we got all the shapes and axes right
    # by doubling the length of y_true and y_pred
    y_true *= 2
    y_pred *= 2
    loss = log_loss(y_true, y_pred, normalize=False)
    assert_almost_equal(loss, 0.6904911 * 6, decimal=6)

    # check eps and handling of absolute zero and one probabilities
    y_pred = np.asarray(y_pred) > .5
    loss = log_loss(y_true, y_pred, normalize=True, eps=.1)
    assert_almost_equal(loss, log_loss(y_true, np.clip(y_pred, .1, .9)))

    # raise error if number of classes are not equal.
    y_true = [1, 0, 2]
    y_pred = [[0.2, 0.7], [0.6, 0.5], [0.4, 0.1]]
    assert_raises(ValueError, log_loss, y_true, y_pred)

    # case when y_true is a string array object
    y_true = ["ham", "spam", "spam", "ham"]
    y_pred = [[0.2, 0.7], [0.6, 0.5], [0.4, 0.1], [0.7, 0.2]]
    loss = log_loss(y_true, y_pred)
    assert_almost_equal(loss, 1.0383217, decimal=6)


def test_brier_score_loss():
    # Check brier_score_loss function
    y_true = np.array([0, 1, 1, 0, 1, 1])
    y_pred = np.array([0.1, 0.8, 0.9, 0.3, 1., 0.95])
    true_score = linalg.norm(y_true - y_pred) ** 2 / len(y_true)

    assert_almost_equal(brier_score_loss(y_true, y_true), 0.0)
    assert_almost_equal(brier_score_loss(y_true, y_pred), true_score)
    assert_almost_equal(brier_score_loss(1. + y_true, y_pred),
                        true_score)
    assert_almost_equal(brier_score_loss(2 * y_true - 1, y_pred),
                        true_score)
    assert_raises(ValueError, brier_score_loss, y_true, y_pred[1:])
    assert_raises(ValueError, brier_score_loss, y_true, y_pred + 1.)
    assert_raises(ValueError, brier_score_loss, y_true, y_pred - 1.)
