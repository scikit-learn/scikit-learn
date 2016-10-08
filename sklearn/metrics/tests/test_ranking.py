from __future__ import division, print_function

import numpy as np
from itertools import product
import warnings
from scipy.sparse import csr_matrix

from sklearn import datasets
from sklearn import svm

from sklearn.datasets import make_multilabel_classification
from sklearn.random_projection import sparse_random_matrix
from sklearn.utils.validation import check_array, check_consistent_length
from sklearn.utils.validation import check_random_state

from sklearn.utils.testing import assert_raises, clean_warning_registry
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_warns

from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import coverage_error
from sklearn.metrics import label_ranking_average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import label_ranking_loss
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

from sklearn.exceptions import UndefinedMetricWarning


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

def _auc(y_true, y_score):
    """Alternative implementation to check for correctness of
    `roc_auc_score`."""
    pos_label = np.unique(y_true)[1]

    # Count the number of times positive samples are correctly ranked above
    # negative samples.
    pos = y_score[y_true == pos_label]
    neg = y_score[y_true != pos_label]
    diff_matrix = pos.reshape(1, -1) - neg.reshape(-1, 1)
    n_correct = np.sum(diff_matrix > 0)

    return n_correct / float(len(pos) * len(neg))


def _average_precision(y_true, y_score):
    """Alternative implementation to check for correctness of
    `average_precision_score`."""
    pos_label = np.unique(y_true)[1]
    n_pos = np.sum(y_true == pos_label)
    order = np.argsort(y_score)[::-1]
    y_score = y_score[order]
    y_true = y_true[order]

    score = 0
    for i in range(len(y_score)):
        if y_true[i] == pos_label:
            # Compute precision up to document i
            # i.e, percentage of relevant documents up to document i.
            prec = 0
            for j in range(0, i + 1):
                if y_true[j] == pos_label:
                    prec += 1.0
            prec /= (i + 1.0)
            score += prec

    return score / n_pos


def test_roc_curve():
    # Test Area under Receiver Operating Characteristic (ROC) curve
    y_true, _, probas_pred = make_prediction(binary=True)
    expected_auc = _auc(y_true, probas_pred)

    for drop in [True, False]:
        fpr, tpr, thresholds = roc_curve(y_true, probas_pred,
                                         drop_intermediate=drop)
        roc_auc = auc(fpr, tpr)
        assert_array_almost_equal(roc_auc, expected_auc, decimal=2)
        assert_almost_equal(roc_auc, roc_auc_score(y_true, probas_pred))
        assert_equal(fpr.shape, tpr.shape)
        assert_equal(fpr.shape, thresholds.shape)


def test_roc_curve_end_points():
    # Make sure that roc_curve returns a curve start at 0 and ending and
    # 1 even in corner cases
    rng = np.random.RandomState(0)
    y_true = np.array([0] * 50 + [1] * 50)
    y_pred = rng.randint(3, size=100)
    fpr, tpr, thr = roc_curve(y_true, y_pred, drop_intermediate=True)
    assert_equal(fpr[0], 0)
    assert_equal(fpr[-1], 1)
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thr.shape)


def test_roc_returns_consistency():
    # Test whether the returned threshold matches up with tpr
    # make small toy dataset
    y_true, _, probas_pred = make_prediction(binary=True)
    fpr, tpr, thresholds = roc_curve(y_true, probas_pred)

    # use the given thresholds to determine the tpr
    tpr_correct = []
    for t in thresholds:
        tp = np.sum((probas_pred >= t) & y_true)
        p = np.sum(y_true)
        tpr_correct.append(1.0 * tp / p)

    # compare tpr and tpr_correct to see if the thresholds' order was correct
    assert_array_almost_equal(tpr, tpr_correct, decimal=2)
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)


def test_roc_curve_multi():
    # roc_curve not applicable for multi-class problems
    y_true, _, probas_pred = make_prediction(binary=False)

    assert_raises(ValueError, roc_curve, y_true, probas_pred)


def test_roc_curve_confidence():
    # roc_curve for confidence scores
    y_true, _, probas_pred = make_prediction(binary=True)

    fpr, tpr, thresholds = roc_curve(y_true, probas_pred - 0.5)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.90, decimal=2)
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)


def test_roc_curve_hard():
    # roc_curve for hard decisions
    y_true, pred, probas_pred = make_prediction(binary=True)

    # always predict one
    trivial_pred = np.ones(y_true.shape)
    fpr, tpr, thresholds = roc_curve(y_true, trivial_pred)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.50, decimal=2)
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)

    # always predict zero
    trivial_pred = np.zeros(y_true.shape)
    fpr, tpr, thresholds = roc_curve(y_true, trivial_pred)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.50, decimal=2)
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)

    # hard decisions
    fpr, tpr, thresholds = roc_curve(y_true, pred)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.78, decimal=2)
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)


def test_roc_curve_one_label():
    y_true = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    y_pred = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    # assert there are warnings
    w = UndefinedMetricWarning
    fpr, tpr, thresholds = assert_warns(w, roc_curve, y_true, y_pred)
    # all true labels, all fpr should be nan
    assert_array_equal(fpr,
                       np.nan * np.ones(len(thresholds)))
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)

    # assert there are warnings
    fpr, tpr, thresholds = assert_warns(w, roc_curve,
                                        [1 - x for x in y_true],
                                        y_pred)
    # all negative labels, all tpr should be nan
    assert_array_equal(tpr,
                       np.nan * np.ones(len(thresholds)))
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)


def test_roc_curve_toydata():
    # Binary classification
    y_true = [0, 1]
    y_score = [0, 1]
    tpr, fpr, _ = roc_curve(y_true, y_score)
    roc_auc = roc_auc_score(y_true, y_score)
    assert_array_almost_equal(tpr, [0, 1])
    assert_array_almost_equal(fpr, [1, 1])
    assert_almost_equal(roc_auc, 1.)

    y_true = [0, 1]
    y_score = [1, 0]
    tpr, fpr, _ = roc_curve(y_true, y_score)
    roc_auc = roc_auc_score(y_true, y_score)
    assert_array_almost_equal(tpr, [0, 1, 1])
    assert_array_almost_equal(fpr, [0, 0, 1])
    assert_almost_equal(roc_auc, 0.)

    y_true = [1, 0]
    y_score = [1, 1]
    tpr, fpr, _ = roc_curve(y_true, y_score)
    roc_auc = roc_auc_score(y_true, y_score)
    assert_array_almost_equal(tpr, [0, 1])
    assert_array_almost_equal(fpr, [0, 1])
    assert_almost_equal(roc_auc, 0.5)

    y_true = [1, 0]
    y_score = [1, 0]
    tpr, fpr, _ = roc_curve(y_true, y_score)
    roc_auc = roc_auc_score(y_true, y_score)
    assert_array_almost_equal(tpr, [0, 1])
    assert_array_almost_equal(fpr, [1, 1])
    assert_almost_equal(roc_auc, 1.)

    y_true = [1, 0]
    y_score = [0.5, 0.5]
    tpr, fpr, _ = roc_curve(y_true, y_score)
    roc_auc = roc_auc_score(y_true, y_score)
    assert_array_almost_equal(tpr, [0, 1])
    assert_array_almost_equal(fpr, [0, 1])
    assert_almost_equal(roc_auc, .5)

    y_true = [0, 0]
    y_score = [0.25, 0.75]
    # assert UndefinedMetricWarning because of no positive sample in y_true
    tpr, fpr, _ = assert_warns(UndefinedMetricWarning, roc_curve, y_true, y_score)
    assert_raises(ValueError, roc_auc_score, y_true, y_score)
    assert_array_almost_equal(tpr, [0., 0.5, 1.])
    assert_array_almost_equal(fpr, [np.nan, np.nan, np.nan])

    y_true = [1, 1]
    y_score = [0.25, 0.75]
    # assert UndefinedMetricWarning because of no negative sample in y_true
    tpr, fpr, _ = assert_warns(UndefinedMetricWarning, roc_curve, y_true, y_score)
    assert_raises(ValueError, roc_auc_score, y_true, y_score)
    assert_array_almost_equal(tpr, [np.nan, np.nan])
    assert_array_almost_equal(fpr, [0.5, 1.])

    # Multi-label classification task
    y_true = np.array([[0, 1], [0, 1]])
    y_score = np.array([[0, 1], [0, 1]])
    assert_raises(ValueError, roc_auc_score, y_true, y_score, average="macro")
    assert_raises(ValueError, roc_auc_score, y_true, y_score,
                  average="weighted")
    assert_almost_equal(roc_auc_score(y_true, y_score, average="samples"), 1.)
    assert_almost_equal(roc_auc_score(y_true, y_score, average="micro"), 1.)

    y_true = np.array([[0, 1], [0, 1]])
    y_score = np.array([[0, 1], [1, 0]])
    assert_raises(ValueError, roc_auc_score, y_true, y_score, average="macro")
    assert_raises(ValueError, roc_auc_score, y_true, y_score,
                  average="weighted")
    assert_almost_equal(roc_auc_score(y_true, y_score, average="samples"), 0.5)
    assert_almost_equal(roc_auc_score(y_true, y_score, average="micro"), 0.5)

    y_true = np.array([[1, 0], [0, 1]])
    y_score = np.array([[0, 1], [1, 0]])
    assert_almost_equal(roc_auc_score(y_true, y_score, average="macro"), 0)
    assert_almost_equal(roc_auc_score(y_true, y_score, average="weighted"), 0)
    assert_almost_equal(roc_auc_score(y_true, y_score, average="samples"), 0)
    assert_almost_equal(roc_auc_score(y_true, y_score, average="micro"), 0)

    y_true = np.array([[1, 0], [0, 1]])
    y_score = np.array([[0.5, 0.5], [0.5, 0.5]])
    assert_almost_equal(roc_auc_score(y_true, y_score, average="macro"), .5)
    assert_almost_equal(roc_auc_score(y_true, y_score, average="weighted"), .5)
    assert_almost_equal(roc_auc_score(y_true, y_score, average="samples"), .5)
    assert_almost_equal(roc_auc_score(y_true, y_score, average="micro"), .5)


def test_roc_curve_drop_intermediate():
    # Test that drop_intermediate drops the correct thresholds
    y_true = [0, 0, 0, 0, 1, 1]
    y_score = [0., 0.2, 0.5, 0.6, 0.7, 1.0]
    tpr, fpr, thresholds = roc_curve(y_true, y_score, drop_intermediate=True)
    assert_array_almost_equal(thresholds, [1., 0.7, 0.])

    # Test dropping thresholds with repeating scores
    y_true = [0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 1, 1, 1]
    y_score = [0., 0.1, 0.6, 0.6, 0.7, 0.8, 0.9,
               0.6, 0.7, 0.8, 0.9, 0.9, 1.0]
    tpr, fpr, thresholds = roc_curve(y_true, y_score, drop_intermediate=True)
    assert_array_almost_equal(thresholds,
                              [1.0, 0.9, 0.7, 0.6, 0.])


def test_auc():
    # Test Area Under Curve (AUC) computation
    x = [0, 1]
    y = [0, 1]
    assert_array_almost_equal(auc(x, y), 0.5)
    x = [1, 0]
    y = [0, 1]
    assert_array_almost_equal(auc(x, y), 0.5)
    x = [1, 0, 0]
    y = [0, 1, 1]
    assert_array_almost_equal(auc(x, y), 0.5)
    x = [0, 1]
    y = [1, 1]
    assert_array_almost_equal(auc(x, y), 1)
    x = [0, 0.5, 1]
    y = [0, 0.5, 1]
    assert_array_almost_equal(auc(x, y), 0.5)


def test_auc_duplicate_values():
    # Test Area Under Curve (AUC) computation with duplicate values

    # auc() was previously sorting the x and y arrays according to the indices
    # from numpy.argsort(x), which was reordering the tied 0's in this example
    # and resulting in an incorrect area computation. This test detects the
    # error.
    x = [-2.0, 0.0, 0.0, 0.0, 1.0]
    y1 = [2.0, 0.0, 0.5, 1.0, 1.0]
    y2 = [2.0, 1.0, 0.0, 0.5, 1.0]
    y3 = [2.0, 1.0, 0.5, 0.0, 1.0]

    for y in (y1, y2, y3):
        assert_array_almost_equal(auc(x, y, reorder=True), 3.0)


def test_auc_errors():
    # Incompatible shapes
    assert_raises(ValueError, auc, [0.0, 0.5, 1.0], [0.1, 0.2])

    # Too few x values
    assert_raises(ValueError, auc, [0.0], [0.1])

    # x is not in order
    assert_raises(ValueError, auc, [1.0, 0.0, 0.5], [0.0, 0.0, 0.0])


def test_auc_score_non_binary_class():
    # Test that roc_auc_score function returns an error when trying
    # to compute AUC for non-binary class values.
    rng = check_random_state(404)
    y_pred = rng.rand(10)
    # y_true contains only one class value
    y_true = np.zeros(10, dtype="int")
    assert_raise_message(ValueError, "ROC AUC score is not defined",
                         roc_auc_score, y_true, y_pred)
    y_true = np.ones(10, dtype="int")
    assert_raise_message(ValueError, "ROC AUC score is not defined",
                         roc_auc_score, y_true, y_pred)
    y_true = -np.ones(10, dtype="int")
    assert_raise_message(ValueError, "ROC AUC score is not defined",
                         roc_auc_score, y_true, y_pred)
    # y_true contains three different class values
    y_true = rng.randint(0, 3, size=10)
    assert_raise_message(ValueError, "multiclass format is not supported",
                         roc_auc_score, y_true, y_pred)

    clean_warning_registry()
    with warnings.catch_warnings(record=True):
        rng = check_random_state(404)
        y_pred = rng.rand(10)
        # y_true contains only one class value
        y_true = np.zeros(10, dtype="int")
        assert_raise_message(ValueError, "ROC AUC score is not defined",
                             roc_auc_score, y_true, y_pred)
        y_true = np.ones(10, dtype="int")
        assert_raise_message(ValueError, "ROC AUC score is not defined",
                             roc_auc_score, y_true, y_pred)
        y_true = -np.ones(10, dtype="int")
        assert_raise_message(ValueError, "ROC AUC score is not defined",
                             roc_auc_score, y_true, y_pred)

        # y_true contains three different class values
        y_true = rng.randint(0, 3, size=10)
        assert_raise_message(ValueError, "multiclass format is not supported",
                             roc_auc_score, y_true, y_pred)


def test_precision_recall_curve():
    y_true, _, probas_pred = make_prediction(binary=True)
    _test_precision_recall_curve(y_true, probas_pred)

    # Use {-1, 1} for labels; make sure original labels aren't modified
    y_true[np.where(y_true == 0)] = -1
    y_true_copy = y_true.copy()
    _test_precision_recall_curve(y_true, probas_pred)
    assert_array_equal(y_true_copy, y_true)

    labels = [1, 0, 0, 1]
    predict_probas = [1, 2, 3, 4]
    p, r, t = precision_recall_curve(labels, predict_probas)
    assert_array_almost_equal(p, np.array([0.5, 0.33333333, 0.5, 1., 1.]))
    assert_array_almost_equal(r, np.array([1., 0.5, 0.5, 0.5, 0.]))
    assert_array_almost_equal(t, np.array([1, 2, 3, 4]))
    assert_equal(p.size, r.size)
    assert_equal(p.size, t.size + 1)


def test_precision_recall_curve_pos_label():
    y_true, _, probas_pred = make_prediction(binary=False)
    pos_label = 2
    p, r, thresholds = precision_recall_curve(y_true,
                                              probas_pred[:, pos_label],
                                              pos_label=pos_label)
    p2, r2, thresholds2 = precision_recall_curve(y_true == pos_label,
                                                 probas_pred[:, pos_label])
    assert_array_almost_equal(p, p2)
    assert_array_almost_equal(r, r2)
    assert_array_almost_equal(thresholds, thresholds2)
    assert_equal(p.size, r.size)
    assert_equal(p.size, thresholds.size + 1)


def _test_precision_recall_curve(y_true, probas_pred):
    # Test Precision-Recall and aread under PR curve
    p, r, thresholds = precision_recall_curve(y_true, probas_pred)
    precision_recall_auc = auc(r, p)
    assert_array_almost_equal(precision_recall_auc, 0.85, 2)
    assert_array_almost_equal(precision_recall_auc,
                              average_precision_score(y_true, probas_pred))
    assert_almost_equal(_average_precision(y_true, probas_pred),
                        precision_recall_auc, 1)
    assert_equal(p.size, r.size)
    assert_equal(p.size, thresholds.size + 1)
    # Smoke test in the case of proba having only one value
    p, r, thresholds = precision_recall_curve(y_true,
                                              np.zeros_like(probas_pred))
    precision_recall_auc = auc(r, p)
    assert_array_almost_equal(precision_recall_auc, 0.75, 3)
    assert_equal(p.size, r.size)
    assert_equal(p.size, thresholds.size + 1)


def test_precision_recall_curve_errors():
    # Contains non-binary labels
    assert_raises(ValueError, precision_recall_curve,
                  [0, 1, 2], [[0.0], [1.0], [1.0]])


def test_precision_recall_curve_toydata():
    with np.errstate(all="raise"):
        # Binary classification
        y_true = [0, 1]
        y_score = [0, 1]
        p, r, _ = precision_recall_curve(y_true, y_score)
        auc_prc = average_precision_score(y_true, y_score)
        assert_array_almost_equal(p, [1, 1])
        assert_array_almost_equal(r, [1, 0])
        assert_almost_equal(auc_prc, 1.)

        y_true = [0, 1]
        y_score = [1, 0]
        p, r, _ = precision_recall_curve(y_true, y_score)
        auc_prc = average_precision_score(y_true, y_score)
        assert_array_almost_equal(p, [0.5, 0., 1.])
        assert_array_almost_equal(r, [1., 0.,  0.])
        assert_almost_equal(auc_prc, 0.25)

        y_true = [1, 0]
        y_score = [1, 1]
        p, r, _ = precision_recall_curve(y_true, y_score)
        auc_prc = average_precision_score(y_true, y_score)
        assert_array_almost_equal(p, [0.5, 1])
        assert_array_almost_equal(r, [1., 0])
        assert_almost_equal(auc_prc, .75)

        y_true = [1, 0]
        y_score = [1, 0]
        p, r, _ = precision_recall_curve(y_true, y_score)
        auc_prc = average_precision_score(y_true, y_score)
        assert_array_almost_equal(p, [1, 1])
        assert_array_almost_equal(r, [1, 0])
        assert_almost_equal(auc_prc, 1.)

        y_true = [1, 0]
        y_score = [0.5, 0.5]
        p, r, _ = precision_recall_curve(y_true, y_score)
        auc_prc = average_precision_score(y_true, y_score)
        assert_array_almost_equal(p, [0.5, 1])
        assert_array_almost_equal(r, [1, 0.])
        assert_almost_equal(auc_prc, .75)

        y_true = [0, 0]
        y_score = [0.25, 0.75]
        assert_raises(Exception, precision_recall_curve, y_true, y_score)
        assert_raises(Exception, average_precision_score, y_true, y_score)

        y_true = [1, 1]
        y_score = [0.25, 0.75]
        p, r, _ = precision_recall_curve(y_true, y_score)
        assert_almost_equal(average_precision_score(y_true, y_score), 1.)
        assert_array_almost_equal(p, [1., 1., 1.])
        assert_array_almost_equal(r, [1, 0.5, 0.])

        # Multi-label classification task
        y_true = np.array([[0, 1], [0, 1]])
        y_score = np.array([[0, 1], [0, 1]])
        assert_raises(Exception, average_precision_score, y_true, y_score,
                      average="macro")
        assert_raises(Exception, average_precision_score, y_true, y_score,
                      average="weighted")
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="samples"), 1.)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="micro"), 1.)

        y_true = np.array([[0, 1], [0, 1]])
        y_score = np.array([[0, 1], [1, 0]])
        assert_raises(Exception, average_precision_score, y_true, y_score,
                      average="macro")
        assert_raises(Exception, average_precision_score, y_true, y_score,
                      average="weighted")
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="samples"), 0.625)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="micro"), 0.625)

        y_true = np.array([[1, 0], [0, 1]])
        y_score = np.array([[0, 1], [1, 0]])
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="macro"), 0.25)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="weighted"), 0.25)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="samples"), 0.25)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="micro"), 0.25)

        y_true = np.array([[1, 0], [0, 1]])
        y_score = np.array([[0.5, 0.5], [0.5, 0.5]])
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="macro"), 0.75)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="weighted"), 0.75)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="samples"), 0.75)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="micro"), 0.75)


def test_score_scale_invariance():
    # Test that average_precision_score and roc_auc_score are invariant by
    # the scaling or shifting of probabilities
    # This test was expanded (added scaled_down) in response to github
    # issue #3864 (and others), where overly aggressive rounding was causing
    # problems for users with very small y_score values
    y_true, _, probas_pred = make_prediction(binary=True)

    roc_auc = roc_auc_score(y_true, probas_pred)
    roc_auc_scaled_up = roc_auc_score(y_true, 100 * probas_pred)
    roc_auc_scaled_down = roc_auc_score(y_true, 1e-6 * probas_pred)
    roc_auc_shifted = roc_auc_score(y_true, probas_pred - 10)
    assert_equal(roc_auc, roc_auc_scaled_up)
    assert_equal(roc_auc, roc_auc_scaled_down)
    assert_equal(roc_auc, roc_auc_shifted)

    pr_auc = average_precision_score(y_true, probas_pred)
    pr_auc_scaled_up = average_precision_score(y_true, 100 * probas_pred)
    pr_auc_scaled_down = average_precision_score(y_true, 1e-6 * probas_pred)
    pr_auc_shifted = average_precision_score(y_true, probas_pred - 10)
    assert_equal(pr_auc, pr_auc_scaled_up)
    assert_equal(pr_auc, pr_auc_scaled_down)
    assert_equal(pr_auc, pr_auc_shifted)


def check_lrap_toy(lrap_score):
    # Check on several small example that it works
    assert_almost_equal(lrap_score([[0, 1]], [[0.25, 0.75]]), 1)
    assert_almost_equal(lrap_score([[0, 1]], [[0.75, 0.25]]), 1 / 2)
    assert_almost_equal(lrap_score([[1, 1]], [[0.75, 0.25]]), 1)

    assert_almost_equal(lrap_score([[0, 0, 1]], [[0.25, 0.5, 0.75]]), 1)
    assert_almost_equal(lrap_score([[0, 1, 0]], [[0.25, 0.5, 0.75]]), 1 / 2)
    assert_almost_equal(lrap_score([[0, 1, 1]], [[0.25, 0.5, 0.75]]), 1)
    assert_almost_equal(lrap_score([[1, 0, 0]], [[0.25, 0.5, 0.75]]), 1 / 3)
    assert_almost_equal(lrap_score([[1, 0, 1]], [[0.25, 0.5, 0.75]]),
                        (2 / 3 + 1 / 1) / 2)
    assert_almost_equal(lrap_score([[1, 1, 0]], [[0.25, 0.5, 0.75]]),
                        (2 / 3 + 1 / 2) / 2)

    assert_almost_equal(lrap_score([[0, 0, 1]], [[0.75, 0.5, 0.25]]), 1 / 3)
    assert_almost_equal(lrap_score([[0, 1, 0]], [[0.75, 0.5, 0.25]]), 1 / 2)
    assert_almost_equal(lrap_score([[0, 1, 1]], [[0.75, 0.5, 0.25]]),
                        (1 / 2 + 2 / 3) / 2)
    assert_almost_equal(lrap_score([[1, 0, 0]], [[0.75, 0.5, 0.25]]), 1)
    assert_almost_equal(lrap_score([[1, 0, 1]], [[0.75, 0.5, 0.25]]),
                        (1 + 2 / 3) / 2)
    assert_almost_equal(lrap_score([[1, 1, 0]], [[0.75, 0.5, 0.25]]), 1)
    assert_almost_equal(lrap_score([[1, 1, 1]], [[0.75, 0.5, 0.25]]), 1)

    assert_almost_equal(lrap_score([[0, 0, 1]], [[0.5, 0.75, 0.25]]), 1 / 3)
    assert_almost_equal(lrap_score([[0, 1, 0]], [[0.5, 0.75, 0.25]]), 1)
    assert_almost_equal(lrap_score([[0, 1, 1]], [[0.5, 0.75, 0.25]]),
                        (1 + 2 / 3) / 2)
    assert_almost_equal(lrap_score([[1, 0, 0]], [[0.5, 0.75, 0.25]]), 1 / 2)
    assert_almost_equal(lrap_score([[1, 0, 1]], [[0.5, 0.75, 0.25]]),
                        (1 / 2 + 2 / 3) / 2)
    assert_almost_equal(lrap_score([[1, 1, 0]], [[0.5, 0.75, 0.25]]), 1)
    assert_almost_equal(lrap_score([[1, 1, 1]], [[0.5, 0.75, 0.25]]), 1)

    # Tie handling
    assert_almost_equal(lrap_score([[1, 0]], [[0.5, 0.5]]), 0.5)
    assert_almost_equal(lrap_score([[0, 1]], [[0.5, 0.5]]), 0.5)
    assert_almost_equal(lrap_score([[1, 1]], [[0.5, 0.5]]), 1)

    assert_almost_equal(lrap_score([[0, 0, 1]], [[0.25, 0.5, 0.5]]), 0.5)
    assert_almost_equal(lrap_score([[0, 1, 0]], [[0.25, 0.5, 0.5]]), 0.5)
    assert_almost_equal(lrap_score([[0, 1, 1]], [[0.25, 0.5, 0.5]]), 1)
    assert_almost_equal(lrap_score([[1, 0, 0]], [[0.25, 0.5, 0.5]]), 1 / 3)
    assert_almost_equal(lrap_score([[1, 0, 1]], [[0.25, 0.5, 0.5]]),
                        (2 / 3 + 1 / 2) / 2)
    assert_almost_equal(lrap_score([[1, 1, 0]], [[0.25, 0.5, 0.5]]),
                        (2 / 3 + 1 / 2) / 2)
    assert_almost_equal(lrap_score([[1, 1, 1]], [[0.25, 0.5, 0.5]]), 1)

    assert_almost_equal(lrap_score([[1, 1, 0]], [[0.5, 0.5, 0.5]]), 2 / 3)

    assert_almost_equal(lrap_score([[1, 1, 1, 0]], [[0.5, 0.5, 0.5, 0.5]]),
                        3 / 4)


def check_zero_or_all_relevant_labels(lrap_score):
    random_state = check_random_state(0)

    for n_labels in range(2, 5):
        y_score = random_state.uniform(size=(1, n_labels))
        y_score_ties = np.zeros_like(y_score)

        # No relevant labels
        y_true = np.zeros((1, n_labels))
        assert_equal(lrap_score(y_true, y_score), 1.)
        assert_equal(lrap_score(y_true, y_score_ties), 1.)

        # Only relevant labels
        y_true = np.ones((1, n_labels))
        assert_equal(lrap_score(y_true, y_score), 1.)
        assert_equal(lrap_score(y_true, y_score_ties), 1.)

    # Degenerate case: only one label
    assert_almost_equal(lrap_score([[1], [0], [1], [0]],
                                   [[0.5], [0.5], [0.5], [0.5]]), 1.)


def check_lrap_error_raised(lrap_score):
    # Raise value error if not appropriate format
    assert_raises(ValueError, lrap_score,
                  [0, 1, 0], [0.25, 0.3, 0.2])
    assert_raises(ValueError, lrap_score, [0, 1, 2],
                  [[0.25, 0.75, 0.0], [0.7, 0.3, 0.0], [0.8, 0.2, 0.0]])
    assert_raises(ValueError, lrap_score, [(0), (1), (2)],
                  [[0.25, 0.75, 0.0], [0.7, 0.3, 0.0], [0.8, 0.2, 0.0]])

    # Check that y_true.shape != y_score.shape raise the proper exception
    assert_raises(ValueError, lrap_score, [[0, 1], [0, 1]], [0, 1])
    assert_raises(ValueError, lrap_score, [[0, 1], [0, 1]], [[0, 1]])
    assert_raises(ValueError, lrap_score, [[0, 1], [0, 1]], [[0], [1]])
    assert_raises(ValueError, lrap_score, [[0, 1]], [[0, 1], [0, 1]])
    assert_raises(ValueError, lrap_score, [[0], [1]], [[0, 1], [0, 1]])
    assert_raises(ValueError, lrap_score, [[0, 1], [0, 1]], [[0], [1]])


def check_lrap_only_ties(lrap_score):
    # Check tie handling in score
    # Basic check with only ties and increasing label space
    for n_labels in range(2, 10):
        y_score = np.ones((1, n_labels))

        # Check for growing number of consecutive relevant
        for n_relevant in range(1, n_labels):
            # Check for a bunch of positions
            for pos in range(n_labels - n_relevant):
                y_true = np.zeros((1, n_labels))
                y_true[0, pos:pos + n_relevant] = 1
                assert_almost_equal(lrap_score(y_true, y_score),
                                    n_relevant / n_labels)


def check_lrap_without_tie_and_increasing_score(lrap_score):
    # Check that Label ranking average precision works for various
    # Basic check with increasing label space size and decreasing score
    for n_labels in range(2, 10):
        y_score = n_labels - (np.arange(n_labels).reshape((1, n_labels)) + 1)

        # First and last
        y_true = np.zeros((1, n_labels))
        y_true[0, 0] = 1
        y_true[0, -1] = 1
        assert_almost_equal(lrap_score(y_true, y_score),
                            (2 / n_labels + 1) / 2)

        # Check for growing number of consecutive relevant label
        for n_relevant in range(1, n_labels):
            # Check for a bunch of position
            for pos in range(n_labels - n_relevant):
                y_true = np.zeros((1, n_labels))
                y_true[0, pos:pos + n_relevant] = 1
                assert_almost_equal(lrap_score(y_true, y_score),
                                    sum((r + 1) / ((pos + r + 1) * n_relevant)
                                        for r in range(n_relevant)))


def _my_lrap(y_true, y_score):
    """Simple implementation of label ranking average precision"""
    check_consistent_length(y_true, y_score)
    y_true = check_array(y_true)
    y_score = check_array(y_score)
    n_samples, n_labels = y_true.shape
    score = np.empty((n_samples, ))
    for i in range(n_samples):
        # The best rank correspond to 1. Rank higher than 1 are worse.
        # The best inverse ranking correspond to n_labels.
        unique_rank, inv_rank = np.unique(y_score[i], return_inverse=True)
        n_ranks = unique_rank.size
        rank = n_ranks - inv_rank

        # Rank need to be corrected to take into account ties
        # ex: rank 1 ex aequo means that both label are rank 2.
        corr_rank = np.bincount(rank, minlength=n_ranks + 1).cumsum()
        rank = corr_rank[rank]

        relevant = y_true[i].nonzero()[0]
        if relevant.size == 0 or relevant.size == n_labels:
            score[i] = 1
            continue

        score[i] = 0.
        for label in relevant:
            # Let's count the number of relevant label with better rank
            # (smaller rank).
            n_ranked_above = sum(rank[r] <= rank[label] for r in relevant)

            # Weight by the rank of the actual label
            score[i] += n_ranked_above / rank[label]

        score[i] /= relevant.size

    return score.mean()


def check_alternative_lrap_implementation(lrap_score, n_classes=5,
                                          n_samples=20, random_state=0):
    _, y_true = make_multilabel_classification(n_features=1,
                                               allow_unlabeled=False,
                                               random_state=random_state,
                                               n_classes=n_classes,
                                               n_samples=n_samples)

    # Score with ties
    y_score = sparse_random_matrix(n_components=y_true.shape[0],
                                   n_features=y_true.shape[1],
                                   random_state=random_state)

    if hasattr(y_score, "toarray"):
        y_score = y_score.toarray()
    score_lrap = label_ranking_average_precision_score(y_true, y_score)
    score_my_lrap = _my_lrap(y_true, y_score)
    assert_almost_equal(score_lrap, score_my_lrap)

    # Uniform score
    random_state = check_random_state(random_state)
    y_score = random_state.uniform(size=(n_samples, n_classes))
    score_lrap = label_ranking_average_precision_score(y_true, y_score)
    score_my_lrap = _my_lrap(y_true, y_score)
    assert_almost_equal(score_lrap, score_my_lrap)


def test_label_ranking_avp():
    for fn in [label_ranking_average_precision_score, _my_lrap]:
        yield check_lrap_toy, fn
        yield check_lrap_without_tie_and_increasing_score, fn
        yield check_lrap_only_ties, fn
        yield check_zero_or_all_relevant_labels, fn
        yield check_lrap_error_raised, label_ranking_average_precision_score

    for n_samples, n_classes, random_state in product((1, 2, 8, 20),
                                                      (2, 5, 10),
                                                      range(1)):
        yield (check_alternative_lrap_implementation,
               label_ranking_average_precision_score,
               n_classes, n_samples, random_state)


def test_coverage_error():
    # Toy case
    assert_almost_equal(coverage_error([[0, 1]], [[0.25, 0.75]]), 1)
    assert_almost_equal(coverage_error([[0, 1]], [[0.75, 0.25]]), 2)
    assert_almost_equal(coverage_error([[1, 1]], [[0.75, 0.25]]), 2)
    assert_almost_equal(coverage_error([[0, 0]], [[0.75, 0.25]]), 0)

    assert_almost_equal(coverage_error([[0, 0, 0]], [[0.25, 0.5, 0.75]]), 0)
    assert_almost_equal(coverage_error([[0, 0, 1]], [[0.25, 0.5, 0.75]]), 1)
    assert_almost_equal(coverage_error([[0, 1, 0]], [[0.25, 0.5, 0.75]]), 2)
    assert_almost_equal(coverage_error([[0, 1, 1]], [[0.25, 0.5, 0.75]]), 2)
    assert_almost_equal(coverage_error([[1, 0, 0]], [[0.25, 0.5, 0.75]]), 3)
    assert_almost_equal(coverage_error([[1, 0, 1]], [[0.25, 0.5, 0.75]]), 3)
    assert_almost_equal(coverage_error([[1, 1, 0]], [[0.25, 0.5, 0.75]]), 3)
    assert_almost_equal(coverage_error([[1, 1, 1]], [[0.25, 0.5, 0.75]]), 3)

    assert_almost_equal(coverage_error([[0, 0, 0]], [[0.75, 0.5, 0.25]]), 0)
    assert_almost_equal(coverage_error([[0, 0, 1]], [[0.75, 0.5, 0.25]]), 3)
    assert_almost_equal(coverage_error([[0, 1, 0]], [[0.75, 0.5, 0.25]]), 2)
    assert_almost_equal(coverage_error([[0, 1, 1]], [[0.75, 0.5, 0.25]]), 3)
    assert_almost_equal(coverage_error([[1, 0, 0]], [[0.75, 0.5, 0.25]]), 1)
    assert_almost_equal(coverage_error([[1, 0, 1]], [[0.75, 0.5, 0.25]]), 3)
    assert_almost_equal(coverage_error([[1, 1, 0]], [[0.75, 0.5, 0.25]]), 2)
    assert_almost_equal(coverage_error([[1, 1, 1]], [[0.75, 0.5, 0.25]]), 3)

    assert_almost_equal(coverage_error([[0, 0, 0]], [[0.5, 0.75, 0.25]]), 0)
    assert_almost_equal(coverage_error([[0, 0, 1]], [[0.5, 0.75, 0.25]]), 3)
    assert_almost_equal(coverage_error([[0, 1, 0]], [[0.5, 0.75, 0.25]]), 1)
    assert_almost_equal(coverage_error([[0, 1, 1]], [[0.5, 0.75, 0.25]]), 3)
    assert_almost_equal(coverage_error([[1, 0, 0]], [[0.5, 0.75, 0.25]]), 2)
    assert_almost_equal(coverage_error([[1, 0, 1]], [[0.5, 0.75, 0.25]]), 3)
    assert_almost_equal(coverage_error([[1, 1, 0]], [[0.5, 0.75, 0.25]]), 2)
    assert_almost_equal(coverage_error([[1, 1, 1]], [[0.5, 0.75, 0.25]]), 3)

    # Non trival case
    assert_almost_equal(coverage_error([[0, 1, 0], [1, 1, 0]],
                                       [[0.1, 10., -3], [0, 1, 3]]),
                        (1 + 3) / 2.)

    assert_almost_equal(coverage_error([[0, 1, 0], [1, 1, 0], [0, 1, 1]],
                                       [[0.1, 10, -3], [0, 1, 3], [0, 2, 0]]),
                        (1 + 3 + 3) / 3.)

    assert_almost_equal(coverage_error([[0, 1, 0], [1, 1, 0], [0, 1, 1]],
                                       [[0.1, 10, -3], [3, 1, 3], [0, 2, 0]]),
                        (1 + 3 + 3) / 3.)


def test_coverage_tie_handling():
    assert_almost_equal(coverage_error([[0, 0]], [[0.5, 0.5]]), 0)
    assert_almost_equal(coverage_error([[1, 0]], [[0.5, 0.5]]), 2)
    assert_almost_equal(coverage_error([[0, 1]], [[0.5, 0.5]]), 2)
    assert_almost_equal(coverage_error([[1, 1]], [[0.5, 0.5]]), 2)

    assert_almost_equal(coverage_error([[0, 0, 0]], [[0.25, 0.5, 0.5]]), 0)
    assert_almost_equal(coverage_error([[0, 0, 1]], [[0.25, 0.5, 0.5]]), 2)
    assert_almost_equal(coverage_error([[0, 1, 0]], [[0.25, 0.5, 0.5]]), 2)
    assert_almost_equal(coverage_error([[0, 1, 1]], [[0.25, 0.5, 0.5]]), 2)
    assert_almost_equal(coverage_error([[1, 0, 0]], [[0.25, 0.5, 0.5]]), 3)
    assert_almost_equal(coverage_error([[1, 0, 1]], [[0.25, 0.5, 0.5]]), 3)
    assert_almost_equal(coverage_error([[1, 1, 0]], [[0.25, 0.5, 0.5]]), 3)
    assert_almost_equal(coverage_error([[1, 1, 1]], [[0.25, 0.5, 0.5]]), 3)


def test_label_ranking_loss():
    assert_almost_equal(label_ranking_loss([[0, 1]], [[0.25, 0.75]]), 0)
    assert_almost_equal(label_ranking_loss([[0, 1]], [[0.75, 0.25]]), 1)

    assert_almost_equal(label_ranking_loss([[0, 0, 1]], [[0.25, 0.5, 0.75]]),
                        0)
    assert_almost_equal(label_ranking_loss([[0, 1, 0]], [[0.25, 0.5, 0.75]]),
                        1 / 2)
    assert_almost_equal(label_ranking_loss([[0, 1, 1]], [[0.25, 0.5, 0.75]]),
                        0)
    assert_almost_equal(label_ranking_loss([[1, 0, 0]], [[0.25, 0.5, 0.75]]),
                        2 / 2)
    assert_almost_equal(label_ranking_loss([[1, 0, 1]], [[0.25, 0.5, 0.75]]),
                        1 / 2)
    assert_almost_equal(label_ranking_loss([[1, 1, 0]], [[0.25, 0.5, 0.75]]),
                        2 / 2)

    # Undefined metrics -  the ranking doesn't matter
    assert_almost_equal(label_ranking_loss([[0, 0]], [[0.75, 0.25]]), 0)
    assert_almost_equal(label_ranking_loss([[1, 1]], [[0.75, 0.25]]), 0)
    assert_almost_equal(label_ranking_loss([[0, 0]], [[0.5, 0.5]]), 0)
    assert_almost_equal(label_ranking_loss([[1, 1]], [[0.5, 0.5]]), 0)

    assert_almost_equal(label_ranking_loss([[0, 0, 0]], [[0.5, 0.75, 0.25]]),
                        0)
    assert_almost_equal(label_ranking_loss([[1, 1, 1]], [[0.5, 0.75, 0.25]]),
                        0)
    assert_almost_equal(label_ranking_loss([[0, 0, 0]], [[0.25, 0.5, 0.5]]),
                        0)
    assert_almost_equal(label_ranking_loss([[1, 1, 1]], [[0.25, 0.5, 0.5]]), 0)

    # Non trival case
    assert_almost_equal(label_ranking_loss([[0, 1, 0], [1, 1, 0]],
                                           [[0.1, 10., -3], [0, 1, 3]]),
                        (0 + 2 / 2) / 2.)

    assert_almost_equal(label_ranking_loss(
        [[0, 1, 0], [1, 1, 0], [0, 1, 1]],
        [[0.1, 10, -3], [0, 1, 3], [0, 2, 0]]),
        (0 + 2 / 2 + 1 / 2) / 3.)

    assert_almost_equal(label_ranking_loss(
        [[0, 1, 0], [1, 1, 0], [0, 1, 1]],
        [[0.1, 10, -3], [3, 1, 3], [0, 2, 0]]),
        (0 + 2 / 2 + 1 / 2) / 3.)

    # Sparse csr matrices
    assert_almost_equal(label_ranking_loss(
        csr_matrix(np.array([[0, 1, 0], [1, 1, 0]])),
        [[0.1, 10, -3], [3, 1, 3]]),
        (0 + 2 / 2) / 2.)


def test_ranking_appropriate_input_shape():
    # Check that y_true.shape != y_score.shape raise the proper exception
    assert_raises(ValueError, label_ranking_loss, [[0, 1], [0, 1]], [0, 1])
    assert_raises(ValueError, label_ranking_loss, [[0, 1], [0, 1]], [[0, 1]])
    assert_raises(ValueError, label_ranking_loss,
                  [[0, 1], [0, 1]], [[0], [1]])

    assert_raises(ValueError, label_ranking_loss, [[0, 1]], [[0, 1], [0, 1]])
    assert_raises(ValueError, label_ranking_loss,
                  [[0], [1]], [[0, 1], [0, 1]])
    assert_raises(ValueError, label_ranking_loss, [[0, 1], [0, 1]], [[0], [1]])


def test_ranking_loss_ties_handling():
    # Tie handling
    assert_almost_equal(label_ranking_loss([[1, 0]], [[0.5, 0.5]]), 1)
    assert_almost_equal(label_ranking_loss([[0, 1]], [[0.5, 0.5]]), 1)
    assert_almost_equal(label_ranking_loss([[0, 0, 1]], [[0.25, 0.5, 0.5]]),
                        1 / 2)
    assert_almost_equal(label_ranking_loss([[0, 1, 0]], [[0.25, 0.5, 0.5]]),
                        1 / 2)
    assert_almost_equal(label_ranking_loss([[0, 1, 1]], [[0.25, 0.5, 0.5]]), 0)
    assert_almost_equal(label_ranking_loss([[1, 0, 0]], [[0.25, 0.5, 0.5]]), 1)
    assert_almost_equal(label_ranking_loss([[1, 0, 1]], [[0.25, 0.5, 0.5]]), 1)
    assert_almost_equal(label_ranking_loss([[1, 1, 0]], [[0.25, 0.5, 0.5]]), 1)
