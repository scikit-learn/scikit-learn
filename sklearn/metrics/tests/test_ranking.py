from __future__ import division, print_function

import pytest
import numpy as np
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
from sklearn.utils.testing import assert_warns_message

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
    `average_precision_score`.

    Note that this implementation fails on some edge cases.
    For example, for constant predictions e.g. [0.5, 0.5, 0.5],
    y_true = [1, 0, 0] returns an average precision of 0.33...
    but y_true = [0, 0, 1] returns 1.0.
    """
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


def _average_precision_slow(y_true, y_score):
    """A second alternative implementation of average precision that closely
    follows the Wikipedia article's definition (see References). This should
    give identical results as `average_precision_score` for all inputs.

    References
    ----------
    .. [1] `Wikipedia entry for the Average precision
       <https://en.wikipedia.org/wiki/Average_precision>`_
    """
    precision, recall, threshold = precision_recall_curve(y_true, y_score)
    precision = list(reversed(precision))
    recall = list(reversed(recall))
    average_precision = 0
    for i in range(1, len(precision)):
        average_precision += precision[i] * (recall[i] - recall[i - 1])
    return average_precision


def _partial_roc_auc_score(y_true, y_predict, max_fpr):
    """Alternative implementation to check for correctness of `roc_auc_score`
    with `max_fpr` set.
    """

    def _partial_roc(y_true, y_predict, max_fpr):
        fpr, tpr, _ = roc_curve(y_true, y_predict)
        new_fpr = fpr[fpr <= max_fpr]
        new_fpr = np.append(new_fpr, max_fpr)
        new_tpr = tpr[fpr <= max_fpr]
        idx_out = np.argmax(fpr > max_fpr)
        idx_in = idx_out - 1
        x_interp = [fpr[idx_in], fpr[idx_out]]
        y_interp = [tpr[idx_in], tpr[idx_out]]
        new_tpr = np.append(new_tpr, np.interp(max_fpr, x_interp, y_interp))
        return (new_fpr, new_tpr)

    new_fpr, new_tpr = _partial_roc(y_true, y_predict, max_fpr)
    partial_auc = auc(new_fpr, new_tpr)

    # Formula (5) from McClish 1989
    fpr1 = 0
    fpr2 = max_fpr
    min_area = 0.5 * (fpr2 - fpr1) * (fpr2 + fpr1)
    max_area = fpr2 - fpr1
    return 0.5 * (1 + (partial_auc - min_area) / (max_area - min_area))


@pytest.mark.parametrize('drop', [True, False])
def test_roc_curve(drop):
    # Test Area under Receiver Operating Characteristic (ROC) curve
    y_true, _, probas_pred = make_prediction(binary=True)
    expected_auc = _auc(y_true, probas_pred)

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
    assert_array_equal(fpr, np.full(len(thresholds), np.nan))
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)

    # assert there are warnings
    fpr, tpr, thresholds = assert_warns(w, roc_curve,
                                        [1 - x for x in y_true],
                                        y_pred)
    # all negative labels, all tpr should be nan
    assert_array_equal(tpr, np.full(len(thresholds), np.nan))
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)


def test_roc_curve_toydata():
    # Binary classification
    y_true = [0, 1]
    y_score = [0, 1]
    tpr, fpr, _ = roc_curve(y_true, y_score)
    roc_auc = roc_auc_score(y_true, y_score)
    assert_array_almost_equal(tpr, [0, 0, 1])
    assert_array_almost_equal(fpr, [0, 1, 1])
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
    assert_array_almost_equal(tpr, [0, 0, 1])
    assert_array_almost_equal(fpr, [0, 1, 1])
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
    assert_array_almost_equal(tpr, [np.nan, np.nan, np.nan])
    assert_array_almost_equal(fpr, [0., 0.5, 1.])

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
    assert_array_almost_equal(thresholds, [2., 1., 0.7, 0.])

    # Test dropping thresholds with repeating scores
    y_true = [0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 1, 1, 1]
    y_score = [0., 0.1, 0.6, 0.6, 0.7, 0.8, 0.9,
               0.6, 0.7, 0.8, 0.9, 0.9, 1.0]
    tpr, fpr, thresholds = roc_curve(y_true, y_score, drop_intermediate=True)
    assert_array_almost_equal(thresholds,
                              [2.0, 1.0, 0.9, 0.7, 0.6, 0.])


def test_roc_curve_fpr_tpr_increasing():
    # Ensure that fpr and tpr returned by roc_curve are increasing.
    # Construct an edge case with float y_score and sample_weight
    # when some adjacent values of fpr and tpr are actually the same.
    y_true = [0, 0, 1, 1, 1]
    y_score = [0.1, 0.7, 0.3, 0.4, 0.5]
    sample_weight = np.repeat(0.2, 5)
    fpr, tpr, _ = roc_curve(y_true, y_score, sample_weight=sample_weight)
    assert_equal((np.diff(fpr) < 0).sum(), 0)
    assert_equal((np.diff(tpr) < 0).sum(), 0)


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


@pytest.mark.filterwarnings("ignore: The 'reorder' parameter")  # 0.22
def test_auc_duplicate_values():
    # Test Area Under Curve (AUC) computation with duplicate values

    # auc() was previously sorting the x and y arrays according to the indices
    # from numpy.argsort(x), which was reordering the tied 0's in this example
    # and resulting in an incorrect area computation. This test detects the
    # error.

    # This will not work again in the future! so regression?
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
    x = [2, 1, 3, 4]
    y = [5, 6, 7, 8]
    error_message = ("x is neither increasing nor decreasing : "
                     "{}".format(np.array(x)))
    assert_raise_message(ValueError, error_message, auc, x, y)


def test_deprecated_auc_reorder():
    depr_message = ("The 'reorder' parameter has been deprecated in version "
                    "0.20 and will be removed in 0.22. It is recommended not "
                    "to set 'reorder' and ensure that x is monotonic "
                    "increasing or monotonic decreasing.")
    assert_warns_message(DeprecationWarning, depr_message, auc,
                         [1, 2], [2, 3], reorder=True)


def test_multi_ovo_auc_toydata():
    # Tests the one-vs-one multiclass ROC AUC algorithm
    # on a small example, representative of an expected use case.
    y_true = np.array([0, 1, 0, 2])
    n_labels = len(np.unique(y_true))
    y_scores = np.array(
        [[0.1, 0.8, 0.1], [0.3, 0.4, 0.3], [0.35, 0.5, 0.15], [0, 0.2, 0.8]])

    # Used to compute the expected output.
    # Consider labels 0 and 1:
    # positive label is 0, negative label is 1
    score_01 = roc_auc_score([1, 0, 1], [0.1, 0.3, 0.35])
    # positive label is 1, negative label is 0
    score_10 = roc_auc_score([0, 1, 0], [0.8, 0.4, 0.5])
    average_score_01 = (score_01 + score_10) / 2.

    # Consider labels 0 and 2:
    score_02 = roc_auc_score([1, 1, 0], [0.1, 0.35, 0])
    score_20 = roc_auc_score([0, 0, 1], [0.1, 0.15, 0.8])
    average_score_02 = (score_02 + score_20) / 2.

    # Consider labels 1 and 2:
    score_12 = roc_auc_score([1, 0], [0.4, 0.2])
    score_21 = roc_auc_score([0, 1], [0.3, 0.8])
    average_score_12 = (score_12 + score_21) / 2.

    # Unweighted, one-vs-one multiclass ROC AUC algorithm
    sum_avg_scores = average_score_01 + average_score_02 + average_score_12
    ovo_unweighted_coefficient = 2. / (n_labels * (n_labels - 1))
    ovo_unweighted_score = ovo_unweighted_coefficient * sum_avg_scores
    assert_almost_equal(
        roc_auc_score(y_true, y_scores, multiclass="ovo"),
        ovo_unweighted_score)

    # Weighted, one-vs-one multiclass ROC AUC algorithm
    # Each term is weighted by the prevalence for the positive label.
    pair_scores = [average_score_01, average_score_02, average_score_12]
    prevalence = [0.75, 0.75, 0.50]
    ovo_weighted_score = np.average(pair_scores, weights=prevalence)
    assert_almost_equal(
        roc_auc_score(y_true, y_scores, multiclass="ovo", average="weighted"),
        ovo_weighted_score)


def test_multi_ovr_auc_toydata():
    # Tests the unweighted, one-vs-rest multiclass ROC AUC algorithm
    # on a small example, representative of an expected use case.
    y_true = np.array([0, 1, 2, 2])
    y_scores = np.array(
        [[1.0, 0.0, 0.0], [0.1, 0.5, 0.4], [0.1, 0.1, 0.8], [0.3, 0.3, 0.4]])
    # Compute the expected result by individually computing the 'one-vs-rest'
    # ROC AUC scores for classes 0, 1, and 2.
    out_0 = roc_auc_score([1, 0, 0, 0], y_scores[:, 0])
    out_1 = roc_auc_score([0, 1, 0, 0], y_scores[:, 1])
    out_2 = roc_auc_score([0, 0, 1, 1], y_scores[:, 2])
    result_unweighted = (out_0 + out_1 + out_2) / 3.

    assert_almost_equal(
        roc_auc_score(y_true, y_scores, multiclass="ovr"),
        result_unweighted)

    # Tests the weighted, one-vs-rest multiclass ROC AUC algorithm
    # on the same input (Provost & Domingos, 2001)
    result_weighted = out_0 * 0.25 + out_1 * 0.25 + out_2 * 0.5
    assert_almost_equal(
        roc_auc_score(y_true, y_scores, multiclass="ovr", average="weighted"),
        result_weighted)


def test_multi_auc_score_under_permutation():
    y_score = np.random.rand(100, 3)
    # Normalize the scores for each row
    row_sums = y_score.sum(axis=1)
    y_score = y_score / row_sums[:, np.newaxis]
    # Generate the true labels
    y_true = np.argmax(y_score, axis=1)
    y_true[np.random.randint(len(y_score), size=20)] = np.random.randint(
        2, size=20)
    for multiclass in ['ovr', 'ovo']:
        for average in ['macro', 'weighted']:
            same_score_under_permutation = None
            for perm in [[0, 1, 2], [0, 2, 1], [1, 0, 2],
                         [1, 2, 0], [2, 0, 1], [2, 1, 0]]:
                inv_perm = np.zeros(3, dtype=int)
                inv_perm[perm] = np.arange(3)
                y_score_perm = y_score[:, inv_perm]
                y_true_perm = np.take(perm, y_true)
                score = roc_auc_score(y_true_perm, y_score_perm,
                                      multiclass=multiclass, average=average)
                if same_score_under_permutation is None:
                    same_score_under_permutation = score
                else:
                    assert_almost_equal(score, same_score_under_permutation)


def test_auc_score_multi_error():
    # Test that roc_auc_score function returns an error when trying
    # to compute multiclass AUC for parameters where an output
    # is not defined.
    rng = check_random_state(404)
    y_pred = rng.rand(10, 3)
    row_sums = y_pred.sum(axis=1)
    y_pred = y_pred / row_sums[:, np.newaxis]
    y_true = rng.randint(0, 3, size=10)
    average_error_msg = ("Parameter 'average' must be one of "
                         "('macro', 'weighted') for multiclass problems.")
    assert_raise_message(ValueError, average_error_msg,
                         roc_auc_score, y_true, y_pred, average="samples")
    assert_raise_message(ValueError, average_error_msg,
                         roc_auc_score, y_true, y_pred, average="micro")
    multiclass_error_msg = ("Parameter multiclass='invalid' is not "
                            "supported for multiclass ROC AUC. 'multiclass' "
                            "must be one of ('ovo', 'ovr').")
    assert_raise_message(ValueError, multiclass_error_msg,
                         roc_auc_score, y_true, y_pred, multiclass="invalid")
    sample_weight_error_msg = ("Parameter 'sample_weight' is not supported "
                               "for multiclass one-vs-one ROC AUC. "
                               "'sample_weight' must be None in this case.")
    assert_raise_message(ValueError, sample_weight_error_msg,
                         roc_auc_score, y_true, y_pred,
                         multiclass="ovo", sample_weight=[])
    partial_comp_error_msg = ("Partial AUC computation not available in "
                              "multiclass setting. Parameter 'max_fpr' must "
                              "be set to `None`. Received `max_fpr=0.5` "
                              "instead.")
    assert_raise_message(ValueError, partial_comp_error_msg,
                         roc_auc_score, y_true, y_pred,
                         multiclass="ovo", max_fpr=0.5)


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
    y_true = np.full(10, -1, dtype="int")
    assert_raise_message(ValueError, "ROC AUC score is not defined",
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
        y_true = np.full(10, -1, dtype="int")
        assert_raise_message(ValueError, "ROC AUC score is not defined",
                             roc_auc_score, y_true, y_pred)


def test_binary_clf_curve():
    rng = check_random_state(404)
    y_true = rng.randint(0, 3, size=10)
    y_pred = rng.rand(10)
    msg = "multiclass format is not supported"
    assert_raise_message(ValueError, msg, precision_recall_curve,
                         y_true, y_pred)

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


def _test_precision_recall_curve(y_true, probas_pred):
    # Test Precision-Recall and aread under PR curve
    p, r, thresholds = precision_recall_curve(y_true, probas_pred)
    precision_recall_auc = _average_precision_slow(y_true, probas_pred)
    assert_array_almost_equal(precision_recall_auc, 0.859, 3)
    assert_array_almost_equal(precision_recall_auc,
                              average_precision_score(y_true, probas_pred))
    assert_almost_equal(_average_precision(y_true, probas_pred),
                        precision_recall_auc, decimal=3)
    assert_equal(p.size, r.size)
    assert_equal(p.size, thresholds.size + 1)
    # Smoke test in the case of proba having only one value
    p, r, thresholds = precision_recall_curve(y_true,
                                              np.zeros_like(probas_pred))
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
        # Here we are doing a terrible prediction: we are always getting
        # it wrong, hence the average_precision_score is the accuracy at
        # chance: 50%
        assert_almost_equal(auc_prc, 0.5)

        y_true = [1, 0]
        y_score = [1, 1]
        p, r, _ = precision_recall_curve(y_true, y_score)
        auc_prc = average_precision_score(y_true, y_score)
        assert_array_almost_equal(p, [0.5, 1])
        assert_array_almost_equal(r, [1., 0])
        assert_almost_equal(auc_prc, .5)

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
        assert_almost_equal(auc_prc, .5)

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
                            average="samples"), 0.75)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="micro"), 0.5)

        y_true = np.array([[1, 0], [0, 1]])
        y_score = np.array([[0, 1], [1, 0]])
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="macro"), 0.5)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="weighted"), 0.5)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="samples"), 0.5)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="micro"), 0.5)

        y_true = np.array([[1, 0], [0, 1]])
        y_score = np.array([[0.5, 0.5], [0.5, 0.5]])
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="macro"), 0.5)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="weighted"), 0.5)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="samples"), 0.5)
        assert_almost_equal(average_precision_score(y_true, y_score,
                            average="micro"), 0.5)


def test_average_precision_constant_values():
    # Check the average_precision_score of a constant predictor is
    # the TPR

    # Generate a dataset with 25% of positives
    y_true = np.zeros(100, dtype=int)
    y_true[::4] = 1
    # And a constant score
    y_score = np.ones(100)
    # The precision is then the fraction of positive whatever the recall
    # is, as there is only one threshold:
    assert_equal(average_precision_score(y_true, y_score), .25)


def test_average_precision_score_pos_label_errors():
    # Raise an error when pos_label is not in binary y_true
    y_true = np.array([0, 1])
    y_pred = np.array([0, 1])
    error_message = ("pos_label=2 is invalid. Set it to a label in y_true.")
    assert_raise_message(ValueError, error_message, average_precision_score,
                         y_true, y_pred, pos_label=2)
    # Raise an error for multilabel-indicator y_true with
    # pos_label other than 1
    y_true = np.array([[1, 0], [0, 1], [0, 1], [1, 0]])
    y_pred = np.array([[0.9, 0.1], [0.1, 0.9], [0.8, 0.2], [0.2, 0.8]])
    error_message = ("Parameter pos_label is fixed to 1 for multilabel"
                     "-indicator y_true. Do not set pos_label or set "
                     "pos_label to 1.")
    assert_raise_message(ValueError, error_message, average_precision_score,
                         y_true, y_pred, pos_label=0)


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


@pytest.mark.parametrize(
        'check',
        (check_lrap_toy,
         check_lrap_without_tie_and_increasing_score,
         check_lrap_only_ties,
         check_zero_or_all_relevant_labels))
@pytest.mark.parametrize(
        'func',
        (label_ranking_average_precision_score, _my_lrap))
def test_label_ranking_avp(check, func):
    check(func)


def test_lrap_error_raised():
    check_lrap_error_raised(label_ranking_average_precision_score)


@pytest.mark.parametrize('n_samples', (1, 2, 8, 20))
@pytest.mark.parametrize('n_classes', (2, 5, 10))
@pytest.mark.parametrize('random_state', range(1))
def test_alternative_lrap_implementation(n_samples, n_classes, random_state):

    check_alternative_lrap_implementation(
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


def test_partial_roc_auc_score():
    # Check `roc_auc_score` for max_fpr != `None`
    y_true = np.array([0, 0, 1, 1])
    assert roc_auc_score(y_true, y_true, max_fpr=1) == 1
    assert roc_auc_score(y_true, y_true, max_fpr=0.001) == 1
    with pytest.raises(ValueError):
        assert roc_auc_score(y_true, y_true, max_fpr=-0.1)
    with pytest.raises(ValueError):
        assert roc_auc_score(y_true, y_true, max_fpr=1.1)
    with pytest.raises(ValueError):
        assert roc_auc_score(y_true, y_true, max_fpr=0)

    y_scores = np.array([0.1,  0,  0.1, 0.01])
    roc_auc_with_max_fpr_one = roc_auc_score(y_true, y_scores, max_fpr=1)
    unconstrained_roc_auc = roc_auc_score(y_true, y_scores)
    assert roc_auc_with_max_fpr_one == unconstrained_roc_auc
    assert roc_auc_score(y_true, y_scores, max_fpr=0.3) == 0.5

    y_true, y_pred, _ = make_prediction(binary=True)
    for max_fpr in np.linspace(1e-4, 1, 5):
        assert_almost_equal(
            roc_auc_score(y_true, y_pred, max_fpr=max_fpr),
            _partial_roc_auc_score(y_true, y_pred, max_fpr))
