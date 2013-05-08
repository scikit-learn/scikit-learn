from __future__ import division

import warnings
import numpy as np

from sklearn import datasets
from sklearn import svm

from sklearn.preprocessing import LabelBinarizer
from sklearn.datasets import make_multilabel_classification
from sklearn.utils import check_random_state, shuffle
from sklearn.utils.multiclass import unique_labels
from sklearn.utils.testing import (assert_true,
                                   assert_raises,
                                   assert_raise_message,
                                   assert_equal,
                                   assert_almost_equal,
                                   assert_not_equal,
                                   assert_array_equal,
                                   assert_array_almost_equal,
                                   assert_greater)


from sklearn.metrics import (accuracy_score,
                             average_precision_score,
                             auc,
                             auc_score,
                             classification_report,
                             confusion_matrix,
                             explained_variance_score,
                             f1_score,
                             fbeta_score,
                             hamming_loss,
                             hinge_loss,
                             jaccard_similarity_score,
                             matthews_corrcoef,
                             mean_squared_error,
                             mean_absolute_error,
                             precision_recall_curve,
                             precision_recall_fscore_support,
                             precision_score,
                             recall_score,
                             r2_score,
                             roc_curve,
                             zero_one,
                             zero_one_score,
                             zero_one_loss)
from sklearn.externals.six.moves import xrange

ALL_METRICS = {
    "accuracy_score": accuracy_score,
    "unormalized_accuracy_score":
    lambda y1, y2: accuracy_score(y1, y2, normalize=False),

    "hamming_loss": hamming_loss,

    "jaccard_similarity_score": jaccard_similarity_score,
    "unormalized_jaccard_similarity_score":
    lambda y1, y2: jaccard_similarity_score(y1, y2, normalize=False),

    "zero_one_loss": zero_one_loss,
    "unnormalized_zero_one_loss":
    lambda y1, y2: zero_one_loss(y1, y2, normalize=False),

    "precision_score": precision_score,
    "recall_score": recall_score,
    "f1_score": f1_score,
    "f2_score": lambda y1, y2: fbeta_score(y1, y2, beta=2),
    "f0.5_score": lambda y1, y2: fbeta_score(y1, y2, beta=0.5),
    "matthews_corrcoef_score": matthews_corrcoef,
    "auc_score": auc_score,
    "average_precision_score": average_precision_score,

    "weighted_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="weighted", beta=0.5),
    "weighted_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="weighted"),
    "weighted_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="weighted", beta=2),
    "weighted_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="weighted"),
    "weighted_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="weighted"),

    "micro_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="micro", beta=0.5),
    "micro_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="micro"),
    "micro_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="micro", beta=2),
    "micro_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="micro"),
    "micro_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="micro"),

    "macro_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="macro", beta=0.5),
    "macro_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="macro"),
    "macro_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="macro", beta=2),
    "macro_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="macro"),
    "macro_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="macro"),

    "mean_absolute_error": mean_absolute_error,
    "mean_squared_error": mean_squared_error,
    "explained_variance_score": explained_variance_score,
    "r2_score": r2_score
}

METRICS_WITH_NORMALIZE_OPTION = {
    "accuracy_score ": lambda y1, y2, normalize:
    accuracy_score(y1, y2, normalize=normalize),
    "jaccard_similarity_score": lambda y1, y2, normalize:
    jaccard_similarity_score(y1, y2, normalize=normalize),
    "zero_one_loss": lambda y1, y2, normalize:
    zero_one_loss(y1, y2, normalize=normalize),
}

MULTILABELS_METRICS = {
    "accuracy_score": accuracy_score,
    "unormalized_accuracy_score":
    lambda y1, y2: accuracy_score(y1, y2, normalize=False),

    "hamming_loss": hamming_loss,

    "jaccard_similarity_score": jaccard_similarity_score,
    "unormalized_jaccard_similarity_score":
    lambda y1, y2: jaccard_similarity_score(y1, y2, normalize=False),

    "zero_one_loss": zero_one_loss,
    "unnormalized_zero_one_loss":
    lambda y1, y2: zero_one_loss(y1, y2, normalize=False),

    "weighted_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="weighted", beta=0.5),
    "weighted_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="weighted"),
    "weighted_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="weighted", beta=2),
    "weighted_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="weighted"),
    "weighted_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="weighted"),

    "micro_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="micro", beta=0.5),
    "micro_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="micro"),
    "micro_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="micro", beta=2),
    "micro_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="micro"),
    "micro_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="micro"),

    "macro_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="macro", beta=0.5),
    "macro_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="macro"),
    "macro_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="macro", beta=2),
    "macro_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="macro"),
    "macro_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="macro"),
}

MULTILABELS_METRICS_WITH_POS_LABELS = {
    "jaccard_similarity_score": jaccard_similarity_score,
    "unormalized_jaccard_similarity_score": lambda y1, y2, pos_label=1:
    jaccard_similarity_score(y1, y2, pos_label=pos_label, normalize=False),

    "weighted_f0.5_score": lambda y1, y2, pos_label=1:
    fbeta_score(y1, y2, pos_label=pos_label, average="weighted", beta=0.5),
    "weighted_f1_score": lambda y1, y2, pos_label=1:
    f1_score(y1, y2, pos_label=pos_label, average="weighted"),
    "weighted_f2_score": lambda y1, y2, pos_label=1:
    fbeta_score(y1, y2, pos_label=pos_label, average="weighted", beta=2),
    "weighted_precision_score": lambda y1, y2, pos_label=1:
    precision_score(y1, y2, pos_label=pos_label, average="weighted"),
    "weighted_recall_score": lambda y1, y2, pos_label=1:
    recall_score(y1, y2, pos_label=pos_label, average="weighted"),

    "micro_f0.5_score": lambda y1, y2, pos_label=1:
    fbeta_score(y1, y2, pos_label=pos_label, average="micro", beta=0.5),
    "micro_f1_score": lambda y1, y2, pos_label=1:
    f1_score(y1, y2, pos_label=pos_label, average="micro"),
    "micro_f2_score": lambda y1, y2, pos_label=1:
    fbeta_score(y1, y2, pos_label=pos_label, average="micro", beta=2),
    "micro_precision_score": lambda y1, y2, pos_label=1:
    precision_score(y1, y2, pos_label=pos_label, average="micro"),
    "micro_recall_score": lambda y1, y2, pos_label=1:
    recall_score(y1, y2, pos_label=pos_label, average="micro"),

    "macro_f0.5_score": lambda y1, y2, pos_label=1:
    fbeta_score(y1, y2, pos_label=pos_label, average="macro", beta=0.5),
    "macro_f1_score": lambda y1, y2, pos_label=1:
    f1_score(y1, y2, pos_label=pos_label, average="macro"),
    "macro_f2_score": lambda y1, y2, pos_label=1:
    fbeta_score(y1, y2, pos_label=pos_label, average="macro", beta=2),
    "macro_precision_score": lambda y1, y2, pos_label=1:
    precision_score(y1, y2, pos_label=pos_label, average="macro"),
    "macro_recall_score": lambda y1, y2, pos_label=1:
    recall_score(y1, y2, pos_label=pos_label, average="macro"),
}

SYMETRIC_METRICS = {
    "accuracy_score": accuracy_score,
    "unormalized_accuracy_score":
    lambda y1, y2: accuracy_score(y1, y2, normalize=False),

    "hamming_loss": hamming_loss,

    "jaccard_similarity_score": jaccard_similarity_score,
    "unormalized_jaccard_similarity_score":
    lambda y1, y2: jaccard_similarity_score(y1, y2, normalize=False),

    "zero_one_loss": zero_one_loss,
    "unnormalized_zero_one_loss":
    lambda y1, y2: zero_one_loss(y1, y2, normalize=False),

    "f1_score": f1_score,
    "weighted_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="weighted"),
    "micro_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="micro"),
    "macro_f1_score":
    lambda y1, y2: f1_score(y1, y2, average="macro"),

    "matthews_corrcoef_score": matthews_corrcoef,
    "mean_absolute_error": mean_absolute_error,
    "mean_squared_error": mean_squared_error
}

NOT_SYMETRIC_METRICS = {
    "explained_variance_score": explained_variance_score,
    "r2_score": r2_score,

    "precision_score": precision_score,
    "recall_score": recall_score,
    "f2_score": lambda y1, y2: fbeta_score(y1, y2, beta=2),
    "f0.5_score": lambda y1, y2: fbeta_score(y1, y2, beta=0.5),

    "weighted_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="weighted", beta=0.5),
    "weighted_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="weighted", beta=2),
    "weighted_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="weighted"),
    "weighted_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="weighted"),

    "micro_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="micro", beta=0.5),
    "micro_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="micro", beta=2),
    "micro_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="micro"),
    "micro_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="micro"),

    "macro_f0.5_score":
    lambda y1, y2: fbeta_score(y1, y2, average="macro", beta=0.5),
    "macro_f2_score":
    lambda y1, y2: fbeta_score(y1, y2, average="macro", beta=2),
    "macro_precision_score":
    lambda y1, y2: precision_score(y1, y2, average="macro"),
    "macro_recall_score":
    lambda y1, y2: recall_score(y1, y2, average="macro"),
}

THRESHOLDED_METRICS = {
    "auc_score": auc_score,
    "average_precision_score": average_precision_score,
}


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
    clf = svm.SVC(kernel='linear', probability=True)
    probas_pred = clf.fit(X[:half], y[:half]).predict_proba(X[half:])

    if binary:
        # only interested in probabilities of the positive case
        # XXX: do we really want a special API for the binary case?
        probas_pred = probas_pred[:, 1]

    y_pred = clf.predict(X[half:])
    y_true = y[half:]
    return y_true, y_pred, probas_pred


def test_roc_curve():
    """Test Area under Receiver Operating Characteristic (ROC) curve"""
    y_true, _, probas_pred = make_prediction(binary=True)

    fpr, tpr, thresholds = roc_curve(y_true, probas_pred)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.90, decimal=2)
    assert_almost_equal(roc_auc, auc_score(y_true, probas_pred))
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)


def test_roc_curve_end_points():
    # Make sure that roc_curve returns a curve start at 0 and ending and
    # 1 even in corner cases
    rng = np.random.RandomState(0)
    y_true = np.array([0] * 50 + [1] * 50)
    y_pred = rng.randint(3, size=100)
    fpr, tpr, thr = roc_curve(y_true, y_pred)
    assert_equal(fpr[0], 0)
    assert_equal(fpr[-1], 1)
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thr.shape)


def test_roc_returns_consistency():
    """Test whether the returned threshold matches up with tpr"""
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
    """roc_curve not applicable for multi-class problems"""
    y_true, _, probas_pred = make_prediction(binary=False)

    assert_raises(ValueError, roc_curve, y_true, probas_pred)


def test_roc_curve_confidence():
    """roc_curve for confidence scores"""
    y_true, _, probas_pred = make_prediction(binary=True)

    fpr, tpr, thresholds = roc_curve(y_true, probas_pred - 0.5)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.90, decimal=2)
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)


def test_roc_curve_hard():
    """roc_curve for hard decisions"""
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
    with warnings.catch_warnings(record=True) as w:
        fpr, tpr, thresholds = roc_curve(y_true, y_pred)
        assert_equal(len(w), 1)
    # all true labels, all fpr should be nan
    assert_array_equal(fpr,
                       np.nan * np.ones(len(thresholds)))
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)

    # assert there are warnings
    with warnings.catch_warnings(record=True) as w:
        fpr, tpr, thresholds = roc_curve([1 - x for x in y_true],
                                         y_pred)
        assert_equal(len(w), 1)
    # all negative labels, all tpr should be nan
    assert_array_equal(tpr,
                       np.nan * np.ones(len(thresholds)))
    assert_equal(fpr.shape, tpr.shape)
    assert_equal(fpr.shape, thresholds.shape)


def test_auc():
    """Test Area Under Curve (AUC) computation"""
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
    """Test that auc_score function returns an error when trying to compute AUC
    for non-binary class values.
    """
    rng = check_random_state(404)
    y_pred = rng.rand(10)
    # y_true contains only one class value
    y_true = np.zeros(10, dtype="int")
    assert_raise_message(ValueError, "AUC is defined for binary "
                         "classification only", auc_score, y_true, y_pred)
    y_true = np.ones(10, dtype="int")
    assert_raise_message(ValueError, "AUC is defined for binary "
                         "classification only", auc_score, y_true, y_pred)
    y_true = -np.ones(10, dtype="int")
    assert_raise_message(ValueError, "AUC is defined for binary "
                         "classification only", auc_score, y_true, y_pred)
    # y_true contains three different class values
    y_true = rng.randint(0, 3, size=10)
    assert_raise_message(ValueError, "AUC is defined for binary "
                         "classification only", auc_score, y_true, y_pred)


def test_precision_recall_f1_score_binary():
    """Test Precision Recall and F1 Score for binary classification task"""
    y_true, y_pred, _ = make_prediction(binary=True)

    # detailed measures for each class
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred, average=None)
    assert_array_almost_equal(p, [0.73, 0.85], 2)
    assert_array_almost_equal(r, [0.88, 0.68], 2)
    assert_array_almost_equal(f, [0.80, 0.76], 2)
    assert_array_equal(s, [25, 25])

    # individual scoring function that can be used for grid search: in the
    # binary class case the score is the value of the measure for the positive
    # class (e.g. label == 1)
    ps = precision_score(y_true, y_pred)
    assert_array_almost_equal(ps, 0.85, 2)

    rs = recall_score(y_true, y_pred)
    assert_array_almost_equal(rs, 0.68, 2)

    fs = f1_score(y_true, y_pred)
    assert_array_almost_equal(fs, 0.76, 2)


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
    """Test confusion matrix - binary classification case"""
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
    test(map(str, y_true), map(str, y_pred))


def test_matthews_corrcoef_nan():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert_equal(matthews_corrcoef([0], [1]), 0.0)


def test_precision_recall_f1_score_multiclass():
    """Test Precision Recall and F1 Score for multiclass classification task"""
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

    # same prediction but with and explicit label ordering
    p, r, f, s = precision_recall_fscore_support(
        y_true, y_pred, labels=[0, 2, 1], average=None)
    assert_array_almost_equal(p, [0.83, 0.41, 0.33], 2)
    assert_array_almost_equal(r, [0.79, 0.90, 0.10], 2)
    assert_array_almost_equal(f, [0.81, 0.57, 0.15], 2)
    assert_array_equal(s, [24, 20, 31])


def test_precision_recall_f1_score_multiclass_pos_label_none():
    """Test Precision Recall and F1 Score for multiclass classification task

    GH Issue #1296
    """
    # initialize data
    y_true = np.array([0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1])
    y_pred = np.array([1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1])

    # compute scores with default labels introspection
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                 pos_label=None,
                                                 average='weighted')


def test_zero_precision_recall():
    """Check that pathological cases do not bring NaNs"""

    try:
        old_error_settings = np.seterr(all='raise')

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
    """Test confusion matrix - multi-class case"""
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
    test(map(str, y_true), map(str, y_pred), string_type=True)


def test_confusion_matrix_multiclass_subset_labels():
    """Test confusion matrix - multi-class case with subset of labels"""
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


def test_classification_report_binary_classification_with_pos_label():
    iris = datasets.load_iris()
    y_true, y_pred, _ = make_prediction(dataset=iris, binary=True)

    print y_true
    expected_report = """\
             precision    recall  f1-score   support

          0       0.73      0.88      0.80        25
          1       0.85      0.68      0.76        25

avg / total       0.79      0.78      0.78        50
"""
    for pos_label in [0, 1]:
        report = classification_report(y_true, y_pred, pos_label=pos_label)
        assert_equal(report, expected_report)


def test_classification_report_multiclass():
    """Test performance report"""
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

          0       0.82      0.92      0.87        25
          1       0.56      0.17      0.26        30
          2       0.47      0.90      0.62        20

avg / total       0.62      0.61      0.56        75
"""
    expected_report = """\
             precision    recall  f1-score   support

          0       0.83      0.79      0.81        24
          1       0.33      0.10      0.15        31
          2       0.42      0.90      0.57        20

avg / total       0.51      0.53      0.47        75
"""
    report = classification_report(y_true, y_pred)
    assert_equal(report, expected_report)


def test_multilabel_classification_report():

    n_classes = 4
    n_samples = 50
    _, y_true_ll = make_multilabel_classification(n_features=1,
                                                  n_classes=n_classes,
                                                  random_state=0,
                                                  n_samples=n_samples)
    _, y_pred_ll = make_multilabel_classification(n_features=1,
                                                  n_classes=n_classes,
                                                  random_state=1,
                                                  n_samples=n_samples)

    expected_report = """\
             precision    recall  f1-score   support

          0       0.39      0.73      0.51        15
          1       0.57      0.75      0.65        28
          2       0.33      0.11      0.17        18
          3       0.44      0.50      0.47        24

avg / total       0.45      0.54      0.47        85
"""

    lb = LabelBinarizer()
    lb.fit([range(4)])
    y_true_bi = lb.transform(y_true_ll)
    y_pred_bi = lb.transform(y_pred_ll)

    for y_true, y_pred in [(y_true_ll, y_pred_ll), (y_true_bi, y_pred_bi)]:
        report = classification_report(y_true, y_pred)
        assert_equal(report, expected_report)

    # With a given pos_label
    pos_label = 5
    y_true_bi = y_true_bi * pos_label
    y_pred_bi = y_pred_bi * pos_label

    expected_report = """\
             precision    recall  f1-score   support

          0       0.39      0.73      0.51        15
          1       0.57      0.75      0.65        28
          2       0.33      0.11      0.17        18
          3       0.44      0.50      0.47        24

avg / total       0.45      0.54      0.47        85
"""
    report = classification_report(y_true_bi, y_pred_bi, pos_label=pos_label)
    assert_equal(report, expected_report)


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
    """Test Precision-Recall and aread under PR curve"""
    p, r, thresholds = precision_recall_curve(y_true, probas_pred)
    precision_recall_auc = auc(r, p)
    assert_array_almost_equal(precision_recall_auc, 0.85, 2)
    assert_array_almost_equal(precision_recall_auc,
                              average_precision_score(y_true, probas_pred))
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


def test_score_scale_invariance():
    # Test that average_precision_score and auc_score are invariant by
    # the scaling or shifting of probabilities
    y_true, _, probas_pred = make_prediction(binary=True)

    roc_auc = auc_score(y_true, probas_pred)
    roc_auc_scaled = auc_score(y_true, 100 * probas_pred)
    roc_auc_shifted = auc_score(y_true, probas_pred - 10)
    assert_equal(roc_auc, roc_auc_scaled)
    assert_equal(roc_auc, roc_auc_shifted)

    pr_auc = average_precision_score(y_true, probas_pred)
    pr_auc_scaled = average_precision_score(y_true, 100 * probas_pred)
    pr_auc_shifted = average_precision_score(y_true, probas_pred - 10)
    assert_equal(pr_auc, pr_auc_scaled)
    assert_equal(pr_auc, pr_auc_shifted)


def test_losses():
    """Test loss functions"""
    y_true, y_pred, _ = make_prediction(binary=True)
    n_samples = y_true.shape[0]
    n_classes = np.size(unique_labels(y_true))

    # Classification
    # --------------
    with warnings.catch_warnings(record=True):
    # Throw deprecated warning
        assert_equal(zero_one(y_true, y_pred), 11)

    assert_almost_equal(zero_one_loss(y_true, y_pred),
                        11 / float(n_samples), 2)
    assert_equal(zero_one_loss(y_true, y_pred, normalize=False), 11)

    assert_almost_equal(zero_one_loss(y_true, y_true), 0.0, 2)
    assert_almost_equal(hamming_loss(y_true, y_pred),
                        2 * 11. / (n_samples * n_classes), 2)

    assert_equal(accuracy_score(y_true, y_pred),
                 1 - zero_one_loss(y_true, y_pred))

    with warnings.catch_warnings(True):
    # Throw deprecated warning
        assert_equal(zero_one_score(y_true, y_pred),
                     1 - zero_one_loss(y_true, y_pred))

    # Regression
    # ----------
    assert_almost_equal(mean_squared_error(y_true, y_pred),
                        10.999 / n_samples, 2)
    assert_almost_equal(mean_squared_error(y_true, y_true),
                        0.00, 2)

    # mean_absolute_error and mean_squared_error are equal because
    # it is a binary problem.
    assert_almost_equal(mean_absolute_error(y_true, y_pred),
                        10.999 / n_samples, 2)
    assert_almost_equal(mean_absolute_error(y_true, y_true), 0.00, 2)

    assert_almost_equal(explained_variance_score(y_true, y_pred), 0.16, 2)
    assert_almost_equal(explained_variance_score(y_true, y_true), 1.00, 2)
    assert_equal(explained_variance_score([0, 0, 0], [0, 1, 1]), 0.0)

    assert_almost_equal(r2_score(y_true, y_pred), 0.12, 2)
    assert_almost_equal(r2_score(y_true, y_true), 1.00, 2)
    assert_equal(r2_score([0, 0, 0], [0, 0, 0]), 1.0)
    assert_equal(r2_score([0, 0, 0], [0, 1, 1]), 0.0)


def test_losses_at_limits():
    # test limit cases
    assert_almost_equal(mean_squared_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(mean_absolute_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(explained_variance_score([0.], [0.]), 1.00, 2)
    assert_almost_equal(r2_score([0., 1], [0., 1]), 1.00, 2)


def test_r2_one_case_error():
    # test whether r2_score raises error given one point
    assert_raises(ValueError, r2_score, [0], [0])


def test_symmetry():
    """Test the symmetry of score and loss functions"""
    y_true, y_pred, _ = make_prediction(binary=True)

    # We shouldn't forget any metrics
    assert_equal(set(SYMETRIC_METRICS).union(NOT_SYMETRIC_METRICS,
                                             THRESHOLDED_METRICS),
                 set(ALL_METRICS))

    assert_equal(set(SYMETRIC_METRICS).intersection(set(NOT_SYMETRIC_METRICS)),
                 set([]))

    # Symmetric metric
    for name, metric in SYMETRIC_METRICS.items():
        assert_equal(metric(y_true, y_pred),
                     metric(y_pred, y_true),
                     msg="%s is not symetric" % name)

    # Not symmetric metrics
    for name, metric in NOT_SYMETRIC_METRICS.items():
        assert_true(metric(y_true, y_pred) != metric(y_pred, y_true),
                    msg="%s seems to be symetric" % name)

    # Deprecated metrics
    with warnings.catch_warnings(record=True):
        # Throw deprecated warning
        assert_equal(zero_one(y_true, y_pred),
                     zero_one(y_pred, y_true))

        assert_equal(zero_one(y_true, y_pred, normalize=False),
                     zero_one(y_pred, y_true, normalize=False))

        assert_equal(zero_one_score(y_true, y_pred),
                     zero_one_score(y_pred, y_true))


def test_sample_order_invariance():
    y_true, y_pred, _ = make_prediction(binary=True)

    y_true_shuffle, y_pred_shuffle = shuffle(y_true, y_pred,
                                             random_state=0)

    for name, metric in ALL_METRICS.items():

        assert_almost_equal(metric(y_true, y_pred),
                            metric(y_true_shuffle, y_pred_shuffle),
                            err_msg="%s is not sample order invariant"
                                    % metric)


def test_format_invariance_with_1d_vectors():
    y1, y2, _ = make_prediction(binary=True)

    y1_list = list(y1)
    y2_list = list(y2)

    y1_1d, y2_1d = np.array(y1), np.array(y2)
    assert_equal(y1_1d.ndim, 1)
    assert_equal(y2_1d.ndim, 1)
    y1_column = np.reshape(y1_1d, (-1, 1))
    y2_column = np.reshape(y2_1d, (-1, 1))
    y1_row = np.reshape(y1_1d, (1, -1))
    y2_row = np.reshape(y2_1d, (1, -1))

    for name, metric in ALL_METRICS.items():

        measure = metric(y1, y2)

        assert_almost_equal(metric(y1_list, y2_list), measure,
                            err_msg="%s is not representation invariant"
                                    "with list" % name)

        assert_almost_equal(metric(y1_1d, y2_1d), measure,
                            err_msg="%s is not representation invariant"
                                    "with np-array-1d" % name)

        assert_almost_equal(metric(y1_column, y2_column), measure,
                            err_msg="%s is not representation invariant "
                                    "with np-array-column" % name)

        assert_almost_equal(metric(y1_row, y2_row), measure,
                            err_msg="%s is not representation invariant "
                                    "with np-array-row" % name)

        # Mix format support
        assert_almost_equal(metric(y1_1d, y2_list), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix np-array-1d and list" % name)

        assert_almost_equal(metric(y1_list, y2_1d), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix np-array-1d and list" % name)

        assert_almost_equal(metric(y1_1d, y2_column), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix np-array-1d and np-array-column"
                                    % name)

        assert_almost_equal(metric(y1_column, y2_1d), measure,
                            err_msg="%s is not representation invariant "
                                    "with mix np-array-1d and np-array-column"
                                    % name)

        assert_almost_equal(metric(y1_list, y2_column), measure,
                            err_msg="%s is not representation invariant"
                                    "with mix list and np-array-column"
                                    % name)

        assert_almost_equal(metric(y1_column, y2_list), measure,
                            err_msg="%s is not representation invariant"
                                    "with mix list and np-array-column"
                                    % name)

        # At the moment, these mix representations aren't allowed
        assert_raises(ValueError, metric, y1_1d, y2_row)
        assert_raises(ValueError, metric, y1_row, y2_1d)
        assert_raises(ValueError, metric, y1_list, y2_row)
        assert_raises(ValueError, metric, y1_row, y2_list)
        assert_raises(ValueError, metric, y1_column, y2_row)
        assert_raises(ValueError, metric, y1_row, y2_column)


def test_hinge_loss_binary():
    y_true = np.array([-1, 1, 1, -1])
    pred_decision = np.array([-8.5, 0.5, 1.5, -0.3])
    assert_equal(hinge_loss(y_true, pred_decision), 1.2 / 4)

    with warnings.catch_warnings():
        # Test deprecated pos_label
        assert_equal(
            hinge_loss(-y_true, pred_decision),
            hinge_loss(y_true, pred_decision, pos_label=-1, neg_label=1))

    y_true = np.array([0, 2, 2, 0])
    pred_decision = np.array([-8.5, 0.5, 1.5, -0.3])

    assert_equal(hinge_loss(y_true, pred_decision), 1.2 / 4)
    with warnings.catch_warnings():
        # Test deprecated pos_label
        assert_equal(hinge_loss(y_true, pred_decision, pos_label=2, neg_label=0),
                     1.2 / 4)


def test_multioutput_regression():
    y_true = np.array([[1, 0, 0, 1], [0, 1, 1, 1], [1, 1, 0, 1]])
    y_pred = np.array([[0, 0, 0, 1], [1, 0, 1, 1], [0, 0, 0, 1]])

    error = mean_squared_error(y_true, y_pred)
    assert_almost_equal(error, (1. / 3 + 2. / 3 + 2. / 3) / 4.)

    # mean_absolute_error and mean_squared_error are equal because
    # it is a binary problem.
    error = mean_absolute_error(y_true, y_pred)
    assert_almost_equal(error, (1. / 3 + 2. / 3 + 2. / 3) / 4.)

    error = r2_score(y_true, y_pred)
    assert_almost_equal(error, 1 - 5. / 2)


def test_multioutput_number_of_output_differ():
    y_true = np.array([[1, 0, 0, 1], [0, 1, 1, 1], [1, 1, 0, 1]])
    y_pred = np.array([[0, 0], [1, 0], [0, 0]])

    assert_raises(ValueError, mean_squared_error, y_true, y_pred)
    assert_raises(ValueError, mean_absolute_error, y_true, y_pred)
    assert_raises(ValueError, r2_score, y_true, y_pred)


def test_multioutput_regression_invariance_to_dimension_shuffling():
    # test invariance to dimension shuffling
    y_true, y_pred, _ = make_prediction()
    n_dims = 3
    y_true = np.reshape(y_true, (-1, n_dims))
    y_pred = np.reshape(y_pred, (-1, n_dims))

    rng = check_random_state(314159)
    for metric in [r2_score, mean_squared_error, mean_absolute_error]:
        error = metric(y_true, y_pred)

        for _ in xrange(3):
            perm = rng.permutation(n_dims)
            assert_almost_equal(metric(y_true[:, perm], y_pred[:, perm]),
                                error)


def test_multilabel_representation_invariance():

    # Generate some data
    n_classes = 4
    n_samples = 50
    _, y1 = make_multilabel_classification(n_features=1, n_classes=n_classes,
                                           random_state=0, n_samples=n_samples)
    _, y2 = make_multilabel_classification(n_features=1, n_classes=n_classes,
                                           random_state=1, n_samples=n_samples)

    # Be sure to have at least one empty label
    y1 += ([], )
    y2 += ([], )

    # NOTE: The "sorted" trick is necessary to shuffle labels, because it
    # allows to return the shuffled tuple.
    rng = check_random_state(42)
    shuffled = lambda x: sorted(x, key=lambda *args: rng.rand())
    y1_shuffle = [shuffled(x) for x in y1]
    y2_shuffle = [shuffled(x) for x in y2]

    # Let's have redundant labels
    y1_redundant = [x * rng.randint(1, 4) for x in y1]
    y2_redundant = [x * rng.randint(1, 4) for x in y2]

    # Binary indicator matrix format
    lb = LabelBinarizer().fit([range(n_classes)])
    y1_binary_indicator = lb.transform(y1)
    y2_binary_indicator = lb.transform(y2)

    y1_shuffle_binary_indicator = lb.transform(y1_shuffle)
    y2_shuffle_binary_indicator = lb.transform(y2_shuffle)

    for name, metric in MULTILABELS_METRICS.items():
        measure = metric(y1, y2)

        # Check representation invariance
        assert_almost_equal(metric(y1_binary_indicator, y2_binary_indicator),
                            measure,
                            err_msg="%s failed representation invariance  "
                                    "between list of list of labels format "
                                    "and dense binary indicator format."
                                    % name)

        # Check invariance with redundant labels with list of labels
        assert_almost_equal(metric(y1, y2_redundant), measure,
                            err_msg="%s failed rendundant label invariance"
                                    % name)

        assert_almost_equal(metric(y1_redundant, y2_redundant), measure,
                            err_msg="%s failed rendundant label invariance"
                                    % name)

        assert_almost_equal(metric(y1_redundant, y2), measure,
                            err_msg="%s failed rendundant label invariance"
                                    % name)

        # Check shuffling invariance with list of labels
        assert_almost_equal(metric(y1_shuffle, y2_shuffle), measure,
                            err_msg="%s failed shuffling invariance "
                                    "with list of list of labels format."
                                    % name)

        # Check shuffling invariance with dense binary indicator matrix
        assert_almost_equal(metric(y1_shuffle_binary_indicator,
                                   y2_shuffle_binary_indicator), measure,
                            err_msg="%s failed shuffling invariance "
                                    " with dense binary indicator format."
                                    % name)

        # Check invariance with mix input representation
        assert_almost_equal(metric(y1, y2_binary_indicator), measure,
                            err_msg="%s failed mix input representation"
                                    "invariance: y_true in list of list of "
                                    "labels format and y_pred in dense binary"
                                    "indicator format"
                                    % name)

        assert_almost_equal(metric(y1_binary_indicator, y2), measure,
                            err_msg="%s failed mix input representation"
                                    "invariance: y_true in dense binary "
                                    "indicator format and y_pred in list of "
                                    "list of labels format."
                                    % name)


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

    # With a given pos_label
    assert_equal(jaccard_similarity_score(y1, y2, pos_label=0), 0.75)
    assert_equal(jaccard_similarity_score(y2, np.zeros(y1.shape),
                                          pos_label=0), 0.5)
    assert_equal(jaccard_similarity_score(y1, y2, pos_label=10), 1)

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


def test_normalize_option_binary_classification():
    # Test in the binary case
    y_true, y_pred, _ = make_prediction(binary=True)
    n_samples = y_true.shape[0]

    for name, metrics in METRICS_WITH_NORMALIZE_OPTION.items():
        measure = metrics(y_true, y_pred, normalize=True)
        assert_greater(measure, 0,
                       msg="We failed to test correctly the normalize option")
        assert_almost_equal(metrics(y_true, y_pred, normalize=False)
                            / n_samples, measure)


def test_normalize_option_multiclasss_classification():
    # Test in the multiclass case
    y_true, y_pred, _ = make_prediction(binary=False)
    n_samples = y_true.shape[0]

    for name, metrics in METRICS_WITH_NORMALIZE_OPTION.items():
        measure = metrics(y_true, y_pred, normalize=True)
        assert_greater(measure, 0,
                       msg="We failed to test correctly the normalize option")
        assert_almost_equal(metrics(y_true, y_pred, normalize=False)
                            / n_samples, measure)


def test_normalize_option_multilabel_classification():
    # Test in the multilabel case
    n_classes = 4
    n_samples = 100
    _, y_true = make_multilabel_classification(n_features=1,
                                               n_classes=n_classes,
                                               random_state=0,
                                               n_samples=n_samples)
    _, y_pred = make_multilabel_classification(n_features=1,
                                               n_classes=n_classes,
                                               random_state=1,
                                               n_samples=n_samples)

    # Be sure to have at least one empty label
    y_true += ([], )
    y_pred += ([], )
    n_samples += 1

    lb = LabelBinarizer().fit([range(n_classes)])
    y_true_binary_indicator = lb.transform(y_true)
    y_pred_binary_indicator = lb.transform(y_pred)

    for name, metrics in METRICS_WITH_NORMALIZE_OPTION.items():
        # List of list of labels
        measure = metrics(y_true, y_pred, normalize=True)
        assert_greater(measure, 0,
                       msg="We failed to test correctly the normalize option")
        assert_almost_equal(metrics(y_true, y_pred, normalize=False)
                            / n_samples, measure,
                            err_msg="Failed with %s" % name)

        # Indicator matrix format
        measure = metrics(y_true_binary_indicator,
                          y_pred_binary_indicator, normalize=True)
        assert_greater(measure, 0,
                       msg="We failed to test correctly the normalize option")
        assert_almost_equal(metrics(y_true_binary_indicator,
                                    y_pred_binary_indicator, normalize=False)
                            / n_samples, measure,
                            err_msg="Failed with %s" % name)


def test_precision_recall_f1_score_multilabel_1():
    """ Test precision_recall_f1_score on a crafted multilabel example
    """
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

        # Check macro
        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="macro")
        assert_almost_equal(p, 1.5 / 4)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 2.5 / 1.5 * 0.25)
        assert_equal(s, None)

        # Check micro
        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="micro")
        assert_almost_equal(p, 0.5)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 0.5)
        assert_equal(s, None)

        # Check weigted
        # |h(x_i) inter y_i | = [0, 1, 1]
        # |y_i| = [1, 1, 2]
        # |h(x_i)| = [1, 1, 2]
        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="weighted")
        assert_almost_equal(p, 0.5)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 0.5)
        assert_equal(s, None)


def test_precision_recall_f1_score_multilabel_2():
    """ Test precision_recall_f1_score on a crafted multilabel example 2
    """
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

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="micro")
        assert_almost_equal(p, 0.25)
        assert_almost_equal(r, 0.25)
        assert_almost_equal(f, 2 * 0.25 * 0.25 / 0.5)
        assert_equal(s, None)

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="macro")
        assert_almost_equal(p, 0.25)
        assert_almost_equal(r, 0.125)
        assert_almost_equal(f, 2 / 12)
        assert_equal(s, None)

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="weighted")
        # Check weigted
        # |h(x_i) inter y_i | = [0, 0, 1]
        # |y_i| = [1, 1, 2]
        # |h(x_i)| = [1, 1, 2]
        assert_almost_equal(p, 1 / 6)
        assert_almost_equal(r, 1 / 6)
        assert_almost_equal(f, 2 / 4 * 1 / 3)
        assert_equal(s, None)


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

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="macro")
        assert_almost_equal(p, 0.5)
        assert_almost_equal(r, 1.5 / 4)
        assert_almost_equal(f, 2.5 / (4 * 1.5))
        assert_equal(s, None)

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="micro")
        assert_almost_equal(p, 2 / 3)
        assert_almost_equal(r, 0.5)
        assert_almost_equal(f, 2 / 3 / (2 / 3 + 0.5))
        assert_equal(s, None)

        p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                     average="weighted")
        # Check weigted
        # |h(x_i) inter y_i | = [0, 0, 2]
        # |y_i| = [1, 1, 2]
        # |h(x_i)| = [0, 1, 2]
        assert_almost_equal(p, 1 / 3)
        assert_almost_equal(r, 2 / 3)
        assert_almost_equal(f, 1 / 3)
        assert_equal(s, None)


def test_precision_recall_f1_no_labels():
    y_true = np.zeros((20, 3))
    y_pred = np.zeros_like(y_true)

    p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                 average=None)
    #tp = [0, 0, 0]
    #fn = [0, 0, 0]
    #fp = [0, 0, 0]

    # Check per class
    assert_array_almost_equal(p, [0, 0, 0], 2)
    assert_array_almost_equal(r, [0, 0, 0], 2)
    assert_array_almost_equal(f, [0, 0, 0], 2)
    assert_array_almost_equal(s, [0, 0, 0], 2)

    # Check macro
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                 average="macro")
    assert_almost_equal(p, 0)
    assert_almost_equal(r, 0)
    assert_almost_equal(f, 0)
    assert_equal(s, None)

    # Check micro
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                 average="micro")
    assert_almost_equal(p, 0)
    assert_almost_equal(r, 0)
    assert_almost_equal(f, 0)
    assert_equal(s, None)

    # # Check weigted
    # |h(x_i) inter y_i | = [0, 0, 0]
    # |y_i| = [0, 0, 0]
    # |h(x_i)| = [1, 1, 2]
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred,
                                                 average="weighted")
    assert_almost_equal(p, 1)
    assert_almost_equal(r, 1)
    assert_almost_equal(f, 1)
    assert_equal(s, None)


def test_multilabel_invariance_with_pos_labels():
    n_classes = 4
    n_samples = 50
    _, y1 = make_multilabel_classification(n_features=1, n_classes=n_classes,
                                           random_state=0, n_samples=n_samples)
    _, y2 = make_multilabel_classification(n_features=1, n_classes=n_classes,
                                           random_state=1, n_samples=n_samples)

    lb = LabelBinarizer().fit([range(n_classes)])
    y1_binary_indicator = lb.transform(y1)
    y2_binary_indicator = lb.transform(y2)

    for name, metric in MULTILABELS_METRICS_WITH_POS_LABELS.items():
        measure = metric(y1, y2)

        for pos_label in [1, 3]:
            assert_almost_equal(metric(y1_binary_indicator * pos_label,
                                       y2_binary_indicator * pos_label,
                                       pos_label=pos_label),
                                measure,
                                err_msg="%s is not representation invariant"
                                        "with pos_label=%s"
                                        % (metric, pos_label))
