from __future__ import division

import random
import warnings
import numpy as np

from sklearn import datasets
from sklearn import svm

from sklearn.utils.testing import (assert_true,
                                   assert_raises,
                                   assert_equal,
                                   assert_almost_equal,
                                   assert_not_equal,
                                   assert_array_equal,
                                   assert_array_almost_equal)


from sklearn.metrics import (accuracy_score,
                             average_precision_score,
                             auc,
                             auc_score,
                             classification_report,
                             confusion_matrix,
                             explained_variance_score,
                             f1_score,
                             hinge_loss,
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
    p = range(n_samples)

    random.seed(0)
    random.shuffle(p)
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
    assert_array_almost_equal(roc_auc, 0.80, decimal=2)
    assert_almost_equal(roc_auc, auc_score(y_true, probas_pred))


def test_roc_curve_end_points():
    # Make sure that roc_curve returns a curve start at 0 and ending and
    # 1 even in corner cases
    rng = np.random.RandomState(0)
    y_true = np.array([0] * 50 + [1] * 50)
    y_pred = rng.randint(3, size=100)
    fpr, tpr, thr = roc_curve(y_true, y_pred)
    assert_equal(fpr[0], 0)
    assert_equal(fpr[-1], 1)


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


def test_roc_curve_multi():
    """roc_curve not applicable for multi-class problems"""
    y_true, _, probas_pred = make_prediction(binary=False)

    assert_raises(ValueError, roc_curve, y_true, probas_pred)


def test_roc_curve_confidence():
    """roc_curve for confidence scores"""
    y_true, _, probas_pred = make_prediction(binary=True)

    fpr, tpr, thresholds = roc_curve(y_true, probas_pred - 0.5)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.80, decimal=2)


def test_roc_curve_hard():
    """roc_curve for hard decisions"""
    y_true, pred, probas_pred = make_prediction(binary=True)

    # always predict one
    trivial_pred = np.ones(y_true.shape)
    fpr, tpr, thresholds = roc_curve(y_true, trivial_pred)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.50, decimal=2)

    # always predict zero
    trivial_pred = np.zeros(y_true.shape)
    fpr, tpr, thresholds = roc_curve(y_true, trivial_pred)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.50, decimal=2)

    # hard decisions
    fpr, tpr, thresholds = roc_curve(y_true, pred)
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.74, decimal=2)


def test_auc():
    """Test Area Under Curve (AUC) computation"""
    x = [0, 1]
    y = [0, 1]
    assert_array_almost_equal(auc(x, y), 0.5)
    x = [1, 0]
    y = [0, 1]
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
    x = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 1.]
    y = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
         1., 1., 1., 1., 1., 1., 1., 1.]
    assert_array_almost_equal(auc(x, y), 1.)


def test_auc_errors():
    # Incompatible shapes
    assert_raises(ValueError, auc, [0.0, 0.5, 1.0], [0.1, 0.2])

    # Too few x values
    assert_raises(ValueError, auc, [0.0], [0.1])


def test_precision_recall_f1_score_binary():
    """Test Precision Recall and F1 Score for binary classification task"""
    y_true, y_pred, _ = make_prediction(binary=True)

    # detailed measures for each class
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred, average=None)
    assert_array_almost_equal(p, [0.73, 0.75], 2)
    assert_array_almost_equal(r, [0.76, 0.72], 2)
    assert_array_almost_equal(f, [0.75, 0.74], 2)
    assert_array_equal(s, [25, 25])

    # individual scoring function that can be used for grid search: in the
    # binary class case the score is the value of the measure for the positive
    # class (e.g. label == 1)
    ps = precision_score(y_true, y_pred)
    assert_array_almost_equal(ps, 0.75, 2)

    rs = recall_score(y_true, y_pred)
    assert_array_almost_equal(rs, 0.72, 2)

    fs = f1_score(y_true, y_pred)
    assert_array_almost_equal(fs, 0.74, 2)


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

    cm = confusion_matrix(y_true, y_pred)
    assert_array_equal(cm, [[19, 6], [7, 18]])

    tp = cm[0, 0]
    tn = cm[1, 1]
    fp = cm[0, 1]
    fn = cm[1, 0]
    num = (tp * tn - fp * fn)
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if den == 0.:
        true_mcc = 0
    else:
        true_mcc = num / den
    mcc = matthews_corrcoef(y_true, y_pred)
    assert_array_almost_equal(mcc, true_mcc, decimal=2)
    assert_array_almost_equal(mcc, 0.48, decimal=2)


def test_matthews_corrcoef_nan():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert_equal(matthews_corrcoef([0], [1]), 0.0)


def test_precision_recall_f1_score_multiclass():
    """Test Precision Recall and F1 Score for multiclass classification task"""
    y_true, y_pred, _ = make_prediction(binary=False)

    # compute scores with default labels introspection
    p, r, f, s = precision_recall_fscore_support(y_true, y_pred, average=None)
    assert_array_almost_equal(p, [0.82, 0.55, 0.47], 2)
    assert_array_almost_equal(r, [0.92, 0.17, 0.90], 2)
    assert_array_almost_equal(f, [0.87, 0.26, 0.62], 2)
    assert_array_equal(s, [25, 30, 20])

    # averaging tests
    ps = precision_score(y_true, y_pred, pos_label=1, average='micro')
    assert_array_almost_equal(ps, 0.61, 2)

    rs = recall_score(y_true, y_pred, average='micro')
    assert_array_almost_equal(rs, 0.61, 2)

    fs = f1_score(y_true, y_pred, average='micro')
    assert_array_almost_equal(fs, 0.61, 2)

    ps = precision_score(y_true, y_pred, average='macro')
    assert_array_almost_equal(ps, 0.62, 2)

    rs = recall_score(y_true, y_pred, average='macro')
    assert_array_almost_equal(rs, 0.66, 2)

    fs = f1_score(y_true, y_pred, average='macro')
    assert_array_almost_equal(fs, 0.58, 2)

    ps = precision_score(y_true, y_pred, average='weighted')
    assert_array_almost_equal(ps, 0.62, 2)

    rs = recall_score(y_true, y_pred, average='weighted')
    assert_array_almost_equal(rs, 0.61, 2)

    fs = f1_score(y_true, y_pred, average='weighted')
    assert_array_almost_equal(fs, 0.55, 2)

    # same prediction but with and explicit label ordering
    p, r, f, s = precision_recall_fscore_support(
        y_true, y_pred, labels=[0, 2, 1], average=None)
    assert_array_almost_equal(p, [0.82, 0.47, 0.55], 2)
    assert_array_almost_equal(r, [0.92, 0.90, 0.17], 2)
    assert_array_almost_equal(f, [0.87, 0.62, 0.26], 2)
    assert_array_equal(s, [25, 20, 30])


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

    # compute confusion matrix with default labels introspection
    cm = confusion_matrix(y_true, y_pred)
    assert_array_equal(cm, [[23, 2, 0],
                            [5, 5, 20],
                            [0, 2, 18]])

    # compute confusion matrix with explicit label ordering
    cm = confusion_matrix(y_true, y_pred, labels=[0, 2, 1])
    assert_array_equal(cm, [[23, 0, 2],
                            [0, 18, 2],
                            [5, 20, 5]])


def test_confusion_matrix_multiclass_subset_labels():
    """Test confusion matrix - multi-class case with subset of labels"""
    y_true, y_pred, _ = make_prediction(binary=False)

    # compute confusion matrix with only first two labels considered
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    assert_array_equal(cm, [[23, 2],
                            [5, 5]])

    # compute confusion matrix with explicit label ordering for only subset
    # of labels
    cm = confusion_matrix(y_true, y_pred, labels=[2, 1])
    assert_array_equal(cm, [[18, 2],
                            [20, 5]])


def test_classification_report():
    """Test performance report"""
    iris = datasets.load_iris()
    y_true, y_pred, _ = make_prediction(dataset=iris, binary=False)

    # print classification report with class names
    expected_report = """\
             precision    recall  f1-score   support

     setosa       0.82      0.92      0.87        25
 versicolor       0.56      0.17      0.26        30
  virginica       0.47      0.90      0.62        20

avg / total       0.62      0.61      0.56        75
"""
    report = classification_report(
        y_true, y_pred, labels=range(len(iris.target_names)),
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
    report = classification_report(y_true, y_pred)
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


def _test_precision_recall_curve(y_true, probas_pred):
    """Test Precision-Recall and aread under PR curve"""
    p, r, thresholds = precision_recall_curve(y_true, probas_pred)
    precision_recall_auc = auc(r, p)
    assert_array_almost_equal(precision_recall_auc, 0.82, 2)
    assert_array_almost_equal(precision_recall_auc,
                              average_precision_score(y_true, probas_pred))
    # Smoke test in the case of proba having only one value
    p, r, thresholds = precision_recall_curve(y_true,
                                              np.zeros_like(probas_pred))
    precision_recall_auc = auc(r, p)
    assert_array_almost_equal(precision_recall_auc, 0.75, 3)


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

    # Classification
    # --------------
    with warnings.catch_warnings(True):
    # Throw deprecated warning
        assert_equal(zero_one(y_true, y_pred), 13)
        assert_almost_equal(zero_one(y_true, y_pred, normalize=True),
                            13 / float(n_samples), 2)

    assert_almost_equal(zero_one_loss(y_true, y_pred),
                        13 / float(n_samples), 2)
    assert_equal(zero_one_loss(y_true, y_pred, normalize=False), 13)
    assert_almost_equal(zero_one_loss(y_true, y_true), 0.0, 2)
    assert_almost_equal(zero_one_loss(y_true, y_true, normalize=False), 0, 2)

    assert_equal(accuracy_score(y_true, y_pred),
                 1 - zero_one_loss(y_true, y_pred))

    with warnings.catch_warnings(True):
    # Throw deprecated warning
        assert_equal(zero_one_score(y_true, y_pred),
                     1 - zero_one_loss(y_true, y_pred))

    # Regression
    # ----------
    assert_almost_equal(mean_squared_error(y_true, y_pred),
                        12.999 / n_samples, 2)
    assert_almost_equal(mean_squared_error(y_true, y_true),
                        0.00, 2)

    # mean_absolute_error and mean_squared_error are equal because
    # it is a binary problem.
    assert_almost_equal(mean_absolute_error(y_true, y_pred),
                        12.999 / n_samples, 2)
    assert_almost_equal(mean_absolute_error(y_true, y_true), 0.00, 2)

    assert_almost_equal(explained_variance_score(y_true, y_pred), -0.04, 2)
    assert_almost_equal(explained_variance_score(y_true, y_true), 1.00, 2)
    assert_equal(explained_variance_score([0, 0, 0], [0, 1, 1]), 0.0)

    assert_almost_equal(r2_score(y_true, y_pred), -0.04, 2)
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

    # symmetric
    assert_equal(accuracy_score(y_true, y_pred),
                 accuracy_score(y_pred, y_true))

    with warnings.catch_warnings(True):
        # Throw deprecated warning
        assert_equal(zero_one(y_true, y_pred),
                     zero_one(y_pred, y_true))

        assert_almost_equal(zero_one(y_true, y_pred, normalize=False),
                            zero_one(y_pred, y_true, normalize=False), 2)

    assert_equal(zero_one_loss(y_true, y_pred),
                 zero_one_loss(y_pred, y_true))

    assert_equal(zero_one_loss(y_true, y_pred, normalize=False),
                 zero_one_loss(y_pred, y_true, normalize=False))

    with warnings.catch_warnings(True):
    # Throw deprecated warning
        assert_equal(zero_one_score(y_true, y_pred),
                     zero_one_score(y_pred, y_true))

    assert_almost_equal(mean_squared_error(y_true, y_pred),
                        mean_squared_error(y_pred, y_true))

    assert_almost_equal(mean_absolute_error(y_true, y_pred),
                        mean_absolute_error(y_pred, y_true))

    # not symmetric
    assert_true(explained_variance_score(y_true, y_pred) !=
                explained_variance_score(y_pred, y_true))
    assert_true(r2_score(y_true, y_pred) !=
                r2_score(y_pred, y_true))
    # FIXME: precision and recall aren't symmetric either


def test_hinge_loss_binary():
    y_true = np.array([-1, 1, 1, -1])
    pred_decision = np.array([-8.5, 0.5, 1.5, -0.3])
    assert_equal(1.2 / 4, hinge_loss(y_true, pred_decision))

    y_true = np.array([0, 2, 2, 0])
    pred_decision = np.array([-8.5, 0.5, 1.5, -0.3])
    assert_equal(1.2 / 4,
                 hinge_loss(y_true, pred_decision, pos_label=2, neg_label=0))


def test_roc_curve_one_label():
    y_true = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    y_pred = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    # assert there are warnings
    with warnings.catch_warnings(True) as w:
        fpr, tpr, thresholds = roc_curve(y_true, y_pred)
        assert_equal(len(w), 1)
    # all true labels, all fpr should be nan
    assert_array_equal(fpr,
                       np.nan * np.ones(len(thresholds) + 1))
    # assert there are warnings
    with warnings.catch_warnings(True) as w:
        fpr, tpr, thresholds = roc_curve([1 - x for x in y_true],
                                         y_pred)
        assert_equal(len(w), 1)
    # all negative labels, all tpr should be nan
    assert_array_equal(tpr,
                       np.nan * np.ones(len(thresholds) + 1))


def test_multioutput_regression():
    y_true = np.array([[1, 0, 0, 1],
                       [0, 1, 1, 1],
                       [1, 1, 0, 1],
                       ])

    y_pred = np.array([[0, 0, 0, 1],
                       [1, 0, 1, 1],
                       [0, 0, 0, 1],
                       ])

    error = mean_squared_error(y_true, y_pred)
    assert_almost_equal(error, (1. / 3 + 2. / 3 + 2. / 3) / 4.)

    # mean_absolute_error and mean_squared_error are equal because
    # it is a binary problem.
    error = mean_absolute_error(y_true, y_pred)
    assert_almost_equal(error, (1. / 3 + 2. / 3 + 2. / 3) / 4.)

    error = r2_score(y_true, y_pred)
    assert_almost_equal(error, 1 - 5. / 2)


def test_multioutput_number_of_output_differ():
    y_true = np.array([[1, 0, 0, 1],
                       [0, 1, 1, 1],
                       [1, 1, 0, 1],
                       ])

    y_pred = np.array([[0, 0],
                       [1, 0],
                       [0, 0],
                       ])

    assert_raises(ValueError, mean_squared_error, y_true, y_pred)
    assert_raises(ValueError, mean_absolute_error, y_true, y_pred)
    assert_raises(ValueError, r2_score, y_true, y_pred)


def test_multioutput_regression_invariance_to_dimension_shuffling():
    # test invariance to dimension shuffling
    y_true, y_pred, _ = make_prediction()
    n_dims = 3
    y_true = np.reshape(y_true, (-1, n_dims))
    y_pred = np.reshape(y_pred, (-1, n_dims))

    for metric in [r2_score, mean_squared_error, mean_absolute_error]:
        error = metric(y_true, y_pred)

        for _ in xrange(3):
            perm = np.random.permutation(n_dims)
            assert_almost_equal(error,
                                metric(y_true[:, perm], y_pred[:, perm]))
