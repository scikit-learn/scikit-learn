import random
import numpy as np

from nose.tools import raises
from nose.tools import assert_true, assert_raises
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_equal, assert_almost_equal

from ... import datasets
from ... import svm
from ..metrics import auc
from ..metrics import classification_report
from ..metrics import confusion_matrix
from ..metrics import explained_variance_score
from ..metrics import r2_score
from ..metrics import f1_score
from ..metrics import matthews_corrcoef
from ..metrics import mean_squared_error
from ..metrics import precision_recall_curve
from ..metrics import precision_recall_fscore_support
from ..metrics import precision_score
from ..metrics import recall_score
from ..metrics import roc_curve
from ..metrics import auc_score
from ..metrics import average_precision_score
from ..metrics import zero_one
from ..metrics import hinge_loss


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


def test_roc_returns_consistency():
    """Test whether the returned threshold matches up with tpr"""
    # make small toy dataset
    y_true, _, probas_pred = make_prediction(binary=True)
    fpr, tpr, thresholds = roc_curve(y_true, probas_pred)

    # use the given thresholds to determine the tpr
    tpr_correct = []
    for t in range(len(thresholds)):
        tp = np.sum((probas_pred >= thresholds[t]) & y_true)
        p = np.sum(y_true)
        tpr_correct.append(1.0 * tp / p)

    # compare tpr and tpr_correct to see if the thresholds' order was correct
    assert_array_almost_equal(tpr, tpr_correct, decimal=2)


@raises(ValueError)
def test_roc_curve_multi():
    """roc_curve not applicable for multi-class problems"""
    y_true, _, probas_pred = make_prediction(binary=False)

    fpr, tpr, thresholds = roc_curve(y_true, probas_pred)


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
    """Test Area Under Curve (AUC) computation with duplicate values

    auc() was previously sorting the x and y arrays according to the indices
    from numpy.argsort(x), which was reordering the tied 0's in this example
    and resulting in an incorrect area computation. This test detects the
    error.
    """
    x = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 1.]
    y = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
         1., 1., 1., 1., 1., 1., 1., 1.]
    assert_array_almost_equal(auc(x, y), 1.)


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
    assert_array_equal(cm, [[23, 2,  0],
                            [5,  5, 20],
                            [0,  2, 18]])

    # compute confusion matrix with explicit label ordering
    cm = confusion_matrix(y_true, y_pred, labels=[0, 2, 1])
    assert_array_equal(cm, [[23, 0,  2],
                            [0, 18,  2],
                            [5, 20,  5]])


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
    """Test Precision-Recall and aread under PR curve"""
    y_true, _, probas_pred = make_prediction(binary=True)

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


def test_losses():
    """Test loss functions"""
    y_true, y_pred, _ = make_prediction(binary=True)
    n = y_true.shape[0]

    assert_equal(zero_one(y_true, y_pred), 13)
    assert_almost_equal(mean_squared_error(y_true, y_pred), 12.999 / n, 2)
    assert_almost_equal(mean_squared_error(y_true, y_true), 0.00, 2)

    assert_almost_equal(explained_variance_score(y_true, y_pred), -0.04, 2)
    assert_almost_equal(explained_variance_score(y_true, y_true), 1.00, 2)

    assert_almost_equal(r2_score(y_true, y_pred), -0.04, 2)
    assert_almost_equal(r2_score(y_true, y_true), 1.00, 2)


def test_losses_at_limits():
    # test limit cases
    assert_almost_equal(mean_squared_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(explained_variance_score([0.], [0.]), 1.00, 2)
    assert_almost_equal(r2_score([0., 1], [0., 1]), 1.00, 2)


def test_r2_one_case_error():
    # test whether r2_score raises error given one point
    assert_raises(ValueError, r2_score, [0], [0])


def test_symmetry():
    """Test the symmetry of score and loss functions"""
    y_true, y_pred, _ = make_prediction(binary=True)

    # symmetric
    assert_equal(zero_one(y_true, y_pred),
                 zero_one(y_pred, y_true))
    assert_almost_equal(mean_squared_error(y_true, y_pred),
                        mean_squared_error(y_pred, y_true))
    # not symmetric
    assert_true(explained_variance_score(y_true, y_pred) != \
            explained_variance_score(y_pred, y_true))
    assert_true(r2_score(y_true, y_pred) != \
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
