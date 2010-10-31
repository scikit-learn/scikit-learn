import random
import numpy as np
import nose

from numpy.testing import assert_
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_equal, assert_almost_equal

from .. import datasets
from .. import svm
from ..metrics import auc
from ..metrics import confusion_matrix
from ..metrics import explained_variance
from ..metrics import f1_score
from ..metrics import mean_square_error
from ..metrics import precision
from ..metrics import precision_recall
from ..metrics import precision_recall_curve
from ..metrics import recall
from ..metrics import roc_curve
from ..metrics import zero_one

# import some data to play with
iris = datasets.load_iris()
X = iris.data
y = iris.target

# restrict to a binary classification task
X, y = X[y != 2], y[y != 2]

n_samples, n_features = X.shape
p = range(n_samples)

random.seed(0)
random.shuffle(p)
X, y = X[p], y[p]
half = int(n_samples / 2)

# Add noisy features
np.random.seed(0)
X = np.c_[X, np.random.randn(n_samples, 200 * n_features)]

# Run classifier
classifier = svm.SVC(kernel='linear', probability=True)
probas_ = classifier.fit(X[:half], y[:half]).predict_proba(X[half:])
y_pred = classifier.predict(X[half:])
y_true = y[half:]


def test_roc_curve():
    """test Receiver operating characteristic (ROC)"""
    fpr, tpr, thresholds = roc_curve(y[half:], probas_[:,1])
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.80, decimal=2)


def test_precision_recall_curve():
    """test Precision-Recall"""
    precision, recall, thresholds = precision_recall_curve(y[half:],
                                                           probas_[:,1])
    precision_recall_auc = auc(precision, recall)
    assert_array_almost_equal(precision_recall_auc, 0.3197, 3)


def test_confusion_matrix():
    """test confusion matrix"""
    cm = confusion_matrix(y_true, y_pred)
    assert_array_equal(cm, [[19, 6],[7, 18]])


def test_losses():
    """test loss functions"""
    assert_equal(zero_one(y_true, y_pred), 13)
    assert_almost_equal(mean_square_error(y_true, y_pred), 12.999, 2)
    assert_almost_equal(mean_square_error(y_true, y_true), 0.00, 2)

    assert_almost_equal(explained_variance(y_true, y_pred), -0.04, 2)
    assert_almost_equal(explained_variance(y_true, y_true), 1.00, 2)


def test_symmetry():
    """test the symmetry of score and loss functions"""
    # symmetric
    assert_equal(zero_one(y_true, y_pred),
                 zero_one(y_pred, y_true))
    assert_almost_equal(mean_square_error(y_true, y_pred),
                        mean_square_error(y_pred, y_true))
    # not symmetric
    assert_(explained_variance(y_true, y_pred) != \
            explained_variance(y_pred, y_true))
    # FIXME: precision and recall aren't symmetric either


def test_precision_recall_multilabel():
    # temporary disabling multilabel support to implement proper multiclass
    # support first
    raise nose.SkipTest("XFailed Test")

    # Y[i, j] = 1 means sample i has label j
    Y_true = np.array([[1, 0, 1, 0],
                       [1, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 1, 1, 1]])

    Y_pred = np.array([[1, 1, 1, 0],
                       [1, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 1, 1]])

    n_pred = 8.0
    n_corr_pred = 6.0
    n_labeled = 7.0
    p = n_corr_pred / n_pred
    r = n_corr_pred / n_labeled
    f1 = 2 * p * r / (p + r)

    assert_equal(p, precision(Y_true, Y_pred))
    assert_equal(r, recall(Y_true, Y_pred))
    assert_equal((p,r), precision_recall(Y_true, Y_pred))
    assert_equal(f1, f1_score(Y_true, Y_pred))


