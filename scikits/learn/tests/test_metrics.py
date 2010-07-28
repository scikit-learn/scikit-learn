import random
import numpy as np
from scikits.learn import svm, datasets
from scikits.learn.metrics import roc, auc, precision_recall, \
            confusion_matrix, zero_one, explained_variance, \
            mean_square_error

from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_equal, assert_almost_equal

# import some data to play with
iris = datasets.load_iris()
X = iris.data
y = iris.target
X, y = X[y!=2], y[y!=2]
n_samples, n_features = X.shape
p = range(n_samples)
random.seed(0)
random.shuffle(p)
X, y = X[p], y[p]
half = int(n_samples/2)

# Add noisy features
np.random.seed(0)
X = np.c_[X,np.random.randn(n_samples, 200*n_features)]

# Run classifier
classifier = svm.SVC(kernel='linear', probability=True)
probas_ = classifier.fit(X[:half],y[:half]).predict_proba(X[half:])
y_ = classifier.predict(X[half:])

def test_roc():
    """test Receiver operating characteristic (ROC)"""
    fpr, tpr, thresholds = roc(y[half:], probas_[:,1])
    roc_auc = auc(fpr, tpr)
    assert_array_almost_equal(roc_auc, 0.80, decimal=2)

def test_precision_recall():
    """test Precision-Recall"""
    precision, recall, thresholds = precision_recall(y[half:], probas_[:,1])
    precision_recall_auc = auc(precision, recall)
    assert_array_almost_equal(precision_recall_auc, 0.3197, 3)

def test_confusion_matrix():
    """test confusion matrix"""
    cm = confusion_matrix(y[half:], y_)
    assert_array_equal(cm, [[19, 6],[7, 18]])

def test_losses():
    """test loss functions"""
    assert_equal(zero_one(y[half:], y_), 13)
    assert_almost_equal(mean_square_error(y[half:], y_), 12.999, 2)
    assert_almost_equal(explained_variance(y[half:], y_), -0.04, 2)
