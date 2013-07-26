"""
=======================================
Receiver operating characteristic (ROC)
=======================================

Example of Receiver operating characteristic (ROC) metric to
evaluate the quality of the output of a classifier.

ROC curves typically feature true positive rate on the Y axis,
and false positive rate on the X axis. This means that the top left corner
of the plot is considered the "ideal" point - a false positive rate of
zero, and a true positive rate of one.
This is not very realistic, but it does mean that more area under the curve
is usually better.

The "steepness" of ROC curves is also important, since it is ideal to
maximize the true positive rate while minimizing the false positive rate.

.. note::

    See also :ref:`example_plot_roc_crossval.py`

"""
print(__doc__)

import numpy as np
import pylab as pl
from sklearn import svm, datasets
from sklearn.utils import shuffle
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import train_test_split


# Import some data to play with
iris = datasets.load_iris()
X = iris.data
y = iris.target

# Make it a binary classification problem by removing the third class
X, y = X[y != 2], y[y != 2]

# Add noisy features to make the problem harder
np.random.seed(0)
n_samples, n_features = X.shape
X = np.c_[X, np.random.randn(n_samples, 200 * n_features)]

# shuffle and split training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5, random_state=0)

# Run classifier
classifier = svm.SVC(kernel='linear', probability=True, random_state=0)
probas_ = classifier.fit(X_train, y_train).predict_proba(X_test)

# Compute ROC curve and area the curve
fpr, tpr, thresholds = roc_curve(y_test, probas_[:, 1])
roc_auc = auc(fpr, tpr)
print("Area under the ROC curve : %f" % roc_auc)

# Plot ROC curve
pl.clf()
pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
pl.plot([0, 1], [0, 1], 'k--')
pl.xlim([0.0, 1.0])
pl.ylim([0.0, 1.0])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('Receiver operating characteristic example')
pl.legend(loc="lower right")
pl.show()
