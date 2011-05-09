"""
============================
Averaged perceptrons
============================

Online learning example using (averaged) perceptrons.

"""

################################################################################
# import some data to play with

from scikits.learn import datasets
digits = datasets.load_digits()

X = digits.data
y = digits.target

################################################################################
# Perceptrons
from scikits.learn.perceptron import Perceptron
import numpy as np

# Although not required for the digits dataset (which is loaded into memory in
# its entirety), this snippet demonstrates the online learning API. We derive
# the number of features and labels/classes from the dataset, but usually you'd
# need to manually specify these when doing online learning.
n_features = X.shape[1]
n_labels   = len(np.unique(y))

clf1 = Perceptron()
clf2 = Perceptron(averaged=True)

clf1.partial_setup(n_features, n_labels).partial_fit(X, y)
clf2.fit(X, y)  # equivalent to above but with batch learning interface

def accuracy(clf):
    y_pred = clf.predict(X)
    return (y == y_pred).sum() / float(len(y))

print('Accuracy of perceptron:          %.2f' % (accuracy(clf1) * 100.))
print('Accuracy of averaged perceptron: %.2f' % (accuracy(clf2) * 100.))
