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
clf2 = Perceptron()

clf1.partial_setup(n_features, n_labels).partial_fit(X, y)
clf2.partial_setup(n_features, n_labels, averaged=True).partial_fit(X, y)

y_pred1 = clf1.predict(X)
y_pred2 = clf2.predict(X)

print "Num. of mislabeled points (perceptron)      : %d" % (y != y_pred1).sum()
print "Num. of mislabeled points (avg. perceptron) : %d" % (y != y_pred2).sum()
