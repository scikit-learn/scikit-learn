"""
=================================================
Outlier detection with Local Outlier Factor (LOF)
=================================================

The Local Outlier Factor (LOF) algorithm is an unsupervised anomaly detection
method which computes the local density deviation of a given data point with
respect to its neighbors. It considers as outliers the samples that have a
substantially lower density than their neighbors. This example shows how to
use LOF for outlier detection which is the default use case of this estimator
in scikit-learn. Note that when LOF is used for outlier detection it has no
predict, decision_function and score_samples methods. See
:ref:`User Guide <outlier_detection>`: for details on the difference between
outlier detection and novelty detection and how to use LOF for novelty
detection.

The number of neighbors considered, (parameter n_neighbors) is typically
set 1) greater than the minimum number of samples a cluster has to contain,
so that other samples can be local outliers relative to this cluster, and 2)
smaller than the maximum number of close by samples that can potentially be
local outliers.
In practice, such informations are generally not available, and taking
n_neighbors=20 appears to work well in general.
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import LocalOutlierFactor

print(__doc__)

np.random.seed(42)

# Generate train data
X_inliers = 0.3 * np.random.randn(100, 2)
X_inliers = np.r_[X_inliers + 2, X_inliers - 2]

# Generate some outliers
X_outliers = np.random.uniform(low=-4, high=4, size=(20, 2))
X = np.r_[X_inliers, X_outliers]

n_outliers = len(X_outliers)
ground_truth = np.ones(len(X), dtype=int)
ground_truth[-n_outliers:] = -1

# fit the model for outlier detection (default)
clf = LocalOutlierFactor(n_neighbors=20)
# use fit_predict to compute the predicted labels of the training samples
# (when LOF is used for outlier detection, the estimator has no predict,
# decision_function and score_samples methods).
y_pred = clf.fit_predict(X)
n_errors = (y_pred != ground_truth).sum()

# IMPORTANT: for illustration purpose we would like to plot the level sets of
# the decision function however decision_function is not available when
# LOF is used for outlier detection. We thus refit the model in a novelty
# detection mode (with novelty=True) so that we can compute the decision
# function on a grid. Note that when using novelty=True, you MUST not use
# predict, decision_function and score_samples on the training set X as this
# would lead to wrong results. You must only use these methods on new unseen
# data (which are not in the training set)

# refit the model with novelty=True
clf = LocalOutlierFactor(n_neighbors=20, novelty=True)
clf.fit(X)
xx, yy = np.meshgrid(np.linspace(-5, 5, 50), np.linspace(-5, 5, 50))
Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
Z = Z.reshape(xx.shape)

plt.title("Local Outlier Factor (LOF)")
plt.contourf(xx, yy, Z, cmap=plt.cm.Blues_r)

a = plt.scatter(X_inliers[:, 0], X_inliers[:, 1], c='white',
                edgecolor='k', s=20)
b = plt.scatter(X_outliers[:, 0], X_outliers[:, 1], c='red',
                edgecolor='k', s=20)
plt.axis('tight')
plt.xlim((-5, 5))
plt.ylim((-5, 5))
plt.legend([a, b],
           ["true inliers",
            "true outliers"],
           loc="upper left")
plt.xlabel("prediction errors: %d" % (n_errors))
plt.show()
