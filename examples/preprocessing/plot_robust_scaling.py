#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
Robust Scaling on Toy Data
=========================================================

Making sure that each Feature has approximately the same scale can be a
crucial preprocessing step. However, when data contains outliers,
:class:`StandardScaler <sklearn.preprocessing.StandardScaler>` can often
be mislead. In such cases, it is better to use a scaler that is robust
against outliers.

Here, we demonstrate this on a toy dataset, where one single datapoint
is a large outlier.
"""
from __future__ import print_function
print(__doc__)


# Code source: Thomas Unterthiner
# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler

# Create training and test data
np.random.seed(42)
n_datapoints = 100
C=[[0.9, 0.0],[0.0, 20.0]]
mu1 = [100.0, -3.0]
mu2 = [101.0, -3.0]
X1 = np.random.multivariate_normal(mean=mu1, cov=C, size=n_datapoints)
X2 = np.random.multivariate_normal(mean=mu2, cov=C, size=n_datapoints)
Y_train = np.hstack([[-1]*n_datapoints, [1]*n_datapoints])
X_train = np.vstack([X1, X2])

X1 = np.random.multivariate_normal(mean=mu1, cov=C, size=n_datapoints)
X2 = np.random.multivariate_normal(mean=mu2, cov=C, size=n_datapoints)
Y_test = np.hstack([[-1]*n_datapoints, [1]*n_datapoints])
X_test = np.vstack([X1, X2])

X_train[0, 0] = -1000  # a fairly large outlier


# Scale data
s = StandardScaler()
Xtr_s = s.fit_transform(X_train)
Xte_s = s.transform(X_test)

r = RobustScaler()
Xtr_r = r.fit_transform(X_train)
Xte_r = r.fit_transform(X_test)


# Plot data (note we zoom in to the data center, so the outlier
# data point can't be seen in this plot).
fig, ax = plt.subplots(1, 2, figsize=(8, 4), sharex=True, sharey=True)
ax[0].scatter(Xtr_s[:, 0], Xtr_s[:, 1], color=np.where(Y_train > 0, 'r', 'b'))
ax[1].scatter(Xtr_r[:, 0], Xtr_r[:, 1], color=np.where(Y_train > 0, 'r', 'b'))
ax[0].set_title("Training set after standard scaling")
ax[1].set_title("Training set after robust scaling")
ax[0].set_xlim(-3, 3)
ax[0].set_ylim(-3, 3)
plt.tight_layout()
plt.show()


# Classify using k-NN
from sklearn.neighbors import KNeighborsClassifier

knn = KNeighborsClassifier()
knn.fit(Xtr_s, Y_train)
acc_s = knn.score(Xte_s, Y_test)
print("Testset accuracy using standard scaler: %.3f" % acc_s)
knn.fit(Xtr_r, Y_train)
acc_r = knn.score(Xte_r, Y_test)
print("Testset accuracy using robust scaler:   %.3f" % acc_r)
