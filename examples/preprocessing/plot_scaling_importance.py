#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
Importance of Feature Scaling
=========================================================

Features scaling though standardization (or Z-score normalization)
can be an importance preprocessing step for many machine learning
algorithms. Standardization involves rescaling the features such
that they’ll have the properties of a standard normal distribution
with a mean of zero and a standard deviation of one.

While many algorithms (such as SVM, K-nearest neighbors and logistic
regression) require features to be normalized, intuitively we can
think of Principle Component Analysis (PCA) as being a prime example
of when normalization is important. In PCA we are interested in the
components that maximize the variance. If there exists components
(e.g human height) that vary less then other components (e.g human
weight) because of their respective scales (meters vs. kilos) it can
be seen how not scaling the features would cause PCA to determine that
the direction of maximal variance more closely corresponds with the
‘weight’ axis. As a change in height of one meter can be considered much
more important than the change in weight of one kilogram, it is easily
seen that this determination is incorrect. In the case of PCA, scaling
features using normalization is preferred over using min-max scaling as
the primary components are computed using the correlation matrix as opposed
to the covariance matrix.

In order to illustrate this in an example, PCA will be performed on a dataset
which has been standardized using :class:`StandardScaler <sklearn.preprocessing.StandardScaler>`,
and a copy which has remained untouched. The results with be visualized and
a clear difference noted.

The results will then be used to train a naive Bayes classifier, and a clear
difference the prediction accuracies will be observed.

"""
from __future__ import print_function
print(__doc__)


# Code source: Tyler Lanigan <tylerlanigan@gmail.com>
# 			   Sebastian Raschka <mail@sebastianraschka.com>

# License: BSD 3 clause

from sklearn.cross_validation import train_test_split
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.datasets import load_wine

# Contants
RAN_STATE = 42
FIG_SIZE = (10, 7)


features, target = load_wine(return_X_y=True)

# Make a train/test split using 30% test size
X_train, X_test, y_train, y_test = train_test_split(features, target,
                                                    test_size=0.30, random_state=RAN_STATE)

# Apply Scaling to X_train and X_test
std_scale = preprocessing.StandardScaler().fit(X_train)
X_train_std = std_scale.transform(X_train)
X_test_std = std_scale.transform(X_test)

# Perform PCA on non-standardized data
pca = PCA(n_components=2).fit(X_train)
X_train = pca.transform(X_train)
X_test = pca.transform(X_test)

# Perform PCA on standardized data
pca_std = PCA(n_components=2).fit(X_train_std)
X_train_std = pca_std.transform(X_train_std)
X_test_std = pca_std.transform(X_test_std)

# Fit GaussianNB on standard and non-standardized data
clf = GaussianNB()
fit = clf.fit(X_train, y_train)
clf_std = GaussianNB()
fit_std = clf_std.fit(X_train_std, y_train)

# Make predictions for standard and non standardized data.
pred_train = clf.predict(X_train)
pred_test = clf.predict(X_test)
pred_train_std = clf_std.predict(X_train_std)
pred_test_std = clf_std.predict(X_test_std)

print('\nPrediction accuracy for the normal test dataset with PCA')
print('{:.2%}\n'.format(metrics.accuracy_score(y_test, pred_test)))

print('\nPrediction accuracy for the standardized test dataset with PCA')
print('{:.2%}\n'.format(metrics.accuracy_score(y_test, pred_test_std)))


# visualize standardized vs. untouched dataset with PCA performed

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=FIG_SIZE)


for l, c, m in zip(range(0, 3), ('blue', 'red', 'green'), ('^', 's', 'o')):
    ax1.scatter(X_train[y_train == l, 0], X_train[y_train == l, 1],
                color=c,
                label='class %s' % l,
                alpha=0.5,
                marker=m
                )

for l, c, m in zip(range(0, 3), ('blue', 'red', 'green'), ('^', 's', 'o')):
    ax2.scatter(X_train_std[y_train == l, 0], X_train_std[y_train == l, 1],
                color=c,
                label='class %s' % l,
                alpha=0.5,
                marker=m
                )

ax1.set_title('Training dataset after PCA')
ax2.set_title('Standardized training dataset after PCA')

for ax in (ax1, ax2):

    ax.set_xlabel('1st principal component')
    ax.set_ylabel('2nd principal component')
    ax.legend(loc='upper right')
    ax.grid()
plt.tight_layout()

plt.show()
