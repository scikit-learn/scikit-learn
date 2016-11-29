"""
=============================================
Generating NMAR / MCAR missing_values in data
=============================================

This example illustrates how the :class:`sklearn.datasets.ValueDropper` can
be used to generate missing values completely at random or conforming to the
given distribution.

The :class`sklearn.datasets.ValueDropper` is a transformer which can be
initialized with a ``missing_proba`` specifying the drop probabilites
for each label (and each feature if needed). It provisions preserving the
missing values of lower scaled ``missing_proba`` in a higher scaled
``missing_proba``. This facilitates benchmarking missing-value
strategies and evaluating the performance of such strategies with
respect to the type and extent of missingness in data.

It allows benchmarking with incremental missing rates (fraction of missing
values to total number of values) without introducing a mismatch in the
missing positions for previous lower rates of missing values.

NMAR or Not Missing At Random refers to the case when the missingness in the
data is distributed not at random. It is either correlated with the target
value(s) or with the data itself.

MCAR or Missing Completely At Random refers to the case when the missingness
in the data is completely random and does not correlate with the classification
target value(s) or the data.

In some references NMAR is sometimes referred to as MNAR (Missing Not At
Random).
"""
# Author: Raghav RV <rvraghav93@gmail.com>
#
# License: BSD 3 clause

from __future__ import print_function

import numpy as np
from sklearn.datasets import ValueDropper

print(__doc__)


X = np.random.RandomState(0).random_sample((20, 3))
y = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2])

# 10% of values in samples of class 1 will have missing values across all
# features
vd = ValueDropper(missing_proba={1: 0.2}, random_state=0)
X_dropped = vd.transform(X, y)

print("\nAfter dropping 10% of values in samples of class 1 across all"
      " features")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop another 10% of values in samples of class 1 across all features
vd = ValueDropper(missing_proba={1: 0.4}, random_state=0)
X_dropped = vd.transform(X, y)

print("\nAfter dropping another 10% of values in samples of class 1 across all"
      " features")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop 30% of values in all features completely at random

vd = ValueDropper(missing_proba=0.3, random_state=0)
X_dropped = vd.transform(X, y)

print("\nAfter dropping 30% of values randomly")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop values based on the given probability distribution

# For samples of class 0, drop 10% of values (evenly across all features)
# For samples of class 1, drop 20% of values in feature 0, 40% in feature 1
#                         and None in feature 2
# Don't drop any values for samples of class 2.
missing_proba = {0: 0.1, 1: [0.2, 0.4, 0]}
vd = ValueDropper(missing_proba=missing_proba, random_state=0)
X_dropped = vd.transform(X, y)

print("The given class wise missing_proba dict is %s " % missing_proba)
print("\nAfter dropping one set of missing values based on the "
      " missing_proba=%s" % missing_proba)
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop twice as many missing values as in previous step.
missing_proba = {0: 0.2, 1: [0.4, 0.6, 0]}
vd = ValueDropper(missing_proba=missing_proba, random_state=0)
X_dropped = vd.transform(X, y)
print("\nAfter dropping another set of missing values based on the new"
      " missing_proba=%s" % missing_proba)
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

missing_proba = {0: 0.3, 1: [0.6, 0.8, 0]}
vd = ValueDropper(missing_proba=missing_proba, random_state=0)
X_dropped = vd.transform(X, y)
print("\nAfter dropping another set of missing values based on the new"
      " missing_proba=%s" % missing_proba)
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")
