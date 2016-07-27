"""
=============================================
Generating NMAR / MCAR missing_values in data
=============================================

This example illustrates how the :class:`sklearn.datasets.ValueDropper` can
be used to generate missing values completely at random or conforming to the
given distribution.

The :class`sklearn.datasets.ValueDropper` is a transformer which can be
initialized with a ``missing_distribution`` specifying the drop probabilites
for each label (and each feature if needed). It provisions preserving the
missing values of lower scaled ``missing_distribution`` in a higher scaled
``missing_distribution``. This facilitates benchmarking missing-value
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


X = np.array([[0, 1, 2],
              [3, 4, 5],
              [6, 7, 8],
              [9, 0, 1],
              [2, 3, 4],
              [8, 9, 8],
              [8, 9, 8],
              [1, 0, 5],
              [5, 4, 3],
              [2, 1, 1],
              [3, 4, 5],
              [2, 3, 4],
              [8, 9, 8],
              [7, 8, 9]], dtype=float)
y = np.array([1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2])

# Drop 10% of values across all features, where all missing values
# come from samples of class 1

vd = ValueDropper(missing_distribution={1: 0.1}, random_state=42)
X_dropped = vd.transform(X, y)

print("\nAfter dropping 10% of values when class label(s) are 1")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop another 10% of values across all features, where all missing values
# come from samples of class 1

vd = ValueDropper(missing_distribution={1: 0.2}, random_state=42)
X_dropped = vd.transform(X, y)

print("\nAfter dropping another 10% of values when class label(s) are 1")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop 30% of values completely at random

vd = ValueDropper(missing_distribution=0.3, random_state=42)
X_dropped = vd.transform(X, y)

print("\nAfter dropping 30% of values randomly")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop 30% of values but according to the given probability distribution
# Incrementally adding 10% each time

# 40% of the dropped values must be from class label 0
# (evenly across all features)
# The rest 60% of the dropped values are from class label 1, distributed in the
# 1:2:0 ratio amongst the features.
# Don't drop any values from samples of class 2
abs_missing_rate = 0.1
missing_distribution = {0: 0.4 * abs_missing_rate,
                        1: np.array([0.2, 0.4, 0]) * abs_missing_rate}

# Also let's use -1 to denote missing values, this time
vd = ValueDropper(missing_values=-1, missing_distribution=missing_distribution,
                  random_state=42)
X_dropped = vd.transform(X, y)

print("The given class wise distribution is %s " % missing_distribution)
print("\nAfter dropping 10% of values according to the distribution")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# NOTE that the relative values of the distribution must not be changed.
abs_missing_rate = 0.2
missing_distribution = {0: 0.4 * abs_missing_rate,
                        1: np.array([0.2, 0.4, 0]) * abs_missing_rate}
vd = ValueDropper(missing_values=-1, missing_distribution=missing_distribution,
                  random_state=42)
X_dropped = vd.transform(X, y)
print("\nAfter dropping another 10% of values according to the distribution")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

abs_missing_rate = 0.3
missing_distribution = {0: 0.3 * abs_missing_rate,
                        1: np.array([0.2, 0.4, 0]) * abs_missing_rate}
vd = ValueDropper(missing_values=-1, missing_distribution=missing_distribution,
                  random_state=42)
X_dropped = vd.transform(X, y)
print("\nAfter dropping another 10% of values according to the distribution")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")
