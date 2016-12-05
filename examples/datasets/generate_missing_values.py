"""
=============================================
Generating NMAR / MCAR missing_values in data
=============================================

This example illustrates how the :class:`sklearn.datasets.ValueDropper` can
be used to generate missing values completely at random or based on the
given drop-probabilities.

The :class`sklearn.datasets.ValueDropper` is a transformer which can be
initialized with a ``missing_proba`` specifying the drop-probabilites
for each class label (and each feature if needed). This facilitates
benchmarking missing-value strategies and evaluating the performance of such
strategies with respect to the type, extent and distribution of missingness in
the data. Importantly, when ``random_state`` is set to an integer, it
provisions preserving the drop-locations as the ``missing_proba`` is upscaled
to study the effect of the increase in missingness. This allows benchmarking
with incremental missing rates without causing variation in the results due to
an inconsistency in the drop-locations between different scales of
``missing_proba``.

NMAR or Not Missing At Random refers to the case when the missingness in the
data is distributed not at random. It is either correlated with the target
value(s) or with the data itself. In some references it is also refered to as
MNAR or Missing Not At Random.

MCAR or Missing Completely At Random refers to the case when the missingness
in the data is completely random and does not correlate with the classification
target value(s) or the data.
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

# For samples from class 1, each feature will be missing 20% of its values
vd = ValueDropper(missing_proba={1: 0.2}, random_state=0)
X_dropped = vd.transform(X, y)

print("\nAfter dropping 20% of values (per feature) in samples of class 1:")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Each feature of samples of class 1 will have a further 20% of its values
# missing. (Old locations will be preserved as random_state is set)
vd = ValueDropper(missing_proba={1: 0.4}, random_state=0)
X_dropped = vd.transform(X, y)

print("\nAfter dropping another 20% of values (per feature) in samples of "
      "class 1:")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop 30% of values in each feature completely at random

vd = ValueDropper(missing_proba=0.3, random_state=0)
X_dropped = vd.transform(X, y)

print("\nAfter dropping 30% of values randomly:")
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop values based on the given drop-probabilities -

# For samples of class 0, drop 10% of values (in each feature)
# For samples of class 2, drop 20% of values in feature 0, 40% in feature 1
#                         and None in feature 2
# Don't drop any values for samples of class 1.
missing_proba = {0: 0.1, 2: [0.2, 0.4, 0]}
vd = ValueDropper(missing_proba=missing_proba, random_state=0)
X_dropped = vd.transform(X, y)

print("\nAfter dropping one set of missing values based on the "
      "missing_proba=%s" % missing_proba)
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop twice as many missing values as in previous step.
missing_proba = {0: 0.2, 2: [0.4, 0.6, 0]}
vd = ValueDropper(missing_proba=missing_proba, random_state=0)
X_dropped = vd.transform(X, y)
print("\nAfter dropping another set of missing values based on the new "
      "missing_proba=%s" % missing_proba)
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")

# Drop more values and also drop 40% of values from samples of class 1
# (in each feature)
missing_proba = {0: 0.3, 1: 0.4, 2: [0.6, 0.8, 0]}
vd = ValueDropper(missing_proba=missing_proba, random_state=0)
X_dropped = vd.transform(X, y)
print("\nAfter dropping another set of missing values based on the new "
      "missing_proba=%s" % missing_proba)
print("y", "X", sep="\t")
print("------------------------")
for i in range(y.shape[0]):
    print(y[i], X_dropped[i], sep="\t")
