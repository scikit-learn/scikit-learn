"""
=================================================================
Test with permutations the significance of a classification score
=================================================================

This example demonstrates the use of
:func:`~sklearn.model_selection.permutation_test_score` to evaluate the
significance of a cross-valdiated score using permutations.

In order to test if a classification score is significative a technique
in repeating the classification procedure after randomizing, permuting,
the labels. The p-value is then given by the percentage of runs for
which the score obtained is greater than the classification score
obtained in the first place.
"""

# Author:  Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause

print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import permutation_test_score
from sklearn.datasets import load_iris


# %%
# Dataset
# -------
#
# We will use the :ref:`iris_dataset`, which consists of measurements taken
# from 3 types of irises.

iris = load_iris()
X = iris.data
y = iris.target
n_classes = np.unique(y).size

# %%
# We will also generate some random data, uncorrelated with the class labels in
# the iris dataset.

rng = np.random.RandomState(seed=0)
X_rand = rng.normal(size=(len(X), 2200))

# %%
# Permutation test score
# ----------------------
#
# Next, we calculate the
# :func:`~sklearn.model_selection.permutation_test_score` using the original
# iris dataset, which has dependency between features and labels, and
# the randomly generated features and iris labels, which should have
# no dependency between features and labels. The
# :class:`~sklearn.svm.svc` classifier and :ref:`accuracy_score` are used.
#
# :func:`~sklearn.model_selection.permutation_test_score` generates a null
# distribution by calculating the accuracy of the classifier
# on 1000 different permutations of the dataset, where features
# remain the same but labels undergo different permutations. This is the
# distribution for the null hypothesis that there is no dependency between
# the features and labels. An empirical p value is then calculated using
# the null distribution and the score obtained using the original data.

clf = SVC(kernel='linear')
cv = StratifiedKFold(2)

score_iris, perm_scores_iris, pvalue_iris = permutation_test_score(
    clf, X, y, scoring="accuracy", cv=cv, n_permutations=1000, n_jobs=1)

score_rand, perm_scores_rand, pvalue_rand = permutation_test_score(
    clf, X, y, scoring="accuracy", cv=cv, n_permutations=1000, n_jobs=1)

# %%
# Original data
# ^^^^^^^^^^^^^
#
# Below we plot a histogram of the ``permutation_scores`` (the null
# distribution). The red line indicates the score obtained by the classifier
# on the original data. The score is much better than those obtained by
# using permuted data and the p value is thus very low. This indicates that
# there is a low likelihood that this good score would be obtained by chance
# alone. It provides evidence that the iris dataset contains real dependency
# between features and labels and the classifier was able to utilise this
# to obtain good results.

fig, ax = plt.subplots()

ax.hist(perm_scores_iris, bins=20)
ax.axvline(score, ls='--', color='r')
score_label = (f"Score on original\ndata: {score_iris:.2f}\n"
               f"(p value: {pvalue_iris:.3f})")
ax.text(0.7, 260, score_label, fontsize=12)
plt.show()

# %%
# Random data
# ^^^^^^^^^^^
#
# Below we plot the null distribution for the randomized data. The permutation
# scores are similar to those obtained using the original iris dataset
# because the permutation always destroys any feature label dependency present.
# The score obtained on the original randomized data in this case though, is
# very poor. This results in a large p value, confirming that there was no
# feature label dependency in the original data.

fig, ax = plt.subplots()

ax.hist(permutation_scores, bins=20)
ax.set_xlim(0.13)
ax.axvline(score, ls='--', color='r')
score_label = (f"Score on original\ndata: {score:.2f}\n"
               f"(p value: {pvalue:.3f})")
ax.text(0.14, 125, score_label, fontsize=12)
