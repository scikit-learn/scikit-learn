"""
=================================================================
Test with permutations the significance of a classification score
=================================================================

This example demonstrates the use of
:func:`~sklearn.model_selection.permutation_test_score` to evaluate the
significance of a cross-valdiated score using permutations.
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
# We will also generate some random feature data, uncorrelated with the class
# labels in the iris dataset.

rng = np.random.RandomState(seed=0)
X_rand = rng.normal(size=(len(X), 2200))

# %%
# Permutation test score
# ----------------------
#
# Next, we calculate the
# :func:`~sklearn.model_selection.permutation_test_score` using the original
# iris dataset, which has strong structure, and
# the randomly generated features and iris labels, which should have
# no dependency between features and labels. The
# :class:`~sklearn.svm.svc` classifier and :ref:`accuracy_score` are used.
#
# :func:`~sklearn.model_selection.permutation_test_score` generates a null
# distribution by calculating the accuracy of the classifier
# on 1000 different permutations of the dataset, where features
# remain the same but labels undergo different permutations. This is the
# distribution for the null hypothesis that there is no dependency between
# the features and labels. An empirical p value is then calculated as
# the percentage of permutations for which the score obtained is greater
# that the score obtained using the original data.

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
# between features and labels and the classifier was able to utilize this
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

# %%
# Another possible reason for a obtaining high p value is that the classifier
# was not able to use the structure in the data. In
# this case the p value would only be low for classifiers that are able to
# utilize the dependency present. In our case above, where the data is random,
# all classifiers would have a high p value as there is no structure present
# in the data.
#
# Finally, note that this test has been shown to produce low p values even
# if there is only weak structure in the data[1]_.
#
# .. topic:: References:
#
#   .. [1] Ojala and Garriga. `Permutation Tests for Studying Classifier
#       Performance
#       <http://www.jmlr.org/papers/volume11/ojala10a/ojala10a.pdf>`_. The
#       Journal of Machine Learning Research (2010) vol. 11
