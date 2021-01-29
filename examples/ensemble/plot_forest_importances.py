"""
=========================================
Feature importances with forest of trees
=========================================

This example shows the use of forest of trees to evaluate the importance of
features on an artificial classification task.

We show two strategies to estimate the feature importances: (i) the
impurity-based feature importances and (ii) the permutation feature
importances on out-of-bag (OOB) samples.

.. warning::
    Impurity-based feature importances can be misleading for high cardinality
    features (many unique values). Check the documentation of the
    `feature_importances` parameter to have more details regarding the
    alternative as the permutation feature importances.
"""
print(__doc__)
import matplotlib.pyplot as plt

# %%
# We generate a synthetic dataset with only 3 informative features. We will
# explicitely not shuffle the dataset to ensure that the informative features
# correspond to the three first columns of `X`.
from sklearn.datasets import make_classification

X, y = make_classification(
    n_samples=1000, n_features=10, n_informative=3,
    n_redundant=0, n_repeated=0, n_classes=2,
    random_state=0, shuffle=False)
feature_names = [f"#{i + 1}" for i in range(X.shape[1])]

# %%
# Impurity-based feature importances
# ----------------------------------
# We start by fitting a random-forest. We explicitely request to compute the
# impurity-based feature importance. Note that this is the default value.
import pandas as pd
import time
from sklearn.ensemble import RandomForestClassifier

forest = RandomForestClassifier(
    n_estimators=250, feature_importances="impurity", random_state=0)

start_time = time.time()
forest.fit(X, y)
forest_importances = pd.Series(
    forest.feature_importances_, index=feature_names)
elapsed_time = time.time() - start_time
print(f"Elapsed time to fit and compute the importances: "
      f"{elapsed_time:.3f} seconds")

# %%
# Impurity-based feature importances is relatively fast to compute. It only
# requires to fit the forest and store information regarding the different
# splits of trees. When the feature importances is requested, the splits of
# all trees are introspected to compute the mean decrease in impurity (MDI).
#
# Let's plot the feature importances ranking.
ax = forest_importances.plot.bar(yerr=forest.importances_.importances_std)
ax.set_title("Feature importances using MDI")
_ = ax.set_ylabel("Mean impurity decrease")

# %%
# We observe that the three important features are reported correctly. We also
# observe that non-informative features do not have a null importance. Indeed,
# theses features have been used by some of the trees that tend to overfit on
# some noisy samples.
#
# Feature permutation importances on OOB samples
# ----------------------------------------------
# We will an alternative to the impurity-based feature importances based on
# feature permutation using the OOB samples. We fit a new random-forest where
# we explicitely specify to compute the permutation feature importances on OOB.
forest = RandomForestClassifier(
    n_estimators=250, feature_importances="permutation_oob", random_state=0)

start_time = time.time()
forest.fit(X, y)
forest_importances = pd.Series(
    forest.feature_importances_, index=feature_names)
elapsed_time = time.time() - start_time
print(f"Elapsed time to fit and compute the importances: "
      f"{elapsed_time:.3f} seconds")

# %%
# The permutation importances is more computationally costly. Indeed, it
# requires to fit the tree and to make additional processing: each tree will
# be evaluated on its OOB sample as well as an OOB sample where features will
# be permuted. This step is costly and explains the time fitting difference
# between of the two forests.
#
# We now plot the feature importance ranking.
ax = forest_importances.plot.bar(yerr=forest.importances_.importances_std)
ax.set_title("Feature importances using permutation on OOB")
ax.set_ylabel("Mean accuracy decrease")
plt.show()

# %%
# As in the impurity-based, the three most important features are detected.
# We see that non-important features have a mean decrease accuracy of zeros.
# Indeed, permuted these features did not have an impact on the score.
# Another difference between the two feature importances is the scale of the
# reported values. The permutation importances corresponds to a difference of
# scores and it is not further normalized. With impurity-based feature
# importances reported are normalized: they sum of the importances across
# features will sum to 1.
