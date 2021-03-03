"""
==========================================
Feature importances with a forest of trees
==========================================

This example shows the use of a forest of trees to evaluate the importance of
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
from sklearn.model_selection import train_test_split

X, y = make_classification(
    n_samples=1000, n_features=10, n_informative=3, n_redundant=0,
    n_repeated=0, n_classes=2, random_state=0, shuffle=False)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, random_state=42)

# %%
# A random forest classifier will be fitted to compute the feature importances.
from sklearn.ensemble import RandomForestClassifier

feature_names = [f'feature {i}' for i in range(X.shape[1])]
forest = RandomForestClassifier(random_state=0)
forest.fit(X_train, y_train)

# %%
# Feature importance based on mean decrease in impurity
# -----------------------------------------------------
# Feature importances are provided by the fitted attribute
# `feature_importances_` and they are computed as the mean and standard
# deviation of accumulation of the impurity decrease within each tree.
#
# .. warning::
#     Impurity-based feature importances can be misleading for high cardinality
#     features (many unique values). See :ref:`permutation_importance` as
#     an alternative below.
import time
import numpy as np

start_time = time.time()
importances = forest.feature_importances_
std = np.std([
    tree.feature_importances_ for tree in forest.estimators_], axis=0)
elapsed_time = time.time() - start_time

print(f"Elapsed time to compute the importances: "
      f"{elapsed_time:.3f} seconds")

# %%
# Let's plot the impurity-based importance.
import pandas as pd
forest_importances = pd.Series(importances, index=feature_names)

fig, ax = plt.subplots()
forest_importances.plot.bar(yerr=std, ax=ax)
ax.set_title("Feature importances using MDI")
ax.set_ylabel("Mean decrease in impurity")
fig.tight_layout()

# %%
# We observe that, as expected, the three first features are found important.
#
# Feature permutation importances on OOB samples
# ----------------------------------------------
# We will an alternative to the impurity-based feature importances based on
# feature permutation using the OOB samples. We fit a new random-forest where
# we explicitely specify to compute the permutation feature importances on OOB.
feature_names = [f'feature {i}' for i in range(X.shape[1])]
forest = RandomForestClassifier(feature_importances="permutation_oob",
                                random_state=0)
start_time = time.time()
forest.fit(X_train, y_train)
elapsed_time = time.time() - start_time

print(f"Elapsed time to compute the importances: "
      f"{elapsed_time:.3f} seconds")

forest_importances = pd.Series(forest.feature_importances_, index=feature_names)

# %%
# The permutation importances is more computationally costly. Indeed, it
# requires to fit the tree and to make additional processing: each tree will
# be evaluated on its OOB sample as well as an OOB sample where features will
# be permuted. This step is costly and explains the time fitting difference
# between of the two forests.
#
# We now plot the feature importance ranking.
fig, ax = plt.subplots()
forest_importances.plot.bar(yerr=forest.importances_.importances_std, ax=ax)
ax.set_title("Feature importances using permutation on full model")
ax.set_ylabel("Mean accuracy decrease")
fig.tight_layout()
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
