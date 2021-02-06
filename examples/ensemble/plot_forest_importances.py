"""
=========================================
Feature importances with forests of trees
=========================================

This examples shows the use of forests of trees to evaluate the importance of
features on an artificial classification task. The red bars are
the impurity-based feature importances of the forest,
along with their inter-trees variability.

As expected, the plot suggests that 3 features are informative, while the
remaining are not.

.. warning::
    Impurity-based feature importances can be misleading for high cardinality
    features (many unique values). See
    :func:`sklearn.inspection.permutation_importance` as an alternative.
"""
print(__doc__)

# %%
import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.inspection import permutation_importance

# %%
# Build a classification task using 3 informative features
X, y = make_classification(n_samples=1000,
                           n_features=10,
                           n_informative=3,
                           n_redundant=0,
                           n_repeated=0,
                           n_classes=2,
                           random_state=0,
                           shuffle=False)

# Build a forest and compute the impurity-based feature importances
forest = ExtraTreesClassifier(n_estimators=250,
                              random_state=0)

forest.fit(X, y)
importances = forest.feature_importances_
std = np.std([tree.feature_importances_ for tree in forest.estimators_],
             axis=0)
indices = np.argsort(importances)

# Print the feature ranking
print("Feature ranking:")
for f in range(X.shape[1]):
    feat_index = indices[::-1][f]
    print(f"{f + 1}. feature {feat_index} ({importances[feat_index]:.4f})")

feature_names = [f'feature {i}' for i in range(X.shape[1])]

# %%
# Plot the impurity-based feature importances of the forest
fig = plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.barh(range(X.shape[1]), importances[indices],
         color="r", xerr=std[indices], align='center')
plt.yticks(range(X.shape[1]), np.array(feature_names)[indices])
plt.ylim([-1, X.shape[1]])
plt.title("Feature importances")

# Plot the feature importances based on permutation importance
result = permutation_importance(forest, X, y, n_repeats=10,
                                random_state=42, n_jobs=2)
sorted_idx = result.importances_mean.argsort()
plt.subplot(1, 2, 2)
plt.boxplot(result.importances[sorted_idx].T,
            vert=False,
            labels=np.array(feature_names)[sorted_idx])
plt.title("Permutation Importances")
fig.tight_layout()
plt.show()
