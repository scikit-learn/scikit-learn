"""
=======================================================
Comparison of feature importances with forests of trees
=======================================================

This examples compares permutation feature importances with the default for
random forest. Feature importances are evaluated on an artificial
classification task. The error bars show the inter-trees variability.

Both methods correctly identify the same 3 informative features. The permuation
method gives a better indication that the remaining features are not
informative rather than suggesting weak importance.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt

from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier

# Build a classification task using 3 informative features
X, y = make_classification(n_samples=1000,
                           n_features=10,
                           n_informative=3,
                           n_redundant=0,
                           n_repeated=0,
                           n_classes=2,
                           random_state=0,
                           shuffle=False)

importances = []
stds = []

# Build a forest and compute the feature importances using both methods
for use_permutation in (False, True):
    forest = RandomForestClassifier(
        n_estimators=250,
        random_state=0,
        permutation_feature_importances=use_permutation
    )

    forest.fit(X, y)
    importances.append(forest.feature_importances_)
    stds.append(
        np.std([tree.feature_importances_ for tree in forest.estimators_],
               axis=0))

# Get the indicies for the importance order in the default method
indices = np.argsort(importances[0])[::-1]

# Plot the feature importances of both methods
plt.figure()
plt.title("Comparison of feature importances")

COLORS = ['r', 'b']
labels = ['default', 'permutation']

bar_offsets = np.arange(len(importances[0])) * 2

for i, (label, importances, std) in enumerate(zip(labels, importances, stds)):
    plt.bar(bar_offsets + i, importances[indices], label=label,
            color=COLORS[i], yerr=std[indices], )

plt.xticks(bar_offsets + 1, indices)
plt.xlim([-0.5, bar_offsets[-1] + 2.5])
plt.xlabel("Feature number")
plt.ylabel("Importance")
plt.legend()
plt.show()
