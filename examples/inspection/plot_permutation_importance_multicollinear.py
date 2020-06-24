"""
=================================================================
Permutation Importance with Multicollinear or Correlated Features
=================================================================

In this example, we compute the permutation importance on the Wisconsin
breast cancer dataset using :func:`~sklearn.inspection.permutation_importance`.
The :class:`~sklearn.ensemble.RandomForestClassifier` can easily get about 97%
accuracy on a test dataset. Because this dataset contains multicollinear
features, the permutation importance will show that none of the features are
important. One approach to handling multicollinearity is by performing
hierarchical clustering on the features' Spearman rank-order correlations,
picking a threshold, and keeping a single feature from each cluster.

.. note::
    See also
    :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`
"""
print(__doc__)
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from scipy.cluster import hierarchy

from sklearn.datasets import load_breast_cancer
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split

##############################################################################
# Random Forest Feature Importance on Breast Cancer Data
# ------------------------------------------------------
# First, we train a random forest on the breast cancer dataset and evaluate
# its accuracy on a test set:
data = load_breast_cancer()
X, y = data.data, data.target
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)
print("Accuracy on test data: {:.2f}".format(clf.score(X_test, y_test)))

##############################################################################
# Next, we plot the tree based feature importance and the permutation
# importance. The permutation importance plot shows that permuting a feature
# drops the accuracy by at most `0.012`, which would suggest that none of the
# features are important. This is in contradiction with the high test accuracy
# computed above: some feature must be important. The permutation importance
# is calculated on the training set to show how much the model relies on each
# feature during training.
result = permutation_importance(clf, X_train, y_train, n_repeats=10,
                                random_state=42)
perm_sorted_idx = result.importances_mean.argsort()

tree_importance_sorted_idx = np.argsort(clf.feature_importances_)
tree_indices = np.arange(0, len(clf.feature_importances_)) + 0.5

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
ax1.barh(tree_indices,
         clf.feature_importances_[tree_importance_sorted_idx], height=0.7)
ax1.set_yticklabels(data.feature_names[tree_importance_sorted_idx])
ax1.set_yticks(tree_indices)
ax1.set_ylim((0, len(clf.feature_importances_)))
ax2.boxplot(result.importances[perm_sorted_idx].T, vert=False,
            labels=data.feature_names[perm_sorted_idx])
fig.tight_layout()
plt.show()

##############################################################################
# Handling Multicollinear Features
# --------------------------------
# When features are collinear, permutating one feature will have little
# effect on the models performance because it can get the same information
# from a correlated feature. One way to handle multicollinear features is by
# performing hierarchical clustering on the Spearman rank-order correlations,
# picking a threshold, and keeping a single feature from each cluster. First,
# we plot a heatmap of the correlated features:
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
corr = spearmanr(X).correlation
corr_linkage = hierarchy.ward(corr)
dendro = hierarchy.dendrogram(
    corr_linkage, labels=data.feature_names.tolist(), ax=ax1, leaf_rotation=90
)
dendro_idx = np.arange(0, len(dendro['ivl']))

ax2.imshow(corr[dendro['leaves'], :][:, dendro['leaves']])
ax2.set_xticks(dendro_idx)
ax2.set_yticks(dendro_idx)
ax2.set_xticklabels(dendro['ivl'], rotation='vertical')
ax2.set_yticklabels(dendro['ivl'])
fig.tight_layout()
plt.show()

##############################################################################
# Next, we manually pick a threshold by visual inspection of the dendrogram
# to group our features into clusters and choose a feature from each cluster to
# keep, select those features from our dataset, and train a new random forest.
# The test accuracy of the new random forest did not change much compared to
# the random forest trained on the complete dataset.
cluster_ids = hierarchy.fcluster(corr_linkage, 1, criterion='distance')
cluster_id_to_feature_ids = defaultdict(list)
for idx, cluster_id in enumerate(cluster_ids):
    cluster_id_to_feature_ids[cluster_id].append(idx)
selected_features = [v[0] for v in cluster_id_to_feature_ids.values()]

X_train_sel = X_train[:, selected_features]
X_test_sel = X_test[:, selected_features]

clf_sel = RandomForestClassifier(n_estimators=100, random_state=42)
clf_sel.fit(X_train_sel, y_train)
print("Accuracy on test data with features removed: {:.2f}".format(
      clf_sel.score(X_test_sel, y_test)))
