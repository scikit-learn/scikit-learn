"""
=================================================================
Permutation Importance with Multicollinear or Correlated Features
=================================================================

In this example, we compute the permutation importance on the Wisconsin
breast cancer dataset using :func:`~sklearn.inspection.permutation_importance`.
The :class:`~sklearn.ensemble.RandomForestClassifier` can easily get about 97%
accuracy on a test dataset. Because this dataset contains multicollinear
features, the permutation importance will show that none of the features are
important. Hierarchical clustering is used to visualize the correlations
between features and :class:`~sklearn.feature_selection.CorrelationThreshold`
to filter the features such that all pairwise correlations are below a certain
threshold.

.. note::
    See also
    :ref:`sphx_glr_auto_examples_inspection_plot_permutation_importance.py`
"""
print(__doc__)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from scipy.cluster import hierarchy

from sklearn.datasets import load_breast_cancer
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import CorrelationThreshold
from sklearn.pipeline import Pipeline
from sklearn.model_selection import RandomizedSearchCV


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
print("Accuracy on test data with 100% of the features: "
      "{:.2f}".format(clf.score(X_test, y_test)))

##############################################################################
# Next, we plot the tree based feature importance and the permutation
# importance. The permutation importance plot shows that permuting a feature
# drops the accuracy by at most `0.012`, which would suggest that none of the
# features are important. This is in contradiction with the high test accuracy
# computed above: some feature must be important. The permutation importance
# is calculated on the training set to show how much the model relies on each
# feature during training.
result = permutation_importance(clf, X_train, y_train, n_repeats=20,
                                random_state=42, n_jobs=2)
perm_sorted_idx = result.importances_mean.argsort()

tree_importance_sorted_idx = np.argsort(clf.feature_importances_)
tree_indices = np.arange(0, len(clf.feature_importances_)) + 0.5

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
ax1.barh(tree_indices,
         clf.feature_importances_[tree_importance_sorted_idx], height=0.7)
ax1.set_yticklabels(data.feature_names)
ax1.set_yticks(tree_indices)
ax1.set_ylim((0, len(clf.feature_importances_)))
ax2.boxplot(result.importances[perm_sorted_idx].T, vert=False,
            labels=data.feature_names)
fig.tight_layout()

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
corr = np.corrcoef(X, rowvar=False)
corr_linkage = hierarchy.ward(corr)
dendro = hierarchy.dendrogram(corr_linkage, labels=data.feature_names, ax=ax1,
                              leaf_rotation=90)
dendro_idx = np.arange(0, len(dendro['ivl']))

ax2.imshow(corr[dendro['leaves'], :][:, dendro['leaves']])
ax2.set_xticks(dendro_idx)
ax2.set_yticks(dendro_idx)
ax2.set_xticklabels(dendro['ivl'], rotation='vertical')
ax2.set_yticklabels(dendro['ivl'])
fig.tight_layout()

##############################################################################
# Next, we use :class:`~sklearn.feature_selection.CorrelationThreshold` to
# select features such that all pairwise correlations are below a certain
# threshold. :class:`~sklearn.model_selection.GridSearchCV` is used to explore
# the different thresholds.
pipe = Pipeline([
    ('correlation', CorrelationThreshold(threshold=0.9)),
    ('forest', RandomForestClassifier(n_estimators=100, random_state=42))])

params = {'correlation__threshold': sp.stats.uniform(0.5, 0.5)}
search = RandomizedSearchCV(pipe, params, refit=False, n_jobs=2, n_iter=30,
                            random_state=42)
search.fit(X_train, y_train)
df = pd.DataFrame(search.cv_results_)
df = df.sort_values(by='rank_test_score')

fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar('rank_test_score', 'mean_test_score', yerr='std_test_score',
            data=df, marker='o')
ax.set_xlabel('rank')
ax.set_ylabel('accuracy')
ax.set_xticks(df['rank_test_score'])

##############################################################################
# When taking into account standard deviation, all the cross validated scores
# overlap. Here we take the threshold with the lowest rank and train a model
# with most of the correlated features removed.
threshold = df.iloc[-1]['param_correlation__threshold']
print("All pairwise correlations above {:4f} will be "
      "removed".format(threshold))

pipe = pipe.set_params(correlation__threshold=threshold)
pipe.fit(X_train, y_train)
support_mask = pipe['correlation'].support_mask_
pct_features = 100 * np.sum(support_mask)/X_train.shape[1]
print("Accuracy on test data with {:.2f}% of the features: {:.2f}".format(
      pct_features, pipe.score(X_test, y_test)))

##############################################################################
# Lastly, we use select features from the test set to ensure the pairwise
# correlations is below the choosen threshold and plot the permutation
# importance with some correlation features removed. This shows that
# "worst symmetry" is important to the model.
X_test_sel = pipe['correlation'].transform(X_test)
feature_names = data.feature_names[support_mask]
result = permutation_importance(pipe['forest'], X_test_sel, y_test,
                                n_repeats=20, random_state=42, n_jobs=2)
perm_sorted_idx = result.importances_mean.argsort()
fig, ax = plt.subplots(figsize=(12, 8))
ax.boxplot(result.importances[perm_sorted_idx].T, vert=False,
           labels=feature_names)
fig.tight_layout()
plt.show()
