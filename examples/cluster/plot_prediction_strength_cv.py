"""
=====================================================================
Selecting the number of clusters with prediction strength grid search
=====================================================================

Prediction strength is a metric that measures the stability of a clustering
algorithm and can be used to determine an optimal number of clusters without
knowing the true cluster assignments.

First, one splits the data into two parts (A) and (B).
One obtains two cluster assignments, the first one using the centroids
derived from the subset (A), and the second one using the centroids
from the subset (B). Prediction strength measures the proportion of observation
pairs that are assigned to the same clusters according to both clusterings.
The overall prediction strength is the minimum of this quantity over all
predicted clusters.

By varying the desired number of clusters from low to high, we can choose the
highest number of clusters for which the prediction strength exceeds some
threshold. This is precisely how
:class:`sklearn.model_selection.PredictionStrengthGridSearchCV` operates,
as illustrated in the example below. We evaluate ``n_clusters`` in the range
2 to 8 via 5-fold cross-validation. While the average prediction strength
is high for 2, 3, and 4, it sharply drops below the threshold of 0.8 if
``n_clusters`` is 5 or higher. Therefore, we can conclude that the optimal
number of clusters is 4.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import sem

from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.model_selection import PredictionStrengthGridSearchCV
from sklearn.model_selection import KFold

# Generating the sample data from make_blobs
# This particular setting has one distinct cluster and 3 clusters placed close
# together.
X, y = make_blobs(n_samples=500,
                  n_features=2,
                  centers=4,
                  cluster_std=1,
                  center_box=(-10.0, 10.0),
                  shuffle=True,
                  random_state=1)  # For reproducibility

# Define list of values for n_clusters we want to explore
range_n_clusters = [2, 3, 4, 5, 6, 7, 8]
param_grid = {'n_clusters': range_n_clusters}

# Determine optimal choice of n_clusters using 5-fold cross-validation.
# The optimal number of clusters k is the largest k such that the
# corresponding prediction strength is above some threshold.
# Tibshirani and Guenther suggest a threshold in the range 0.8 to 0.9
# for well separated clusters.
clusterer = KMeans(random_state=10)
n_splits = 5
grid_search = PredictionStrengthGridSearchCV(clusterer, threshold=0.8,
                                             param_grid=param_grid,
                                             cv=KFold(n_splits))
grid_search.fit(X)

# Retrieve the best configuration
print(grid_search.best_params_, grid_search.best_score_)

# Retrieve the results stored in the cv_results_ attribute
n_parameters = len(range_n_clusters)
param_n_clusters = grid_search.cv_results_["param_n_clusters"]
mean_test_score = grid_search.cv_results_["mean_test_score"]

# plot average prediction strength for each value for n_clusters
points = np.empty((n_parameters, 2), dtype=np.float_)
for i, values in enumerate(zip(param_n_clusters, mean_test_score)):
    points[i, :] = values
plt.plot(points[:, 0], points[:, 1], marker='o', markerfacecolor='none')
plt.xlabel("n_clusters")
plt.ylabel("average prediction strength")

# plot the standard error of the prediction strength as error bars
test_score_keys = ["split%d_test_score" % split_i
                   for split_i in range(n_splits)]
test_scores = [grid_search.cv_results_[key] for key in test_score_keys]
se = np.fromiter((sem(values) for values in zip(*test_scores)),
                 dtype=np.float_)
plt.errorbar(points[:, 0], points[:, 1], se)

plt.hlines(grid_search.threshold, min(range_n_clusters), max(range_n_clusters),
           linestyles='dashed')

plt.show()
