"""
====================================================================
Selection of the optimal number of clusters based on agreement score
====================================================================

V-Measure is a score that grows when two cluster labeling assignments
agree on how to split a dataset, independently of the number of
splits. This measure is furthermore symmetric::

    v_measure_score(label_a, label_b) == v_measure_score(label_b, label_a)

The following demonstrates how it is possible to randomly split a dataset
into two overlapping subsets, perform clustering on each split (here
using k-means) and then measure their agreement on the intersection of
the two subsets.

This process can be repeated for various number of the parameter k
of k-means. The best V-Measure score can hence give us a hint on the
optimal value for the parameter k in a completely unsupervied fashion.

In order to find a robust estimate of the optimal agreement score, the
process can be repeated and averaged several times by boostrapping the
overlapping sets used for clustering.
"""
# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD
print __doc__

import numpy as np
from scikits.learn.utils import shuffle
from scikits.learn.utils import check_random_state
from scikits.learn.cross_val import Bootstrap
from scikits.learn.grid_search import GridSearchCV
from scikits.learn.cluster import KMeans
from scikits.learn.metrics import v_measure_score
from scikits.learn.datasets.samples_generator import make_blobs

n_samples = 500
n_features = 2
n_centers = 3
n_groups = 4
cluster_std = 0.3
n_bootstraps = 5
k_range = np.arange(2, 20)

random_state = check_random_state(0)


# make a dataset of 4 top level clusters, all of them being shifted copies of
# the same base cluster, which is turn is composed of smaller clusters
base_cluster, _ = make_blobs(n_samples=n_samples / 4, n_features=n_features,
                             centers=n_centers, cluster_std=cluster_std,
                             center_box=(-8, 8), random_state=0)

# shift the clusters more horizontally than vertically so that there is a
# natural vertical grouping of the clusters
samples_1 = base_cluster.copy()
samples_1[:, 0] -= 10
samples_1[:, 1] -= 3

samples_2 = base_cluster.copy()
samples_2[:, 0] -= 10
samples_2[:, 1] += 3

samples_3 = base_cluster.copy()
samples_3[:, 0] += 10
samples_3[:, 1] -= 3

samples_4 = base_cluster.copy()
samples_4[:, 0] += 10
samples_4[:, 1] += 3

samples = np.concatenate((samples_1, samples_2, samples_3, samples_4))
samples = shuffle(samples, random_state=random_state)

scores = np.zeros((n_bootstraps, len(k_range)))

cv = Bootstrap(n_samples, n_bootstraps=3, random_state=random_state)
gs = GridSearchCV(
    KMeans(init='k-means++', random_state=random_state),
    {'k': k_range},
    score_func=v_measure_score,
    cv=cv,
    verbose=1)
gs.fit(samples)
scores = gs.scores_

# compute a bunch of statistics on the results of the bootstrapped runs so as
# to find the optimal value(s) for k
mean_scores = scores.mean(axis=1)
std_scores = scores.std(axis=1)
best_mean_score = mean_scores.max()

# filter out parameters that have highly varying scores
admissible_scores = std_scores <= std_scores.mean() * 2

# retain only the best scores including a half standard deviation of tolerance
admissible_scores *= mean_scores >= best_mean_score - std_scores / 2

admissible_k = [k for i, k in enumerate(k_range) if admissible_scores[i]]

print "Optimal values for k: " + ", ".join(str(k) for k in admissible_k)

# plot the runs
import pylab as pl
pl.subplot(121)
pl.plot(k_range, mean_scores)
pl.plot(k_range, mean_scores + std_scores / 2, 'b--')
pl.plot(k_range, mean_scores - std_scores / 2, 'b--')
pl.xlim(xmin=k_range.min() - 1, xmax=k_range.max() + 1)
pl.ylim(ymin=0.0, ymax=1.0)
for k in admissible_k:
    pl.vlines(k, 0.0, 1.0, 'g')
pl.title("V-Measure agreement of K-Means on the intersection of two\n"
         "overlapping splits for various values of k")
pl.xlabel("Number of centers 'k' for each run of k-means")
pl.ylabel("V-Measure")

pl.subplot(122)
pl.plot(samples[:, 0], samples[:, 1], '.')
pl.xlim(xmin=-15.0, xmax=15.0)
pl.ylim(ymin=-15.0, ymax=15.0)
pl.show()
