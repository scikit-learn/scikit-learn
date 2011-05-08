"""
=================================================
Spectral clustering for non convex cluster shapes
=================================================

"""
print __doc__

# Authors:  Olivier Grisel
# License: BSD

from time import time
import numpy as np
import pylab as pl

from scikits.learn.cluster import k_means
from scikits.learn.cluster import affinity_propagation
from scikits.learn.cluster import mean_shift
from scikits.learn.cluster import spectral_clustering
from scikits.learn.cluster import Ward
from scikits.learn.cluster import power_iteration_clustering
from scikits.learn.metrics.pairwise import euclidean_distances
from scikits.learn.metrics import homogeneity_completeness_v_measure
from scikits.learn.neighbors import kneighbors_graph

# Generate random samples roughly arranged as nested circles

circle_parameters = (
    # (center_x, center_y, radius, n_points)
    (0, 0, 10, 100),
    (8, 0, 25, 200),
    (8, 4, 55, 300),
)
noise_level = 0.05
random_state = np.random.RandomState(42)
circles = []

labels = []
for i, (center_x, center_y, radius, n_points) in enumerate(circle_parameters):
    t = random_state.uniform(12 * np.pi, size=n_points)

    circle_x = center_x + radius * np.cos(t)
    circle_y = center_y + radius * np.sin(t)
    circle = np.array([circle_x, circle_y]).T
    noise = random_state.normal(scale=noise_level * radius, size=(n_points, 2))

    circles.append(circle + noise)
    labels += [i] * n_points

X = np.concatenate(circles)
labels_true = np.array(labels)

# Shuffle the samples to ensure that the algo has no way of cheating
indices = np.arange(X.shape[0])
random_state.shuffle(indices)
X = X[indices]
labels_true = labels_true[indices]


# Utility functions to report on the results of the various strategies

def plot_labels(title, labels):
    """Visual clustering port as 2D plot"""
    unique_labels = np.unique(labels)
    for l in unique_labels:
        X_l = X[labels == l, :]
        color = pl.cm.hsv(float(l) / unique_labels.shape[0])
        pl.scatter(X_l[:, 0], X_l[:, 1], color=color)
    pl.title(title)
    pl.xticks(())
    pl.yticks(())


def report(title, labels_true, labels_pred, duration, do_plot=True):
    """Print lustering report on stdout"""
    h, c, v = homogeneity_completeness_v_measure(labels_true, labels_pred)
    print title
    print "Homogeneity: %0.3f" % h
    print "Completeness: %0.3f" % c
    print "V-Measure: %0.3f" % v
    print "Duration: %0.3fs" % duration
    print
    if do_plot:
        title = "%s\nv=%0.2f (%0.3fs)" % (title, v, duration)
        plot_labels(title, labels)

pl.figure()

# Random assignment
t0 = time()
labels = random_state.randint(0, np.unique(labels_true).shape[0],
                              size=labels_true.shape)
duration = time() - t0
pl.subplot(331)
report("Random", labels_true, labels, duration)

# K-Means
t0 = time()
_, labels, inertia = k_means(X, k=3)
duration = time() - t0
pl.subplot(332)
report("K-Means", labels_true, labels, duration)

# Mean Shift
t0 = time()
_, labels = mean_shift(X, bandwidth=28.0)
duration = time() - t0
pl.subplot(333)
report("Mean Shift", labels_true, labels, duration)

# Build a knn graph as affinity matrix
t0 = time()
affinity = kneighbors_graph(X, n_neighbors=10)
affinity = 0.5 * (affinity + affinity.T)  # make affinity symmetric
duration_affinity = time() - t0

# Affinity propagation
# XXX: I cannot get it to work as expected
#_, labels = affinity_propagation(affinity.toarray(), p=0.5)
#pl.subplot(334)
#plot_labels(labels, "Affinity propagation")

# Ward clustering
t0 = time()
labels = Ward(n_clusters=3, connectivity=affinity).fit(X).labels_
duration = time() - t0
pl.subplot(335)
report("Ward", labels_true, labels, duration + duration_affinity)

# Spectral Clustering
# XXX: the spectral clustering results is unstable with the amg-based method
# XXX: we should implement the fast_svd method too
t0 = time()
labels = spectral_clustering(affinity, k=3, mode='arpack',
                             random_state=random_state)
duration = time() - t0
pl.subplot(337)
report("Spectral", labels_true, labels, duration + duration_affinity)

# Power iteration
t0 = time()
labels = power_iteration_clustering(
    affinity, k=3, n_vectors=10, tol=1e-5, random_state=random_state,
    verbose=False, plot_vector=False)
duration = time() - t0
pl.subplot(338)
report("Power Iteration", labels_true, labels, duration + duration_affinity)

pl.show()
