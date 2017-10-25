"""
========================================================================
Comparing different anomaly/novelty detection algorithms on toy datasets
========================================================================

This example shows characteristics of different
anomaly/novelty detection algorithms on 2D datasets. Datasets contain
one or two modes (regions of high density) to illustrate the ability
of algorithms to cope with multimodal data.

While these examples give some intuition about the
algorithms, this intuition might not apply to very high
dimensional data.
"""
print(__doc__)  # noqa

import time

import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt

from sklearn import svm
from sklearn.datasets import make_moons
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

# Example settings
n_samples = 300
outliers_fraction = 0.15
clusters_separation = (0, 1, 2)
n_outliers = int(outliers_fraction * n_samples)
n_inliers = n_samples - n_outliers


def make_dataset(n_inliers, offset, stds=(0.5, 0.5), random_state=0):
    rng = np.random.RandomState(random_state)
    # Data generation
    X1 = stds[0] * rng.randn(n_inliers // 2, 2) - offset
    X2 = stds[1] * rng.randn(n_inliers // 2, 2) + offset
    X = np.concatenate([X1, X2], axis=0)
    return X

# define outlier/anomaly detection methods to be compared
anomaly_algorithms = [
    ("Robust covariance", EllipticEnvelope(contamination=outliers_fraction)),
    ("One-Class SVM", svm.OneClassSVM(nu=0.95 * outliers_fraction + 0.05,
                                      kernel="rbf", gamma=0.1)),
    ("Isolation Forest", IsolationForest(contamination=outliers_fraction,
                                         random_state=42)),
    ("Local Outlier Factor", LocalOutlierFactor(
        n_neighbors=35, contamination=outliers_fraction))]

# Define datasets
datasets = [
    make_dataset(n_inliers, stds=(0.5, 0.5), offset=0.),
    make_dataset(n_inliers, stds=(0.3, 1.5), offset=2.),
    4. * (make_moons(n_samples=n_samples, noise=.05)[0] -
          np.array([0.5, 0.25])),
    14. * (np.random.RandomState(42).rand(n_samples, 2) - 0.5)]

# Compare given classifiers under given settings
xx, yy = np.meshgrid(np.linspace(-7, 7, 150),
                     np.linspace(-7, 7, 150))

plt.figure(figsize=(len(anomaly_algorithms) * 2 + 3, 12.5))
plt.subplots_adjust(left=.02, right=.98, bottom=.001, top=.96, wspace=.05,
                    hspace=.01)

plot_num = 1

for i_dataset, X in enumerate(datasets):
    # Add outliers
    rng = np.random.RandomState(42)
    X = np.concatenate([X, rng.uniform(low=-6, high=6,
                       size=(n_outliers, 2))], axis=0)

    for name, algorithm in anomaly_algorithms:
        t0 = time.time()
        algorithm.fit(X)
        t1 = time.time()
        plt.subplot(len(datasets), len(anomaly_algorithms), plot_num)
        if i_dataset == 0:
            plt.title(name, size=18)

        # fit the data and tag outliers
        if name == "Local Outlier Factor":
            y_pred = algorithm.fit_predict(X)
            anomaly_scores = algorithm.negative_outlier_factor_
        else:
            y_pred = algorithm.fit(X).predict(X)
            anomaly_scores = algorithm.decision_function(X)

        threshold = stats.scoreatpercentile(anomaly_scores,
                                            100 * outliers_fraction)

        # plot the levels lines and the points
        if hasattr(algorithm, 'predict'):
            Z = algorithm.predict(np.c_[xx.ravel(), yy.ravel()])
            Z = Z.reshape(xx.shape)
            plt.contour(xx, yy, Z, levels=[threshold],
                        linewidths=2, colors='black')

        colors = np.array(['#377eb8', '#ff7f00'])
        plt.scatter(X[:, 0], X[:, 1], s=10, color=colors[(y_pred + 1) // 2])

        plt.xlim(-7, 7)
        plt.ylim(-7, 7)
        plt.xticks(())
        plt.yticks(())
        plt.text(.99, .01, ('%.2fs' % (t1 - t0)).lstrip('0'),
                 transform=plt.gca().transAxes, size=15,
                 horizontalalignment='right')
        plot_num += 1

plt.show()
