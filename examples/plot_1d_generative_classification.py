"""
Generative Bayesian Classification
==================================
This example shows a 1-dimensional, two-class generative classification
using a Gaussian naive Bayes classifier, and some extensions which drop
the naive Gaussian assumption.

In generative Bayesian classification, each class is separately modeled,
and the class yielding the highest posterior probability is selected in
the classification.
"""

# Author: Jake Vanderplas <jakevdp@cs.washington.edu>
# License: BSD 3 Clause

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from sklearn.cross_validation import cross_val_score
from sklearn.naive_bayes import GenerativeBayes
from sklearn.neighbors.kde import KernelDensity
from sklearn.mixture import GMM

# Generate some two-class data with slight overlap
np.random.seed(0)
X1 = np.vstack([stats.laplace.rvs(2.0, 1, size=(1000, 1)),
                stats.laplace.rvs(0.3, 0.2, size=(300,1))])
X2 = np.vstack([stats.laplace.rvs(-2.5, 1, size=(300, 1)),
                stats.laplace.rvs(-1.0, 0.5, size=(200, 1))])
X = np.vstack([X1, X2])
y = np.hstack([np.ones(X1.size), np.zeros(X2.size)])
x_plot = np.linspace(-6, 6, 200)

# Test three density estimators
density_estimators = ['normal_approximation',
                      GMM(3),
                      KernelDensity(0.25)]
names = ['Normal Approximation',
         'Gaussian Mixture Model',
         'Kernel Density Estimation']
linestyles = [':', '--', '-']
colors = []

# Plot histograms of the two input distributions
fig, ax = plt.subplots()
for j in range(2):
    h = ax.hist(X[y == j, 0], bins=np.linspace(-6, 6, 80),
                histtype='stepfilled', normed=False,
                alpha=0.3)
    colors.append(h[2][0].get_facecolor())
binsize = h[1][1] - h[1][0]


for i in range(3):
    clf = GenerativeBayes(density_estimator=density_estimators[i])
    clf.fit(X, y)
    L = np.exp(clf._joint_log_likelihood(x_plot[:, None]))

    for j in range(2):
        ax.plot(x_plot,
                L[:, j] * np.sum(y == j) * binsize / clf.class_prior_[j],
                linestyle=linestyles[i],
                color=colors[j],
                alpha=1)

    # Trick the legend into showing what we want
    scores = cross_val_score(clf, X, y, scoring="accuracy", cv=10)
    ax.plot([], [], linestyle=linestyles[i], color='black',
            label="{0}:\n  {1:.1f}% accuracy.".format(names[i],
                                                      100 * scores.mean()))

ax.set_xlabel('$x$')
ax.set_ylabel('$N(x)$')
ax.legend(loc='upper left', fontsize=12)

plt.show()
