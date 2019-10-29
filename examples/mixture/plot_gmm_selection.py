"""
================================
Gaussian Mixture Model Selection
================================

This example shows that model selection can be performed with
Gaussian Mixture Models using information-theoretic criteria (BIC).
Model selection concerns both the covariance type
and the number of components in the model.
In that case, AIC also provides the right result (not shown to save time),
but BIC is better suited if the problem is to identify the right model.
Unlike Bayesian procedures, such inferences are prior-free.

In this case, the model with 3 components and full covariance is selected.
"""

import numpy as np
import itertools

from scipy import linalg
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn import mixture
from sklearn.datasets import make_blobs

print(__doc__)

n_samples = 500
skew_points = True
n_components = 3

# Generate random sample, two components
np.random.seed(1)
n_features = 2
cluster_std = np.random.random(n_components) + .5
X = make_blobs(n_samples=n_samples,
               n_features=n_features,
               centers=n_components,
               cluster_std=cluster_std)[0]
if skew_points:
    transformation = [[1, .8], [-.2, 1]]
    X = np.dot(X, transformation)

lowest_bic = np.infty
bic = []
n_components_range = range(1, 7)
cv_types = ['spherical', 'tied', 'diag', 'full']
for cv_type in cv_types:
    for n_components_trial in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = mixture.GaussianMixture(n_components=n_components_trial,
                                      covariance_type=cv_type)
        gmm.fit(X)
        bic.append(gmm.bic(X))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            best_gmm = gmm

bic = np.array(bic)
color_iter = itertools.cycle(['navy', 'turquoise', 'cornflowerblue',
                              'darkorange'])
clf = best_gmm
bars = []

# Plot the BIC scores
plt.figure(figsize=(8, 6))
spl = plt.subplot(2, 1, 1)
for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
    xpos = np.array(n_components_range) + .2 * (i - 2)
    bars.append(plt.bar(xpos, bic[i * len(n_components_range):
                                  (i + 1) * len(n_components_range)],
                        width=.2, color=color))
plt.xticks(n_components_range)
plt.ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
plt.title('BIC score per model')
xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
    .2 * np.floor(bic.argmin() / len(n_components_range))
plt.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
spl.set_xlabel('Number of components')
spl.legend([b[0] for b in bars], cv_types)

# Plot the winner
splot = plt.subplot(2, 1, 2)
Y_ = clf.predict(X)

# convert covariances_ attribute to shape:
#  (n_components, n_features, n_features)

if clf.covariance_type == 'spherical':
    identity_mx = np.identity(n_features)
    covariances = identity_mx[:, :, np.newaxis] * clf.covariances_
    covariances = covariances.T
elif clf.covariance_type == 'tied':
    covariances = clf.covariances_[np.newaxis, :, :]
    covariances = np.repeat(covariances, clf.n_components, axis=0)
elif clf.covariance_type == 'diag':
    identity_mx = np.identity(n_features)
    covariances = np.repeat(identity_mx[:, :, np.newaxis],
                            clf.n_components,
                            axis=2)
    covariances = covariances.T * clf.covariances_[:, :, np.newaxis]
elif clf.covariance_type == 'full':
    covariances = clf.covariances_

for i, (mean, cov, color) in enumerate(zip(clf.means_, covariances,
                                           color_iter)):
    v, w = linalg.eigh(cov)
    if not np.any(Y_ == i):
        continue
    plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)

    # Plot an ellipse to show the Gaussian component
    angle = np.arctan2(w[0][1], w[0][0])
    angle = 180. * angle / np.pi  # convert to degrees
    v = 2. * np.sqrt(2.) * np.sqrt(v)
    ell = mpl.patches.Ellipse(mean, v[0], v[1], 180. + angle, color=color)
    ell.set_clip_box(splot.bbox)
    ell.set_alpha(.5)
    splot.add_artist(ell)

plt.xticks(())
plt.yticks(())
plt.title('Selected GMM: {cv_type} model, {n_components} components'.format(
    cv_type=clf.covariance_type,
    n_components=clf.n_components))
plt.subplots_adjust(hspace=.35, bottom=.02)
plt.show()
