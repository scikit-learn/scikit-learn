#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=================================================================
Model selection with Probabilistic (PCA) and Factor Analysis (FA)
=================================================================

Probabilistic PCA and Factor Analysis are probabilistic models.
The consequence is that the likelihood of new data can be used
for model selection. Here we compare PCA and FA with cross-validation
on low rank data corrupted with homoscedastic noise (noise variance
is the same for each feature) or heteroscedastic noise (noise variance
is the different for each feature).

One can observe that with homoscedastic noise both FA and PCA succeed
in recovering the size of the low rank subspace. The likelihood with PCA
is higher than FA in this case. However PCA fails and overestimates
the rank when heteroscedastic noise is present. The automatic estimation from
Automatic Choice of Dimensionality for PCA. NIPS 2000: 598-604
by Thomas P. Minka is also compared.

"""
print(__doc__)

# Authors: Alexandre Gramfort
# License: BSD 3 clause

import numpy as np
import pylab as pl
from scipy import linalg

from sklearn.decomposition import PCA, FactorAnalysis
from sklearn.cross_validation import cross_val_score

###############################################################################
# Create the data

n_samples, n_features, rank = 1000, 50, 10
sigma = 1.
rng = np.random.RandomState(42)
U, _, _ = linalg.svd(rng.randn(n_features, n_features))
X = np.dot(rng.randn(n_samples, rank), U[:, :rank].T)

# Adding homoscedastic noise
X_homo = X + sigma * rng.randn(n_samples, n_features)

# Adding heteroscedastic noise
sigmas = sigma * rng.rand(n_features) + sigma / 2.
X_hetero = X + rng.randn(n_samples, n_features) * sigmas

###############################################################################
# Fit the models

n_components = np.arange(0, n_features, 5)  # options for n_components


def compute_scores(X):
    pca = PCA()
    fa = FactorAnalysis()

    pca_scores, fa_scores = [], []
    for n in n_components:
        pca.n_components = n
        fa.n_components = n
        pca_scores.append(np.mean(cross_val_score(pca, X)))
        fa_scores.append(np.mean(cross_val_score(fa, X)))

    return pca_scores, fa_scores

for X, title in [(X_homo, 'Homoscedastic Noise'),
                 (X_hetero, 'Heteroscedastic Noise')]:
    pca_scores, fa_scores = compute_scores(X)
    n_components_pca = n_components[np.argmax(pca_scores)]
    n_components_fa = n_components[np.argmax(fa_scores)]

    pca = PCA(n_components='mle')
    pca.fit(X)
    n_components_pca_mle = pca.n_components_

    print("best n_components by PCA CV = %d" % n_components_pca)
    print("best n_components by FactorAnalysis CV = %d" % n_components_fa)
    print("best n_components by PCA MLE = %d" % n_components_pca_mle)

    pl.figure()
    pl.plot(n_components, pca_scores, 'b', label='PCA scores')
    pl.plot(n_components, fa_scores, 'r', label='FA scores')
    pl.axvline(rank, color='g', label='TRUTH: %d' % rank, linestyle='-')
    pl.axvline(n_components_pca, color='b',
               label='PCA CV: %d' % n_components_pca, linestyle='--')
    pl.axvline(n_components_fa, color='r',
               label='FactorAnalysis CV: %d' % n_components_fa, linestyle='--')
    pl.axvline(n_components_pca_mle, color='k',
               label='PCA MLE: %d' % n_components_pca_mle, linestyle='--')
    pl.xlabel('nb of components')
    pl.ylabel('CV scores')
    pl.legend(loc='lower right')
    pl.title(title)

pl.show()
