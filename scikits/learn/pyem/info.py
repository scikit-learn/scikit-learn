"""
Routines for Gaussian Mixture Models
and learning with Expectation Maximization 
==========================================

This module contains classes and function to compute multivariate Gaussian densities
(diagonal and full covariance matrices), Gaussian mixtures, Gaussian mixtures models
and an Em trainer.

More specifically, the module contains the following file:

- densities.py: functions to compute multivariate Gaussian pdf and ellipsoides of
confidence (gauss_den)
- gauss_mix.py: defines the GM (Gaussian Mixture) class. A Gaussian Mixture can be
created from its parameters weights, mean and variances, or from its meta parameters
d (dimension of the Gaussian) and k (number of components in the mixture). A Gaussian
Model can then be sampled or plot (if d>1, plot confidence ellipsoids projected on 
2 chosen dimensions, if d == 1, plot the pdf of each component and fill the zone
of confidence for a given level)
- gmm_em.py: defines a class GMM (Gaussian Mixture Model). This class is constructed
from a GM model gm, and can be used to train gm. The GMM can be initiated by
kmean or at random, and can compute sufficient statistics, and update its parameters
from the sufficient statistics.
- kmean.py: implements a kmean algorithm. We cannot use scipy.cluster.vq kmeans, since
its does not give membership of observations.

Example of use: simply execute gmm_em.py for an example of training a GM and plotting
the results compared to true values

Copyright: David Cournapeau 2006
License: BSD-style (see LICENSE.txt in main source directory)
"""
version = '0.5.3'

depends = ['linalg', 'stats']
ignore  = False
