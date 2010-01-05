"""
Routines for Gaussian Mixture Models
and learning with Expectation Maximization 
==========================================

This module contains classes and function to compute multivariate Gaussian densities
(diagonal and full covariance matrices), Gaussian mixtures, Gaussian mixtures models
and an Em trainer.

More specifically, the module defines the following classes, functions:

- densities.gauss_den: function to compute multivariate Gaussian pdf 
- gauss_mix.GM: defines the GM (Gaussian Mixture) class. A Gaussian Mixture can be
created from its parameters weights, mean and variances, or from its meta parameters
d (dimension of the Gaussian) and k (number of components in the mixture). A Gaussian
Model can then be sampled or plot (if d>1, plot confidence ellipsoids projected on 
2 chosen dimensions, if d == 1, plot the pdf of each component and fill the zone
of confidence for a given level)
- gmm_em.GMM: defines a class GMM (Gaussian Mixture Model). This class is constructed
from a GM model gm, and can be used to train gm. The GMM can be initiated by
kmean or at random, and can compute sufficient statistics, and update its parameters
from the sufficient statistics.
- kmean.kmean: implements a kmean algorithm. We cannot use scipy.cluster.vq kmeans, since
its does not give membership of observations.

Example of use: 
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Create an artificial 2 dimension, 3 clusters GM model, samples it
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    w, mu, va   = GM.gen_param(2, 3, 'diag', spread = 1.5)
    gm          = GM.fromvalues(w, mu, va)

    # Sample 1000 frames  from the model
    data    = gm.sample(1000)

    #++++++++++++++++++++++++
    # Learn the model with EM
    #++++++++++++++++++++++++
    # Init the model
    lgm = GM(d, k, mode)
    gmm = GMM(lgm, 'kmean')

    # The actual EM, with likelihood computation. The threshold
    # is compared to the (linearly appromixated) derivative of the likelihood
    em      = EM()
    like    = em.train(data, gmm, maxiter = 30, thresh = 1e-8)

Files example.py and example2.py show more capabilities of the toolbox, including
plotting capabilities (using matplotlib) and model selection using Bayesian 
Information Criterion (BIC).

Bibliography:
    * Maximum likelihood from incomplete data via the EM algorithm in Journal of 
    the Royal Statistical Society, Series B, 39(1):1--38, 1977, by A. P. Dempster, 
    N. M. Laird, and D. B. Rubin
    * Bayesian Approaches to Gaussian Mixture Modelling (1998) by 
    Stephen J. Roberts, Dirk Husmeier, Iead Rezek, William Penny in 
    IEEE Transactions on Pattern Analysis and Machine Intelligence
     
Copyright: David Cournapeau 2006
License: BSD-style (see LICENSE.txt in main source directory)
"""
version = '0.5.6'

depends = ['linalg', 'stats']
ignore  = False
