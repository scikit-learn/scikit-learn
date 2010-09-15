.. _gmm:

===================================================
Gaussian mixture models
===================================================

.. contents:: Tables of contents

`scikits.learn.gmm` is a package which enables to create Gaussian
Mixture Models (diagonal, spherical, tied and full covariance matrices
supported), to sample them, and to estimate them from data using
Expectation Maximization algorithm.  It can also draw confidence
ellipsoides for multivariate models, and compute the Bayesian
Information Criterion to assess the number of clusters in the data. In
a near future, I hope to add so-called online EM (ie recursive EM) and
variational Bayes implementation.

It is implemented in python, and uses the excellent numpy and scipy
packages. Numpy is a python packages which gives python a fast 
multi-dimensional array capabilities (ala matlab and the likes); scipy
leverages numpy to build common scientific features for signal processing,
linear algebra, statistics, etc...
     

GMM classifier
==============

.. autoclass:: scikits.learn.gmm.GMM
    :members:


Examples
--------

See :ref:`example_gmm_plot_gmm.py` for an example on plotting the
confidence ellipsoids.

See :ref:`example_gmm_plot_gmm_pdf.py` for an example on plotting the density
estimation.

