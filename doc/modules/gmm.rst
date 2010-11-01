.. _gmm:

===================================================
Gaussian mixture models
===================================================

`scikits.learn.gmm` is a package which enables to create Gaussian
Mixture Models (diagonal, spherical, tied and full covariance matrices
supported), to sample them, and to estimate them from data using
Expectation Maximization algorithm.  It can also draw confidence
ellipsoides for multivariate models, and compute the Bayesian
Information Criterion to assess the number of clusters in the data.
    

GMM classifier
==============

.. autoclass:: scikits.learn.gmm.GMM
    :members:


.. topic:: Examples:

    * See :ref:`example_gmm_plot_gmm.py` for an example on plotting the
      confidence ellipsoids.

    * See :ref:`example_gmm_plot_gmm_pdf.py` for an example on plotting the 
      density estimation.

