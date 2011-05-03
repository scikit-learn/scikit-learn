.. _mixture:

===================================================
Gaussian mixture models
===================================================

`scikits.learn.mixture` is a package which enables to create Mixture
Models (diagonal, spherical, tied and full covariance matrices
supported), to sample them, and to estimate them from data using
Expectation Maximization algorithm.  It can also draw confidence
ellipsoides for multivariate models, and compute the Bayesian
Information Criterion to assess the number of clusters in the data.
    
For the moment, only Gaussian Mixture Models (GMM) are
implemented. These are a class of probabilistic models describing the
data as drawn from a mixture of Gaussian probability
distributions. The challenge that is GMM tackles is to learn the
parameters of these Gaussians from the data.

GMM classifier
==============

.. currentmodule:: scikits.learn.mixture

The :class:`GMM` object implements a :meth:`GMM.fit` method to learn a
Gaussian Mixture Models from train data. Given test data, it can assign
to each sample the class of the Gaussian it mostly probably belong to
using the :meth:`GMM.predict` method. 

..  
    Alternatively, the probability of each
    sample beloning to the various Gaussians may be retrieved using the
    :meth:`GMM.predict_proba` method.

.. figure:: ../auto_examples/mixture/images/plot_gmm_classifier_1.png
   :target: ../auto_examples/cluster/plot_gmm_classifier.html
   :align: center
   :scale: 75%

.. topic:: Examples:

    * See :ref:`example_mixture_plot_gmm_classifier.py` for an example of
      using a GMM as a classifier on the iris dataset.

    * See :ref:`example_mixture_plot_gmm.py` for an example on plotting the
      confidence ellipsoids.

    * See :ref:`example_mixture_plot_gmm_pdf.py` for an example on plotting the 
      density estimation.

