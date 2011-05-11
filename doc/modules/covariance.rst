.. _covariance:

===================================================
Covariance estimation
===================================================

.. currentmodule:: scikits.learn.covariance


Many statistical problems require at some point the estimation of a
population's covariance matrix, which can be seen as an estimation of
data set scatter plot shape. Most of the time, such an estimation has
to be done on a sample whose properties (size, structure, homogeneity)
has a large influence on the estimation's quality. The
`scikits.learn.covariance` package aims at providing tools affording
an accurate estimation of a population's covariance matrix under
various settings.

The package does not include robust tools yet, so we assume that the
data sets do not contain any outlying data. We also assume thah the
observations are independant and identically distributed.

Empirical covariance
====================


The covariance matrix of a data set is known to be well approximated
with the classical `Maximum Likelihood Estimator` (or `empirical
covariance`), provided the number of observations is large enough
compared to the number of features (the variables describing the
observations). More precisely, the Maximum Likelihood Estimator of a
sample is an unbiased estimator of the corresponding population
covariance matrix.

The empirical covariance matrix of a sample can be computed using the
:meth:`empirical_covariance` function of the package, or by fitting an
:class:`EmpiricalCovariance` object to the data sample with the
:meth:`EmpiricalCovariance.fit` method.  Be careful that depending
whether the data are centered or not, the result will be different, so
one may want to use the `assume_centered` parameter accurately.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_covariance_estimation.py` for
     an example on how to fit an :class:`EmpiricalCovariance` object
     to data.

Shrunk Covariance
=================


Basic shrinkage
---------------

Despite it is an unbiased estimator of the covariance matrix, the
Maximum Likelihood Estimator is not a good estimator of the
eigenvalues of the covariance matrix, so the precision matrix obtained
from its inversion is not accurate. Sometimes, it even occurs that the
empirical covariance matrix cannot be inverted for numerical
reasons. To avoid such an inversion problem, a transformation of the
empirical covariance matrix has been introduced: the `shrinkage`. It
consists in reducing the ratio between the smallest and the largest
eigenvalue of the empirical covariance matrix. This can be done by
simply shifting every eigenvalue according to a given offset, which is
equivalent of finding the l2-Penalized Maximum Likelihood Estimator of
the covariance matrix, or by reducing the highest eigenvalue while
increasing the smallest with the help of a convex transformation :
:math:`\Sigma_{\rm shrunk} = (1-\alpha)\hat{\Sigma} +
\alpha\frac{{\rm Tr}\hat{\Sigma}}{p}\rm Id`.  The latter approach has been
implemented in scikit-learn.

A convex transformation (with a user-defined shrinkage coefficient)
can be directly applied to a pre-computed covariance with the
:meth:`shrunk_covariance` method. Also, a shrunk estimator of the
covariance can be fitted to data with a :class:`ShrunkCovariance`
object and its :meth:`ShrunkCovariance.fit` method.  Again, depending
whether the data are centered or not, the result will be different, so
one may want to use the `assume_centered` parameter accurately.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_covariance_estimation.py` for
     an example on how to fit a :class:`ShrunkCovariance` object
     to data.


Ledoit-Wolf shrinkage
---------------------

In their 2004 paper [1], O. Ledoit and M. Wolf propose a formula so as
to compute the optimal shrinkage coefficient :math:`\alpha` that
minimizes the Mean Squared Error between the estimated and the real
covariance matrix in terms of Frobenius norm.

The Ledoit-Wolf estimator of the covariance matrix can be computed on
a sample with the :meth:`ledoit_wolf` function of the
`scikits.learn.covariance` package, or it can be otherwise obtained by
fitting a :class:`LedoitWolf` object to the same sample.

[1] "A Well-Conditioned Estimator for Large-Dimensional Covariance
    Matrices", Ledoit and Wolf, Journal of Multivariate Analysis,
    Volume 88, Issue 2, February 2004, pages 365-411.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_covariance_estimation.py` for
     an example on how to fit a :class:`LedoitWolf` object to data and
     for visualizing the performances of the Ledoit-Wolf estimator in
     terms of likelihood.

.. figure:: ../auto_examples/covariance/images/plot_covariance_estimation_-1.png
   :target: ../auto_examples/covariance/plot_covariance_estimation.html
   :align: center
   :scale: 75%


.. _oracle_apprroximating_shrinkage:

Oracle Approximating Shrinkage
------------------------------

Under the assumption that the data are Gaussian distributed, Chen et
al. [2] derived a formula aimed at choosing a shrinkage coefficient that
yields a smaller Mean Squared Error than the one given by Ledoit and
Wolf's formula. The resulting estimator is known as the Oracle
Shrinkage Approximating estimator of the covariance.

The OAS estimator of the covariance matrix can be computed on a sample
with the :meth:`oas` function of the `scikits.learn.covariance`
package, or it can be otherwise obtained by fitting an :class:`OAS`
object to the same sample.  The formula we used to implement the OAS
does not correspond to the one given in the article. It has been taken
from the matlab programm available from the authors webpage
(https://tbayes.eecs.umich.edu/yilun/covestimation).


[2] "Shrinkage Algorithms for MMSE Covariance Estimation" Chen et al.,
    IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_covariance_estimation.py` for
     an example on how to fit an :class:`OAS` object
     to data.

   * See :ref:`example_covariance_plot_lw_vs_oas.py` to visualize the
     Mean Squared Error difference between a :class:`LedoitWolf` and
     an :class:`OAS` estimator of the covariance.


.. figure:: ../auto_examples/covariance/images/plot_lw_vs_oas_1.png
   :target: ../auto_examples/covariance/plot_lw_vs_oas.html
   :align: center
   :scale: 75%
