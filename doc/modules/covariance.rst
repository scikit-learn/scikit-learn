.. _covariance:

===================================================
Covariance estimation
===================================================

.. currentmodule:: sklearn.covariance


Many statistical problems require at some point the estimation of a
population's covariance matrix, which can be seen as an estimation of
data set scatter plot shape. Most of the time, such an estimation has
to be done on a sample whose properties (size, structure, homogeneity)
has a large influence on the estimation's quality. The
`sklearn.covariance` package aims at providing tools affording
an accurate estimation of a population's covariance matrix under
various settings.

We assume that the observations are independent and identically
distributed (i.i.d.).


Empirical covariance
====================

The covariance matrix of a data set is known to be well approximated
with the classical *maximum likelihood estimator* (or "empirical
covariance"), provided the number of observations is large enough
compared to the number of features (the variables describing the
observations). More precisely, the Maximum Likelihood Estimator of a
sample is an unbiased estimator of the corresponding population
covariance matrix.

The empirical covariance matrix of a sample can be computed using the
:func:`empirical_covariance` function of the package, or by fitting an
:class:`EmpiricalCovariance` object to the data sample with the
:meth:`EmpiricalCovariance.fit` method.  Be careful that depending
whether the data are centered or not, the result will be different, so
one may want to use the ``assume_centered`` parameter accurately. More precisely
if one uses ``assume_centered=False``, then the test set is supposed to have the
same mean vector as the training set. If not so, both should be centered by the 
user, and ``assume_centered=True`` should be used.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_covariance_estimation.py` for
     an example on how to fit an :class:`EmpiricalCovariance` object
     to data.


.. _shrunk_covariance:

Shrunk Covariance
=================

Basic shrinkage
---------------

Despite being an unbiased estimator of the covariance matrix, the
Maximum Likelihood Estimator is not a good estimator of the
eigenvalues of the covariance matrix, so the precision matrix obtained
from its inversion is not accurate. Sometimes, it even occurs that the
empirical covariance matrix cannot be inverted for numerical
reasons. To avoid such an inversion problem, a transformation of the
empirical covariance matrix has been introduced: the ``shrinkage``.

In the scikit-learn, this transformation (with a user-defined shrinkage
coefficient) can be directly applied to a pre-computed covariance with
the :func:`shrunk_covariance` method. Also, a shrunk estimator of the
covariance can be fitted to data with a :class:`ShrunkCovariance` object
and its :meth:`ShrunkCovariance.fit` method.  Again, depending whether
the data are centered or not, the result will be different, so one may
want to use the ``assume_centered`` parameter accurately.


Mathematically, this shrinkage consists in reducing the ratio between the
smallest and the largest eigenvalue of the empirical covariance matrix.
It can be done by simply shifting every eigenvalue according to a given
offset, which is equivalent of finding the l2-penalized Maximum
Likelihood Estimator of the covariance matrix. In practice, shrinkage
boils down to a simple a convex transformation : :math:`\Sigma_{\rm
shrunk} = (1-\alpha)\hat{\Sigma} + \alpha\frac{{\rm
Tr}\hat{\Sigma}}{p}\rm Id`.

Choosing the amount of shrinkage, :math:`\alpha` amounts to setting a
bias/variance trade-off, and is discussed below.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_covariance_estimation.py` for
     an example on how to fit a :class:`ShrunkCovariance` object
     to data.


Ledoit-Wolf shrinkage
---------------------

In their 2004 paper [1], O. Ledoit and M. Wolf propose a formula so as
to compute the optimal shrinkage coefficient :math:`\alpha` that
minimizes the Mean Squared Error between the estimated and the real
covariance matrix.

The Ledoit-Wolf estimator of the covariance matrix can be computed on
a sample with the :meth:`ledoit_wolf` function of the
`sklearn.covariance` package, or it can be otherwise obtained by
fitting a :class:`LedoitWolf` object to the same sample.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_covariance_estimation.py` for
     an example on how to fit a :class:`LedoitWolf` object to data and
     for visualizing the performances of the Ledoit-Wolf estimator in
     terms of likelihood.


[1] O. Ledoit and M. Wolf, "A Well-Conditioned Estimator for Large-Dimensional
    Covariance Matrices", Journal of Multivariate Analysis, Volume 88, Issue 2,
    February 2004, pages 365-411.

.. _oracle_approximating_shrinkage:

Oracle Approximating Shrinkage
------------------------------

Under the assumption that the data are Gaussian distributed, Chen et
al. [2] derived a formula aimed at choosing a shrinkage coefficient that
yields a smaller Mean Squared Error than the one given by Ledoit and
Wolf's formula. The resulting estimator is known as the Oracle
Shrinkage Approximating estimator of the covariance.

The OAS estimator of the covariance matrix can be computed on a sample
with the :meth:`oas` function of the `sklearn.covariance`
package, or it can be otherwise obtained by fitting an :class:`OAS`
object to the same sample.

.. figure:: ../auto_examples/covariance/images/plot_covariance_estimation_001.png
   :target: ../auto_examples/covariance/plot_covariance_estimation.html
   :align: center
   :scale: 65%

   Bias-variance trade-off when setting the shrinkage: comparing the
   choices of Ledoit-Wolf and OAS estimators

[2] Chen et al., "Shrinkage Algorithms for MMSE Covariance Estimation",
    IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_covariance_estimation.py` for
     an example on how to fit an :class:`OAS` object
     to data.

   * See :ref:`example_covariance_plot_lw_vs_oas.py` to visualize the
     Mean Squared Error difference between a :class:`LedoitWolf` and
     an :class:`OAS` estimator of the covariance.


.. figure:: ../auto_examples/covariance/images/plot_lw_vs_oas_001.png
   :target: ../auto_examples/covariance/plot_lw_vs_oas.html
   :align: center
   :scale: 75%


.. _sparse_inverse_covariance:

Sparse inverse covariance
==========================

The matrix inverse of the covariance matrix, often called the precision
matrix, is proportional to the partial correlation matrix. It gives the
partial independence relationship. In other words, if two features are
independent conditionally on the others, the corresponding coefficient in
the precision matrix will be zero. This is why it makes sense to estimate
a sparse precision matrix: by learning independence relations from the
data, the estimation of the covariance matrix is better conditioned. This
is known as *covariance selection*.

In the small-samples situation, in which ``n_samples`` is on the order
of ``n_features`` or smaller, sparse inverse covariance estimators tend to work
better than shrunk covariance estimators. However, in the opposite
situation, or for very correlated data, they can be numerically unstable.
In addition, unlike shrinkage estimators, sparse estimators are able to
recover off-diagonal structure.

The :class:`GraphLasso` estimator uses an l1 penalty to enforce sparsity on
the precision matrix: the higher its ``alpha`` parameter, the more sparse
the precision matrix. The corresponding :class:`GraphLassoCV` object uses
cross-validation to automatically set the ``alpha`` parameter.

.. figure:: ../auto_examples/covariance/images/plot_sparse_cov_001.png
   :target: ../auto_examples/covariance/plot_sparse_cov.html
   :align: center
   :scale: 60%

   *A comparison of maximum likelihood, shrinkage and sparse estimates of
   the covariance and precision matrix in the very small samples
   settings.*

.. note:: **Structure recovery**

   Recovering a graphical structure from correlations in the data is a
   challenging thing. If you are interested in such recovery keep in mind
   that:

   * Recovery is easier from a correlation matrix than a covariance
     matrix: standardize your observations before running :class:`GraphLasso`

   * If the underlying graph has nodes with much more connections than
     the average node, the algorithm will miss some of these connections.

   * If your number of observations is not large compared to the number
     of edges in your underlying graph, you will not recover it.

   * Even if you are in favorable recovery conditions, the alpha
     parameter chosen by cross-validation (e.g. using the
     :class:`GraphLassoCV` object) will lead to selecting too many edges.
     However, the relevant edges will have heavier weights than the
     irrelevant ones.

The mathematical formulation is the following:

.. math::

    \hat{K} = \mathrm{argmin}_K \big(
                \mathrm{tr} S K - \mathrm{log} \mathrm{det} K
                + \alpha \|K\|_1
                \big)

Where :math:`K` is the precision matrix to be estimated, and :math:`S` is the
sample covariance matrix. :math:`\|K\|_1` is the sum of the absolute values of
off-diagonal coefficients of :math:`K`. The algorithm employed to solve this
problem is the GLasso algorithm, from the Friedman 2008 Biostatistics
paper. It is the same algorithm as in the R ``glasso`` package.


.. topic:: Examples:

   * :ref:`example_covariance_plot_sparse_cov.py`: example on synthetic
     data showing some recovery of a structure, and comparing to other
     covariance estimators.

   * :ref:`example_applications_plot_stock_market.py`: example on real
     stock market data, finding which symbols are most linked.

.. topic:: References:

   * Friedman et al, `"Sparse inverse covariance estimation with the
     graphical lasso" <http://biostatistics.oxfordjournals.org/content/9/3/432.short>`_,
     Biostatistics 9, pp 432, 2008

.. _robust_covariance:

Robust Covariance Estimation
============================

Real data set are often subjects to measurement or recording
errors. Regular but uncommon observations may also appear for a variety
of reason. Every observation which is very uncommon is called an
outlier.
The empirical covariance estimator and the shrunk covariance
estimators presented above are very sensitive to the presence of
outlying observations in the data. Therefore, one should use robust
covariance estimators to estimate the covariance of its real data
sets. Alternatively, robust covariance estimators can be used to
perform outlier detection and discard/downweight some observations
according to further processing of the data.

The ``sklearn.covariance`` package implements a robust estimator of covariance,
the Minimum Covariance Determinant [3].


Minimum Covariance Determinant
------------------------------

The Minimum Covariance Determinant estimator is a robust estimator of
a data set's covariance introduced by P.J. Rousseeuw in [3].  The idea
is to find a given proportion (h) of "good" observations which are not
outliers and compute their empirical covariance matrix.  This
empirical covariance matrix is then rescaled to compensate the
performed selection of observations ("consistency step").  Having
computed the Minimum Covariance Determinant estimator, one can give
weights to observations according to their Mahalanobis distance,
leading to a reweighted estimate of the covariance matrix of the data
set ("reweighting step").

Rousseeuw and Van Driessen [4] developed the FastMCD algorithm in order
to compute the Minimum Covariance Determinant. This algorithm is used
in scikit-learn when fitting an MCD object to data. The FastMCD
algorithm also computes a robust estimate of the data set location at
the same time.

Raw estimates can be accessed as ``raw_location_`` and ``raw_covariance_``
attributes of a :class:`MinCovDet` robust covariance estimator object.

[3] P. J. Rousseeuw. Least median of squares regression.
    J. Am Stat Ass, 79:871, 1984.
[4] A Fast Algorithm for the Minimum Covariance Determinant Estimator,
    1999, American Statistical Association and the American Society
    for Quality, TECHNOMETRICS.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_robust_vs_empirical_covariance.py` for
     an example on how to fit a :class:`MinCovDet` object to data and see how
     the estimate remains accurate despite the presence of outliers.

   * See :ref:`example_covariance_plot_mahalanobis_distances.py` to
     visualize the difference between :class:`EmpiricalCovariance` and
     :class:`MinCovDet` covariance estimators in terms of Mahalanobis distance
     (so we get a better estimate of the precision matrix too).

.. |robust_vs_emp| image:: ../auto_examples/covariance/images/plot_robust_vs_empirical_covariance_001.png
   :target: ../auto_examples/covariance/plot_robust_vs_empirical_covariance.html
   :scale: 49%

.. |mahalanobis| image:: ../auto_examples/covariance/images/plot_mahalanobis_distances_001.png
   :target: ../auto_examples/covariance/plot_mahalanobis_distances.html
   :scale: 49%



____

.. list-table::
    :header-rows: 1

    * - Influence of outliers on location and covariance estimates
      - Separating inliers from outliers using a Mahalanobis distance

    * - |robust_vs_emp|
      - |mahalanobis|

