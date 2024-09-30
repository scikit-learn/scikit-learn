.. _lda_qda:

==========================================
Linear and Quadratic Discriminant Analysis
==========================================

.. currentmodule:: sklearn

Linear Discriminant Analysis
(:class:`~discriminant_analysis.LinearDiscriminantAnalysis`) and Quadratic
Discriminant Analysis
(:class:`~discriminant_analysis.QuadraticDiscriminantAnalysis`) are two classic
classifiers, with, as their names suggest, a linear and a quadratic decision
surface, respectively.

These classifiers are attractive because they have closed-form solutions that
can be easily computed, are inherently multiclass, have proven to work well in
practice, and have no hyperparameters to tune.

.. |ldaqda| image:: ../auto_examples/classification/images/sphx_glr_plot_lda_qda_001.png
        :target: ../auto_examples/classification/plot_lda_qda.html
        :scale: 80

.. centered:: |ldaqda|

The plot shows decision boundaries for Linear Discriminant Analysis and
Quadratic Discriminant Analysis. The bottom row demonstrates that Linear
Discriminant Analysis can only learn linear boundaries, while Quadratic
Discriminant Analysis can learn quadratic boundaries and is therefore more
flexible.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_classification_plot_lda_qda.py`: Comparison of LDA and
  QDA on synthetic data.

Dimensionality reduction using Linear Discriminant Analysis
===========================================================

:class:`~discriminant_analysis.LinearDiscriminantAnalysis` can be used to
perform supervised dimensionality reduction, by projecting the input data to a
linear subspace consisting of the directions which maximize the separation
between classes (in a precise sense discussed in the mathematics section
below). The dimension of the output is necessarily less than the number of
classes, so this is in general a rather strong dimensionality reduction, and
only makes sense in a multiclass setting.

This is implemented in the `transform` method. The desired dimensionality can
be set using the ``n_components`` parameter. This parameter has no influence
on the `fit` and `predict` methods.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_decomposition_plot_pca_vs_lda.py`: Comparison of LDA and
  PCA for dimensionality reduction of the Iris dataset

.. _lda_qda_math:

Mathematical formulation of the LDA and QDA classifiers
=======================================================

Both LDA and QDA can be derived from simple probabilistic models which model
the class conditional distribution of the data :math:`P(X|y=k)` for each class
:math:`k`. Predictions can then be obtained by using Bayes' rule, for each
training sample :math:`x \in \mathcal{R}^d`:

.. math::
    P(y=k | x) = \frac{P(x | y=k) P(y=k)}{P(x)} = \frac{P(x | y=k) P(y = k)}{ \sum_{l} P(x | y=l) \cdot P(y=l)}

and we select the class :math:`k` which maximizes this posterior probability.

More specifically, for linear and quadratic discriminant analysis,
:math:`P(x|y)` is modeled as a multivariate Gaussian distribution with
density:

.. math:: P(x | y=k) = \frac{1}{(2\pi)^{d/2} |\Sigma_k|^{1/2}}\exp\left(-\frac{1}{2} (x-\mu_k)^t \Sigma_k^{-1} (x-\mu_k)\right)

where :math:`d` is the number of features.

QDA
---

According to the model above, the log of the posterior is:

.. math::

    \log P(y=k | x) &= \log P(x | y=k) + \log P(y = k) + Cst \\
    &= -\frac{1}{2} \log |\Sigma_k| -\frac{1}{2} (x-\mu_k)^t \Sigma_k^{-1} (x-\mu_k) + \log P(y = k) + Cst,

where the constant term :math:`Cst` corresponds to the denominator
:math:`P(x)`, in addition to other constant terms from the Gaussian. The
predicted class is the one that maximises this log-posterior.

.. note:: **Relation with Gaussian Naive Bayes**

	  If in the QDA model one assumes that the covariance matrices are diagonal,
	  then the inputs are assumed to be conditionally independent in each class,
	  and the resulting classifier is equivalent to the Gaussian Naive Bayes
	  classifier :class:`naive_bayes.GaussianNB`.

LDA
---

LDA is a special case of QDA, where the Gaussians for each class are assumed
to share the same covariance matrix: :math:`\Sigma_k = \Sigma` for all
:math:`k`. This reduces the log posterior to:

.. math:: \log P(y=k | x) = -\frac{1}{2} (x-\mu_k)^t \Sigma^{-1} (x-\mu_k) + \log P(y = k) + Cst.

The term :math:`(x-\mu_k)^t \Sigma^{-1} (x-\mu_k)` corresponds to the
`Mahalanobis Distance <https://en.wikipedia.org/wiki/Mahalanobis_distance>`_
between the sample :math:`x` and the mean :math:`\mu_k`. The Mahalanobis
distance tells how close :math:`x` is from :math:`\mu_k`, while also
accounting for the variance of each feature. We can thus interpret LDA as
assigning :math:`x` to the class whose mean is the closest in terms of
Mahalanobis distance, while also accounting for the class prior
probabilities.

The log-posterior of LDA can also be written [3]_ as:

.. math::

    \log P(y=k | x) = \omega_k^t x + \omega_{k0} + Cst.

where :math:`\omega_k = \Sigma^{-1} \mu_k` and :math:`\omega_{k0} =
-\frac{1}{2} \mu_k^t\Sigma^{-1}\mu_k + \log P (y = k)`. These quantities
correspond to the `coef_` and `intercept_` attributes, respectively.

From the above formula, it is clear that LDA has a linear decision surface.
In the case of QDA, there are no assumptions on the covariance matrices
:math:`\Sigma_k` of the Gaussians, leading to quadratic decision surfaces.
See [1]_ for more details.

Mathematical formulation of LDA dimensionality reduction
========================================================

First note that the K means :math:`\mu_k` are vectors in
:math:`\mathcal{R}^d`, and they lie in an affine subspace :math:`H` of
dimension at most :math:`K - 1` (2 points lie on a line, 3 points lie on a
plane, etc.).

As mentioned above, we can interpret LDA as assigning :math:`x` to the class
whose mean :math:`\mu_k` is the closest in terms of Mahalanobis distance,
while also accounting for the class prior probabilities. Alternatively, LDA
is equivalent to first *sphering* the data so that the covariance matrix is
the identity, and then assigning :math:`x` to the closest mean in terms of
Euclidean distance (still accounting for the class priors).

Computing Euclidean distances in this d-dimensional space is equivalent to
first projecting the data points into :math:`H`, and computing the distances
there (since the other dimensions will contribute equally to each class in
terms of distance). In other words, if :math:`x` is closest to :math:`\mu_k`
in the original space, it will also be the case in :math:`H`.
This shows that, implicit in the LDA
classifier, there is a dimensionality reduction by linear projection onto a
:math:`K-1` dimensional space.

We can reduce the dimension even more, to a chosen :math:`L`, by projecting
onto the linear subspace :math:`H_L` which maximizes the variance of the
:math:`\mu^*_k` after projection (in effect, we are doing a form of PCA for the
transformed class means :math:`\mu^*_k`). This :math:`L` corresponds to the
``n_components`` parameter used in the
:func:`~discriminant_analysis.LinearDiscriminantAnalysis.transform` method. See
[1]_ for more details.

Shrinkage and Covariance Estimator
==================================

Shrinkage is a form of regularization used to improve the estimation of
covariance matrices in situations where the number of training samples is
small compared to the number of features.
In this scenario, the empirical sample covariance is a poor
estimator, and shrinkage helps improving the generalization performance of
the classifier.
Shrinkage LDA can be used by setting the ``shrinkage`` parameter of
the :class:`~discriminant_analysis.LinearDiscriminantAnalysis` class to 'auto'.
This automatically determines the optimal shrinkage parameter in an analytic
way following the lemma introduced by Ledoit and Wolf [2]_. Note that
currently shrinkage only works when setting the ``solver`` parameter to 'lsqr'
or 'eigen'.

The ``shrinkage`` parameter can also be manually set between 0 and 1. In
particular, a value of 0 corresponds to no shrinkage (which means the empirical
covariance matrix will be used) and a value of 1 corresponds to complete
shrinkage (which means that the diagonal matrix of variances will be used as
an estimate for the covariance matrix). Setting this parameter to a value
between these two extrema will estimate a shrunk version of the covariance
matrix.

The shrunk Ledoit and Wolf estimator of covariance may not always be the
best choice. For example if the distribution of the data
is normally distributed, the
Oracle Approximating Shrinkage estimator :class:`sklearn.covariance.OAS`
yields a smaller Mean Squared Error than the one given by Ledoit and Wolf's
formula used with shrinkage="auto". In LDA, the data are assumed to be gaussian
conditionally to the class. If these assumptions hold, using LDA with
the OAS estimator of covariance will yield a better classification
accuracy than if Ledoit and Wolf or the empirical covariance estimator is used.

The covariance estimator can be chosen using with the ``covariance_estimator``
parameter of the :class:`discriminant_analysis.LinearDiscriminantAnalysis`
class. A covariance estimator should have a :term:`fit` method and a
``covariance_`` attribute like all covariance estimators in the
:mod:`sklearn.covariance` module.


.. |shrinkage| image:: ../auto_examples/classification/images/sphx_glr_plot_lda_001.png
        :target: ../auto_examples/classification/plot_lda.html
        :scale: 75

.. centered:: |shrinkage|

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_classification_plot_lda.py`: Comparison of LDA classifiers
  with Empirical, Ledoit Wolf and OAS covariance estimator.

Estimation algorithms
=====================

Using LDA and QDA requires computing the log-posterior which depends on the
class priors :math:`P(y=k)`, the class means :math:`\mu_k`, and the
covariance matrices.

The 'svd' solver is the default solver used for
:class:`~sklearn.discriminant_analysis.LinearDiscriminantAnalysis`, and it is
the only available solver for
:class:`~sklearn.discriminant_analysis.QuadraticDiscriminantAnalysis`.
It can perform both classification and transform (for LDA).
As it does not rely on the calculation of the covariance matrix, the 'svd'
solver may be preferable in situations where the number of features is large.
The 'svd' solver cannot be used with shrinkage.
For QDA, the use of the SVD solver relies on the fact that the covariance
matrix :math:`\Sigma_k` is, by definition, equal to :math:`\frac{1}{n - 1}
X_k^tX_k = \frac{1}{n - 1} V S^2 V^t` where :math:`V` comes from the SVD of the (centered)
matrix: :math:`X_k = U S V^t`. It turns out that we can compute the
log-posterior above without having to explicitly compute :math:`\Sigma`:
computing :math:`S` and :math:`V` via the SVD of :math:`X` is enough. For
LDA, two SVDs are computed: the SVD of the centered input matrix :math:`X`
and the SVD of the class-wise mean vectors.

The 'lsqr' solver is an efficient algorithm that only works for
classification. It needs to explicitly compute the covariance matrix
:math:`\Sigma`, and supports shrinkage and custom covariance estimators.
This solver computes the coefficients
:math:`\omega_k = \Sigma^{-1}\mu_k` by solving for :math:`\Sigma \omega =
\mu_k`, thus avoiding the explicit computation of the inverse
:math:`\Sigma^{-1}`.

The 'eigen' solver is based on the optimization of the between class scatter to
within class scatter ratio. It can be used for both classification and
transform, and it supports shrinkage. However, the 'eigen' solver needs to
compute the covariance matrix, so it might not be suitable for situations with
a high number of features.

.. rubric:: References

.. [1] "The Elements of Statistical Learning", Hastie T., Tibshirani R.,
    Friedman J., Section 4.3, p.106-119, 2008.

.. [2] Ledoit O, Wolf M. Honey, I Shrunk the Sample Covariance Matrix.
    The Journal of Portfolio Management 30(4), 110-119, 2004.

.. [3] R. O. Duda, P. E. Hart, D. G. Stork. Pattern Classification
    (Second Edition), section 2.6.2.
