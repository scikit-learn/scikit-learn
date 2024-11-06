.. _cross_decomposition:

===================
Cross decomposition
===================

.. currentmodule:: sklearn.cross_decomposition

The cross decomposition module contains **supervised** estimators for
dimensionality reduction and regression, belonging to the "Partial Least
Squares" family.

.. figure:: ../auto_examples/cross_decomposition/images/sphx_glr_plot_compare_cross_decomposition_001.png
   :target: ../auto_examples/cross_decomposition/plot_compare_cross_decomposition.html
   :scale: 75%
   :align: center


Cross decomposition algorithms find the fundamental relations between two
matrices (X and Y). They are latent variable approaches to modeling the
covariance structures in these two spaces. They will try to find the
multidimensional direction in the X space that explains the maximum
multidimensional variance direction in the Y space. In other words, PLS
projects both `X` and `Y` into a lower-dimensional subspace such that the
covariance between `transformed(X)` and `transformed(Y)` is maximal.

PLS draws similarities with `Principal Component Regression
<https://en.wikipedia.org/wiki/Principal_component_regression>`_ (PCR), where
the samples are first projected into a lower-dimensional subspace, and the
targets `y` are predicted using `transformed(X)`. One issue with PCR is that
the dimensionality reduction is unsupervised, and may lose some important
variables: PCR would keep the features with the most variance, but it's
possible that features with a small variances are relevant from predicting
the target. In a way, PLS allows for the same kind of dimensionality
reduction, but by taking into account the targets `y`. An illustration of
this fact is given in the following example:
* :ref:`sphx_glr_auto_examples_cross_decomposition_plot_pcr_vs_pls.py`.

Apart from CCA, the PLS estimators are particularly suited when the matrix of
predictors has more variables than observations, and when there is
multicollinearity among the features. By contrast, standard linear regression
would fail in these cases unless it is regularized.

Classes included in this module are :class:`PLSRegression`,
:class:`PLSCanonical`, :class:`CCA` and :class:`PLSSVD`

PLSCanonical
------------

We here describe the algorithm used in :class:`PLSCanonical`. The other
estimators use variants of this algorithm, and are detailed below.
We recommend section [1]_ for more details and comparisons between these
algorithms. In [1]_, :class:`PLSCanonical` corresponds to "PLSW2A".

Given two centered matrices :math:`X \in \mathbb{R}^{n \times d}` and
:math:`Y \in \mathbb{R}^{n \times t}`, and a number of components :math:`K`,
:class:`PLSCanonical` proceeds as follows:

Set :math:`X_1` to :math:`X` and :math:`Y_1` to :math:`Y`. Then, for each
:math:`k \in [1, K]`:

- a) compute :math:`u_k \in \mathbb{R}^d` and :math:`v_k \in \mathbb{R}^t`,
  the first left and right singular vectors of the cross-covariance matrix
  :math:`C = X_k^T Y_k`.
  :math:`u_k` and :math:`v_k` are called the *weights*.
  By definition, :math:`u_k` and :math:`v_k` are
  chosen so that they maximize the covariance between the projected
  :math:`X_k` and the projected target, that is :math:`\text{Cov}(X_k u_k,
  Y_k v_k)`.
- b) Project :math:`X_k` and :math:`Y_k` on the singular vectors to obtain
  *scores*: :math:`\xi_k = X_k u_k` and :math:`\omega_k = Y_k v_k`
- c) Regress :math:`X_k` on :math:`\xi_k`, i.e. find a vector :math:`\gamma_k
  \in \mathbb{R}^d` such that the rank-1 matrix :math:`\xi_k \gamma_k^T`
  is as close as possible to :math:`X_k`. Do the same on :math:`Y_k` with
  :math:`\omega_k` to obtain :math:`\delta_k`. The vectors
  :math:`\gamma_k` and :math:`\delta_k` are called the *loadings*.
- d) *deflate* :math:`X_k` and :math:`Y_k`, i.e. subtract the rank-1
  approximations: :math:`X_{k+1} = X_k - \xi_k \gamma_k^T`, and
  :math:`Y_{k + 1} = Y_k - \omega_k \delta_k^T`.

At the end, we have approximated :math:`X` as a sum of rank-1 matrices:
:math:`X = \Xi \Gamma^T` where :math:`\Xi \in \mathbb{R}^{n \times K}`
contains the scores in its columns, and :math:`\Gamma^T \in \mathbb{R}^{K
\times d}` contains the loadings in its rows. Similarly for :math:`Y`, we
have :math:`Y = \Omega \Delta^T`.

Note that the scores matrices :math:`\Xi` and :math:`\Omega` correspond to
the projections of the training data :math:`X` and :math:`Y`, respectively.

Step *a)* may be performed in two ways: either by computing the whole SVD of
:math:`C` and only retain the singular vectors with the biggest singular
values, or by directly computing the singular vectors using the power method (cf section 11.3 in [1]_),
which corresponds to the `'nipals'` option of the `algorithm` parameter.

.. dropdown:: Transforming data

  To transform :math:`X` into :math:`\bar{X}`, we need to find a projection
  matrix :math:`P` such that :math:`\bar{X} = XP`. We know that for the
  training data, :math:`\Xi = XP`, and :math:`X = \Xi \Gamma^T`. Setting
  :math:`P = U(\Gamma^T U)^{-1}` where :math:`U` is the matrix with the
  :math:`u_k` in the columns, we have :math:`XP = X U(\Gamma^T U)^{-1} = \Xi
  (\Gamma^T U) (\Gamma^T U)^{-1} = \Xi` as desired. The rotation matrix
  :math:`P` can be accessed from the `x_rotations_` attribute.

  Similarly, :math:`Y` can be transformed using the rotation matrix
  :math:`V(\Delta^T V)^{-1}`, accessed via the `y_rotations_` attribute.

.. dropdown:: Predicting the targets `Y`

  To predict the targets of some data :math:`X`, we are looking for a
  coefficient matrix :math:`\beta \in R^{d \times t}` such that :math:`Y =
  X\beta`.

  The idea is to try to predict the transformed targets :math:`\Omega` as a
  function of the transformed samples :math:`\Xi`, by computing :math:`\alpha
  \in \mathbb{R}` such that :math:`\Omega = \alpha \Xi`.

  Then, we have :math:`Y = \Omega \Delta^T = \alpha \Xi \Delta^T`, and since
  :math:`\Xi` is the transformed training data we have that :math:`Y = X \alpha
  P \Delta^T`, and as a result the coefficient matrix :math:`\beta = \alpha P
  \Delta^T`.

  :math:`\beta` can be accessed through the `coef_` attribute.

PLSSVD
------

:class:`PLSSVD` is a simplified version of :class:`PLSCanonical`
described earlier: instead of iteratively deflating the matrices :math:`X_k`
and :math:`Y_k`, :class:`PLSSVD` computes the SVD of :math:`C = X^TY`
only *once*, and stores the `n_components` singular vectors corresponding to
the biggest singular values in the matrices `U` and `V`, corresponding to the
`x_weights_` and `y_weights_` attributes. Here, the transformed data is
simply `transformed(X) = XU` and `transformed(Y) = YV`.

If `n_components == 1`, :class:`PLSSVD` and :class:`PLSCanonical` are
strictly equivalent.

PLSRegression
-------------

The :class:`PLSRegression` estimator is similar to
:class:`PLSCanonical` with `algorithm='nipals'`, with 2 significant
differences:

- at step a) in the power method to compute :math:`u_k` and :math:`v_k`,
  :math:`v_k` is never normalized.
- at step c), the targets :math:`Y_k` are approximated using the projection
  of :math:`X_k` (i.e. :math:`\xi_k`) instead of the projection of
  :math:`Y_k` (i.e. :math:`\omega_k`). In other words, the loadings
  computation is different. As a result, the deflation in step d) will also
  be affected.

These two modifications affect the output of `predict` and `transform`,
which are not the same as for :class:`PLSCanonical`. Also, while the number
of components is limited by `min(n_samples, n_features, n_targets)` in
:class:`PLSCanonical`, here the limit is the rank of :math:`X^TX`, i.e.
`min(n_samples, n_features)`.

:class:`PLSRegression` is also known as PLS1 (single targets) and PLS2
(multiple targets). Much like :class:`~sklearn.linear_model.Lasso`,
:class:`PLSRegression` is a form of regularized linear regression where the
number of components controls the strength of the regularization.

Canonical Correlation Analysis
------------------------------

Canonical Correlation Analysis was developed prior and independently to PLS.
But it turns out that :class:`CCA` is a special case of PLS, and corresponds
to PLS in "Mode B" in the literature.

:class:`CCA` differs from :class:`PLSCanonical` in the way the weights
:math:`u_k` and :math:`v_k` are computed in the power method of step a).
Details can be found in section 10 of [1]_.

Since :class:`CCA` involves the inversion of :math:`X_k^TX_k` and
:math:`Y_k^TY_k`, this estimator can be unstable if the number of features or
targets is greater than the number of samples.

.. rubric:: References

.. [1] `A survey of Partial Least Squares (PLS) methods, with emphasis on the two-block
  case <https://stat.uw.edu/sites/default/files/files/reports/2000/tr371.pdf>`_,
  JA Wegelin

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_cross_decomposition_plot_compare_cross_decomposition.py`
* :ref:`sphx_glr_auto_examples_cross_decomposition_plot_pcr_vs_pls.py`
