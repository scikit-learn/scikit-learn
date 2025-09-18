.. _polynomial_chaos:

===========================
Polynomial Chaos expansions
===========================

.. currentmodule:: sklearn.polynomial_chaos

The **Polynomial Chaos Expansion (PCE)** is a supervised learning method
commonly used in uncertainty quantification (UQ) to represent and propagate
uncertainties in mathematical models.

The advantages of Polynomial Chaos expansions are:

    - Efficient quantification of uncertainty in mathematical models, by
      representing input features as random variables and expressing the
      solution as a polynomial expansion in terms of the input features.

    - Fast convergence using only a limited number of model evaluation in case
      the model output depends smoothly on the input features.

    - Insights into the sensitivity of the model output with respect to the
      different input features. Sensitivity indices can easily be calculated in
      a post-processing step to identify which input features contribute most to
      the overall uncertainty.

The disadvantages of Polynomial Chaos expansions include:

    - Not well suited for highly nonlinear or discontinuous models, as
      these models cannot be well-approximated by polynomials.

    - Lose efficiency in high-dimensional spaces, when the number
      of input features exceeds a few dozens.


.. _pcr:

Polynomial Chaos Regression (PCR)
=================================

.. currentmodule:: sklearn.polynomial_chaos

The :class:`PolynomialChaosRegressor` implements the Polynomial Chaos method for
regression purposes. Given the assumed probability distributions for each
feature, PCR will transform the input features into
:class:`~sklearn.preprocessing.OrthogonalPolynomialFeatures`, where the type of
orthogonal polynomial for each feature depends on the given input distribution.
The coefficients in the PC expansion are then computed by solving a linear
regression problem.

.. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_simple_1d_001.png
   :target: ../auto_examples/polynomial_chaos/plot_simple_1d.html
   :align: center

The figure above illustrates how Polynomial Chaos methods work when there is
only one feature. In case there is more than one feature, we need to specify
which (combinations of) polynomial terms we want in the basis of the PCE. We
allow two different ways to specify the basis terms, following the
specifications for multiindex sets in
:ref:`generating_orthogonal_polynomial_features`.

* In a first approach, we construct a multiindex set based on a given `degree`,
  `truncation` rule and (optional) `weights`. There are 4 different types of
  predefind multiindex set shapes: `full_tensor`, `total_degree`,
  `hyperbolic_cross` and `Zaremba_cross`. These multiindex set shapes are
  illustrated below.

  .. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_index_sets_002.png
     :target: ../auto_examples/polynomial_chaos/plot_index_sets.html
     :align: center

  The weights allow one to select certain preferential features where the
  polynomial order should be higher.

* In a second approach, we construct a custom multiindex set based on a given
  set of multiindices. This approach allows for even more flexibility in
  specifying the basis terms to retain.

The type of orthogonal polynomial in the expansion depends on the given distributions
for the input features:

   * uniform distributions (`scipy.stats.uniform`) are associated with `Legendre
     polynomials <https://en.wikipedia.org/wiki/Legendre_polynomials>`_
   * normal distributions (`scipy.stats.norm`) are associated with `Hermite
     polynomials <https://en.wikipedia.org/wiki/Hermite_polynomials>`_
   * exponential distributions (`scipy.stats.expon`) are associated with
     `Laguerre polynomials
     <https://en.wikipedia.org/wiki/Laguerre_polynomials>`_
   * Beta distributions (`scipy.stats.beta`) are associated with `Jacobi
     polynomials <https://en.wikipedia.org/wiki/Jacobi_polynomials>`_

After the transformation into
:class:`~sklearn.preprocessing.OrthogonalPolynomialFeatures`, the coefficients of
the Polynomial Chaos expansion are found by solving a linear regression problem.

In real-world problems, the number of uncertain inputs (i.e., the number of
input features, or the dimension of the polynomial basis) and the required
degree of the orthogonal polynomial basis can be very large. When only a limited
number of model evaluations are available, the resulting linear system we need
to solve to obtain the Polynomial Chaos coefficients can be underdetermined
(that is, we have more unknowns than measurements). In such cases,
`regularization <https://en.wikipedia.org/wiki/Regularization_(mathematics)>`_ of
the problem may be required to find a good Polynomial Chaos surrogate. We can
take advantage of the various linear models, such as
:class:`~sklearn.linear_model.LassoCV`, by specifying the `solver` argument in
the :class:`PolynomialChaosRegressor`.

For more details and a technical description of Polynomial Chaos expansions, we
refer to, e.g., [Sudret2008]_.

Usage tips
==========

* Polynomial Chaos expansions work best for a moderate amount of features (say,
  a few dozen). However, always keep in mind that the complexity of the method
  scales with the number of polynomial basis terms. For the default
  `total_degree` multiindex set shape, the number of polynomial terms is

  .. math::
    P = \frac{(d + k)!}{d!k!}

  where :math:`d` is the dimension (i.e., the number of input features), and
  :math:`k` is the degree of the polynomial. The table below shows the value of
  :math:`P` for various combinations of :math:`d` (columns) and :math:`k`
  (rows). Note that the values in the table are symmetric (i.e., the role of
  :math:`d` and :math:`k` can be interchanged).

  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |   :math:`k` \\ :math:`d` |   1 |   2 |   3 |    4 |    5 |    6 |     7 |     8 |     9 |     10 |
  +==========================+=====+=====+=====+======+======+======+=======+=======+=======+========+
  |                        1 |   2 |   3 |   4 |    5 |    6 |    7 |     8 |     9 |    10 |     11 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                        2 |   3 |   6 |  10 |   15 |   21 |   28 |    36 |    45 |    55 |     66 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                        3 |   4 |  10 |  20 |   35 |   56 |   84 |   120 |   165 |   220 |    286 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                        4 |   5 |  15 |  35 |   70 |  126 |  210 |   330 |   495 |   715 |   1001 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                        5 |   6 |  21 |  56 |  126 |  252 |  462 |   792 |  1287 |  2002 |   3003 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                        6 |   7 |  28 |  84 |  210 |  462 |  924 |  1716 |  3003 |  5005 |   8008 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                        7 |   8 |  36 | 120 |  330 |  792 | 1716 |  3432 |  6435 | 11440 |  19448 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                        8 |   9 |  45 | 165 |  495 | 1287 | 3003 |  6435 | 12870 | 24310 |  43758 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                        9 |  10 |  55 | 220 |  715 | 2002 | 5005 | 11440 | 24310 | 48620 |  92378 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+
  |                       10 |  11 |  66 | 286 | 1001 | 3003 | 8008 | 19448 | 43758 | 92378 | 184756 |
  +--------------------------+-----+-----+-----+------+------+------+-------+-------+-------+--------+

* When using the default solver
  (:class:`~sklearn.linear_model.LinearRegression`), and for large :math:`d` or
  :math:`k`, the resulting linear system for the Polynomial Chaos coefficients
  can become undetermined, meaning that we have fewer model evaluations than
  polynomial basis terms. In such cases, make sure to use regularization by
  choosing a different solver, such as :class:`~sklearn.linear_model.LassoCV`.
  See also :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_simple_1d.py`.

* When the output data is noisy, as is the case in almost all real-world
  problems, use regularization to avoid overfitting. See also
  :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_noisy_data.py`.

* When using a sparse solver such as :class:`~sklearn.linear_model.LassoCV` or
  :class:`~sklearn.linear_model.ElasticNetCV`, many of the basis terms might be
  small, meaning that they can be removed from the expansion without impacting
  the prediction accuracy. This can be done by pruning the basis terms using,
  for example, :class:`~sklearn.feature_selection.SelectFromModel`. See also
  :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_feature_selection.py`.

Polynomial Chaos examples
=========================

.. topic:: Overview

   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_simple_1d.py`
   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_index_sets.py`
   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_ishigami.py`
   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_sobol_g.py`
   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_noisy_data.py`
   * :ref:`sphx_glr_auto_examples_polynomial_chaos_plot_pcr_feature_selection.py`

Simple one-dimensional model
----------------------------

This example illustrates how Polynomial Chaos expansions work by constructing a
surrogate for a one-dimensional test problem. The first figure below shows the
approximation of the model output for various polynomial degrees. As expected,
the approximation improves as the maximum degree of the orthogonal polynomial
basis increases.

However, when the polynomial order increases too much, we risk overfitting the
model. This is illustrated in the second figure, where we plot the training and
test error as a function of the polynomial degree. As the polynomial degree
increases, the training error continues to decrease. However, the test error
starts to increase when the polynomial order becomes too high. This example also
illustrates how to use regularization to avoid overfitting.

.. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_simple_1d_001.png
   :target: ../auto_examples/polynomial_chaos/plot_simple_1d.html
   :align: center

.. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_simple_1d_002.png
   :target: ../auto_examples/polynomial_chaos/plot_simple_1d.html
   :align: center

Choosing appropriate multiindex set shapes
------------------------------------------

This example discusses how to select the basis terms in a Polynomial Chaos
expansion. We consider a simple polynomial model described in [Saltelli2000]_.
The figure below show 4 different multiindex set shapes of degree 6. Since it is
hard to know a priori which multiindex set shape (and which maximum polynomial
degree) will give the best result for a given set of training points, this
example also describes how we can use the
:class:`~sklearn.model_selection.GridSearchCV` method to help us select the best
combination by minimizing the root mean squared error.

.. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_index_sets_002.png
   :target: ../auto_examples/polynomial_chaos/plot_index_sets.html
   :align: center

Global sensitivity analysis of the Ishigami function
----------------------------------------------------

This example illustrates how to use Polynomial Chaos regression to perform a
Global Sensitivity Analysis of the Ishigami function. The Ishigami function is a
well-known test problem in uncertainty quantification and sensitivity analysis
that exhibits strong nonlinearity and nonmonotonicity [Saltelli2000]_. The
figure below shows the convergence of the main-effect sensitivity indices as a
function of the number of training points. Observe that the values of the
sensitivity indices converge rapidly to the exact values (indicated by the
dashed lines). We also compare these results to a classic *pick-and-freeze*
approach to perform Global Sensitivity Analysis, and demonstrate how the
sensitivity indices computed from the polynomial expansion are far superior in
terms of accuracy compared to the classic, sampling-based approach.

.. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_pcr_ishigami_004.png
   :target: ../auto_examples/polynomial_chaos/plot_pcr_ishigami.html
   :align: center

Adaptive basis construction for the Sobol :math:`G` function
------------------------------------------------------------

This example illustrates how to use Polynomial Chaos regression with adaptive
basis growth to perform a Global Sensitivity Analysis of the Sobol function. The
Sobol function (also known as the :math:`G` function, because it was denoted as
such in the original paper by Sobol), is another well-known test problem in
uncertainty quantification and sensitivity analysis [Saltelli2000]_.  We
consider the model problem in :math:`d = 4` dimensions. The figure below shows
the values of the main sensitivity indices extracted from a second-order
accurate Polynomial Chaos surrogate model. Despite the relatively low dimension
of the problem, this example demonstrates the challenges faced when constructing
higher-order accurate Polynomial Chaos surrogates. We then discuss the adaptive
basis construction approach from [Gerstner2003]_., that somewhat alleviates this
issue.

.. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_pcr_sobol_g_001.png
   :target: ../auto_examples/polynomial_chaos/plot_pcr_sobol_g.html
   :align: center

Polynomial Chaos expansions with noisy measurements
---------------------------------------------------

This example illustrates how to use Polynomial Chaos regression with noisy data,
using a model taken from [Storlie2009]_. We demonstrate the challenges faced
when constructing accurate surrogate models when only noisy model evaluations
are available. By using regularization, we are still able to construct a
reasonably accurate Polynomial Chaos surrogate model, as illustrated by the
parity plot shown in the figure below.

.. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_pcr_noisy_data_004.png
   :target: ../auto_examples/polynomial_chaos/plot_pcr_noisy_data.html
   :align: center

Pruning basis terms in Polynomial Chaos regression
--------------------------------------------------

This example illustrates how to use Polynomial Chaos regression with basis
pruning. By default, a Polynomial Chaos expansion of degree :math:`k` includes
all multivariate polynomial terms up to that degree, which can lead to large
models with many small terms. Using a feature selection method like
:class:`~sklearn.feature_selection.SelectFromModel`, we can automatically prune
away irrelevant polynomial terms during fitting. This results in a sparse and
more interpretable expansion that can be faster to evaluate during prediction.

.. figure:: ../auto_examples/polynomial_chaos/images/sphx_glr_plot_pcr_feature_selection_001.png
   :target: ../auto_examples/polynomial_chaos/plot_pcr_feature_selection.html
   :align: center

References
----------

.. [Gerstner2003] `Gerstner, T., and Griebel, M. "Dimension-adaptive
    tensor-product quadrature." Computing 71, 2003.
    <https://link.springer.com/content/pdf/10.1007/s00607-003-0015-5.pdf>`_

.. [Saltelli2000] `Saltelli, A., et al. "Global sensitivity analysis: the
   primer." John Wiley & Sons, 2008.
   <https://onlinelibrary.wiley.com/doi/book/10.1002/9780470725184>`_

.. [Storlie2009] `Storlie, C. B., Swiler, L. P., Helton, J. C., and Sallaberry,
   C. J. "Implementation and evaluation of nonparametric regression procedures
   for sensitivity analysis of computationally demanding models." Reliability
   Engineering and System Safety, 94(11), 2009.
   <https://www.sciencedirect.com/science/article/abs/pii/S0951832009001112>`_

.. [Sudret2008] `Sudret, B. "Global sensitivity analysis using polynomial
   chaos expansions." Reliability engineering & system safety 93(7), 2008.
   <https://www.sciencedirect.com/science/article/pii/S0951832007001329>`_

.. currentmodule:: sklearn.polynomial_chaos
