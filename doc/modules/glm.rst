=========================
Generalized Linear Models
=========================

.. currentmodule:: scikits.learn.glm

The following are a set of methods intended for regression in which
the target value is expected to be a linear combination of the input
variables. In mathematical notion, if :math:`\hat{y}` is the predicted
value.

.. math::    \hat{y}(\beta, x) = \beta_0 + \beta_1 x_1 + ... + \beta_D x_D

Across the module, we designate the vector :math:`\beta = (\beta_1,
..., \beta_D)` as ``coef_`` and :math:`\beta_0` as ``intercept_``.


.. _ordinary_least_squares:

Ordinary Least Squares
======================

:class:`LinearRegression` fits a linear model with coefficients
:math:`\beta = (\beta_1, ..., \beta_D)` to minimize the residual sum
of squares between the observed responses in the dataset, and the
responses predicted by the linear approximation.


.. figure:: ../auto_examples/glm/images/plot_ols.png
   :target: ../auto_examples/glm/plot_ols.html
   :scale: 50%
   :align: center

:class:`LinearRegression` will take in its `fit` method arrays X, y
and will store the coefficients :math:`w` of the linear model in its
`coef\_` member.


    >>> from scikits.learn import glm
    >>> clf = glm.LinearRegression()
    >>> clf.fit ([[0, 0], [1, 1], [2, 2]], [0, 1, 2])
    LinearRegression(fit_intercept=True)
    >>> clf.coef_
    array([ 0.5,  0.5])


However, coefficient estimates for Ordinary Least Squares rely on the
independence of the model terms. When terms are correlated and the
columns of the design matrix :math:`X` have an approximate linear
dependence, the matrix :math:`X(X^T X)^{-1}` becomes close to singular
and as a result, the least-squares estimate becomes highly sensitive
to random errors in the observed response, producing a large
variance. This situation of *multicollinearity* can arise, for
example, when data are collected without an experimental design.


Complexity
----------

This method computes the least squares solution using a singular value
decomposition of X. If X is a matrix of size (n, p ) this method has a
cost of :math:`O(n p^2)`, assuming that :math:`n \geq p`.


Examples
--------
:ref:`example_glm_plot_ols.py`


Ridge Regression
================

:class:`Ridge` regression adresses some of the problems of
:ref:`ordinary_least_squares` by imposing a penalty on the size of
coefficients. The ridge coefficients minimize a penalized residual sum
of squares,


.. math::

   \beta^{ridge} = \underset{\beta}{argmin} { \sum_{i=1}{N} (y_i -
                 \beta_0 - \sum_{j=1}{p} x_ij \beta_j)^2 + \alpha
                 \sum_{j=1}{p} \beta_{j}^2}

Here, :math:`\alpha \geq 0` is a complexity parameter that controls
the amount of shrinkage: the larger the value of :math:`\alpha`, the
greater the amount of shrinkage.


    >>> from scikits.learn import glm
    >>> clf = glm.Ridge (alpha = .5)
    >>> clf.fit ([[0, 0], [0, 0], [1, 1]], [0, .1, 1])
    Ridge(alpha=0.5, fit_intercept=True)
    >>> clf.coef_
    array([ 0.34545455,  0.34545455])
    >>> clf.intercept_
    0.13636363636363638

Complexity
----------

This method has the same order of complexity than an
:ref:`ordinary_least_squares`.

Lasso
=====

The :class:`Lasso` is a linear model trained with L1 prior as
regularizer. The objective function to minimize is:

.. math::  0.5 * ||y - X w||_2 ^ 2 + \alpha * ||w||_1

The lasso estimate solves thus solves the minization of the
least-squares penalty with :math:`\alpha * ||w||_1` added, where
:math:`\alpha` is a constant and :math:`||w||_1` is the L1-norm of the
parameter vector.


This formulation is useful in some context due to its tendency to
prefer solutions with fewer parameter values, effectively reducing the
number of variables upon which the given solution is dependent. For
this reason, the Lasso and its variants are fundamental to the field
of compressed sensing.

This implementation uses coordinate descent as the algorithm to fit
the coeffcients. See :ref:`lars_algorithm` for another implementation.

    >>> clf = glm.Lasso(alpha = 0.1)
    >>> clf.fit ([[0, 0], [1, 1]], [0, 1])
    Lasso(alpha=0.1, coef_=array([ 0.6,  0. ]), fit_intercept=True)
    >>> clf.predict ([[1, 1]])
    array([ 0.8])

The function lasso_path computes the coefficients along the full path
of possible values.

Examples
--------
:ref:`example_glm_lasso_and_elasticnet.py`,
:ref:`example_glm_lasso_path_with_crossvalidation.py`


Elastic Net
===========
:class:`ElasticNet` is a linear model trained with L1 and L2 prior as
regularizer.


The objective function to minize is in this case

.. math::        0.5 * ||y - X w||_2 ^ 2 + \alpha * \rho * ||w||_1 + \alpha * (1-\rho) * 0.5 * ||w||_2 ^ 2


Examples
--------

:ref:`example_glm_lasso_and_elasticnet.py`
:ref:`example_plot_lasso_coordinate_descent_path.py`


.. _lars_algorithm:

LARS algorithm and its variants
===============================

Least-angle regression (LARS) is a regression algorithm for
high-dimensional data, developed by Bradley Efron, Trevor Hastie, Iain
Johnstone and Robert Tibshirani.

The advantages of LARS are:

  - It is computationally just as fast as forward selection and has
    the same order of complexity as an ordinary least squares.

  - It produces a full piecewise linear solution path, which is
    useful in cross-validation or similar attempts to tune the model.

  - If two variables are almost equally correlated with the response,
    then their coefficients should increase at approximately the same
    rate. The algorithm thus behaves as intuition would expect, and
    also is more stable.

  - It is easily modified to produce solutions for other estimators,
    like the Lasso. 

  - It is effective in contexts where p >> n (IE, when the number of
    dimensions is significantly greater than the number of points)

The disadvantages of the LARS method include:

  - Because LARS is based upon an iterative refitting of the
    residuals, it would appear to be especially sensitive to the
    effects of noise. This problem is discussed in detail by Weisberg
    in the discussion section of the Efron et al. (2004) Annals of
    Statistics article.

It is implemented using the LARS algorithm in class :class:`glm.LARS`.


LARS Lasso
==========

This implementation is based on the LARS algorithm, and unlike the
implementation based on coordinate_descent, this yields the exact
solution, which is piecewise linear as a function of the norm of its
coefficients.


   >>> clf = glm.LassoLARS(alpha=.1)
   >>> clf.fit ([[0, 0], [1, 1]], [0, 1])
   LassoLARS(normalize=True, alpha=0.1, max_iter=None)
   >>> clf.coef_
   array([ 0.50710678,  0.        ])


Getting the full path
---------------------
See function scikits.learn.glm.lars_path.


Mathematical formulation
------------------------

The algorithm is similar to forward stepwise regression, but instead
of including variables at each step, the estimated parameters are
increased in a direction equiangular to each one's correlations with
the residual.

Instead of giving a vector result, the LARS solution consists of a
curve denoting the solution for each value of the L1 norm of the
parameter vector. The full coeffients path is stored in the array
``coef_path_``, which has size (n_features, max_features+1). The first
column is always zero.


Examples
--------
:ref:`example_glm_plot_lar.py`, :ref:`example_glm_plot_lasso_lars.py`

References 
----------
Original Algorithm is detailed in the `paper
<http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf>`_
by Hastie et al.



Bayesian Regression
===================

Bayesian regression techniques can be used to include regularization parameters
in the estimation procedure. This can be done by introducing some prior
knowledge over the parameters. 
For example, penalization by weighted :math:`\ell_{2}` norm is equivalent to
setting Gaussian priors on the weights. 

The advantages of *Bayesian Regression* are:

    - It adapts to the data at hand.

    - It can be used to include regularization parameters in the
      estimation procedure.

The dissadvantages of *Bayesian Regression* include:

    - Inference of the model can be time consuming.


Bayesian Ridge Regression
-------------------------

:class:`BayesianRidge` tries to avoid the overfit issue of
:ref:`ordinary_least_squares`, by adding the following prior on
:math:`\beta`:

.. math:: p(\beta|\lambda) =  
    \mathcal{N}(\beta|0,\lambda^{-1}\bold{I_{p}})

The resulting model is called *Bayesian Ridge Regression*, it is
similar to the classical :class:`Ridge`.  :math:`\lambda` is an
*hyper-parameter* and the prior over :math:`\beta` performs a
shrinkage or regularization, by constraining the values of the weights
to be small. Indeed, with a large value of :math:`\lambda`, the
Gaussian is narrowed around 0 which does not allow large values of
:math:`\beta`, and with low value of :math:`\lambda`, the Gaussian is
very flattened which allows values of :math:`\beta`.  Here, we use a
*non-informative* prior for :math:`\lambda`.


The parameters are estimated by maximizing the *marginal log likelihood*.


.. figure:: ../auto_examples/glm/images/plot_bayesian_ridge.png
   :target: ../auto_examples/glm/plot_bayesian_ridge.html
   :align: center


*Bayesian Ridge Regression* is used for regression:

    >>> from scikits.learn import glm
    >>> X = [[0., 0.], [1., 2.], [2., -1.], [3., -1.]]
    >>> Y = [0., 1., 2., 3.]
    >>> clf = glm.BayesianRidge()
    >>> clf.fit (X, Y)
    BayesianRidge(n_iter=300, th_w=1e-12, compute_ll=False, fit_intercept=True)

After being fitted, the model can then be used to predict new values::

    >>> clf.predict ([[-1.5, 0.]])
    array([ 1.5])


The weights :math:`\beta` of the model can be access:

    >>> clf.coef_
    array([ 0.93688528, -0.03034525])

Due to the Bayesian framework, the weights found are slightly different to the
ones found by :ref:`ordinary_least_squares`. However, *Bayesian Ridge
Regression* is more robust to ill-posed problem.


Examples
--------
:ref:`example_glm_plot_bayesian_ridge.py`

References
----------
More details can be found in the article
`paper
<http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.27.9072&rep=rep1&type=
pdf>`_ by MacKay, David J. C.





Automatic Relevance Determination - ARD
=======================================

:class:`ARDRegression` adds a more sophisticated prior :math:`\beta`,
where we assume that each weight :math:`\beta_{i}` is drawn in a
Gaussian distribution, centered on zero and with a precision
:math:`\lambda_{i}`:

.. math:: p(\beta|\lambda) = \mathcal{N}(\beta|0,A^{-1})

with :math:`diag \; (A) = lambda = \{\lambda_{1},...,\lambda_{p}\}`.
We use a *non-informative* prior for :math:`\lambda`.


.. figure:: ../auto_examples/glm/images/plot_ard.png
   :target: ../auto_examples/glm/plot_ard.html
   :align: center


Examples
--------
:ref:`example_glm_plot_ard.py`

Mathematical formulation
------------------------

A prior is introduced as a distribution :math:`p(\theta)` over the parameters.
This distribution is set before processing the data. The parameters of a prior
distribution are called *hyper-parameters*. This description is based on the
Bayes theorem :

.. math:: p(\theta|\{X,y\})
   = \frac{p(\{X,y\}|\theta)p(\theta)}{p(\{X,y\})}

With :
    - :math:`p({X, y}|\theta)` the likelihood : it expresses how probable it is
      to observe :math:`{X,y}` given :math:`\theta`.

    - :math:`p({X, y})` the marginal probability of the data : it can be
      considered as a normalizing constant, and is computed by integrating
      :math:`p({X, y}|\theta)` with respect to :math:`\theta`.

    - :math:`p(\theta)` the prior over the parameters : it expresses the
      knowledge that we can have about :math:`\theta` before processing the
      data.

    - :math:`p(\theta|{X, y})` the conditional probability (or posterior
      probability) : it expresses the uncertainty in :math:`\theta`  after
      observing the data.


All the following regressions are based on the following Gaussian
assumption:

.. math::  p(y|X,w,\alpha) = \mathcal{N}(y|X w,\alpha)

where :math:`\alpha` is the precision of the noise.



References
----------
Original Algorithm is detailed in the  book *Bayesian learning for neural
networks* by Radford M. Neal

