
.. _linear_model:

=========================
Generalized Linear Models
=========================

.. currentmodule:: scikits.learn.linear_model

The following are a set of methods intended for regression in which
the target value is expected to be a linear combination of the input
variables. In mathematical notion, if :math:`\hat{y}` is the predicted
value.

.. math::    \hat{y}(\beta, x) = \beta_0 + \beta_1 x_1 + ... + \beta_D x_D

Across the module, we designate the vector :math:`\beta = (\beta_1,
..., \beta_D)` as ``coef_`` and :math:`\beta_0` as ``intercept_``.

To perform classification with generalized linear models, see
:ref:`Logistic_regression`.

.. _ordinary_least_squares:

Ordinary Least Squares (OLS)
==============================

:class:`LinearRegression` fits a linear model with coefficients
:math:`\beta = (\beta_1, ..., \beta_D)` to minimize the residual sum
of squares between the observed responses in the dataset, and the
responses predicted by the linear approximation.

.. figure:: ../auto_examples/linear_model/images/plot_ols_1.png
   :target: ../auto_examples/linear_model/plot_ols.html
   :scale: 50%

:class:`LinearRegression` will take in its `fit` method arrays X, y
and will store the coefficients :math:`w` of the linear model in its
`coef\_` member.


    >>> from scikits.learn import linear_model
    >>> clf = linear_model.LinearRegression()
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

.. topic:: Examples:

   * :ref:`example_linear_model_plot_ols.py`


OLS Complexity
--------------

This method computes the least squares solution using a singular value
decomposition of X. If X is a matrix of size (n, p ) this method has a
cost of :math:`O(n p^2)`, assuming that :math:`n \geq p`.

Ridge Regression
================

:class:`Ridge` regression addresses some of the problems of
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


    >>> from scikits.learn import linear_model
    >>> clf = linear_model.Ridge (alpha = .5)
    >>> clf.fit ([[0, 0], [0, 0], [1, 1]], [0, .1, 1])
    Ridge(alpha=0.5, fit_intercept=True)
    >>> clf.coef_
    array([ 0.34545455,  0.34545455])
    >>> clf.intercept_ #doctest: +ELLIPSIS
    0.13636...

Ridge Complexity
--------------------

This method has the same order of complexity than an
:ref:`ordinary_least_squares`.

Generalized Cross-Validation
----------------------------

:class:`RidgeCV` implements ridge regression with built-in cross-validation of the alpha parameter.
The object works in the same way as GridSearchCV except that it defaults to Generalized Cross-Validation (GCV), an efficient form of leave-one-out cross-validation.

    >>> from scikits.learn import linear_model
    >>> clf = linear_model.RidgeCV(alphas=[0.1, 1.0, 10.0])
    >>> clf.fit ([[0, 0], [0, 0], [1, 1]], [0, .1, 1])
    RidgeCV(alphas=[0.10000000000000001, 1.0, 10.0], loss_func=None, cv=None,
        score_func=None, fit_intercept=True)
    >>> clf.best_alpha
    0.10000000000000001

.. topic:: References

    * "Notes on Regularized Least Squares", Rifkin & Lippert (`technical report
      <http://cbcl.mit.edu/projects/cbcl/publications/ps/MIT-CSAIL-TR-2007-025.pdf>`_,
      `course slides
      <http://www.mit.edu/~9.520/spring07/Classes/rlsslides.pdf>`_).

Lasso
=====

The :class:`Lasso` is a linear model trained with L1 prior as
regularizer. The objective function to minimize is:

.. math::  0.5 * ||y - X w||_2 ^ 2 + \alpha * ||w||_1

The lasso estimate thus solves the minimization of the
least-squares penalty with :math:`\alpha * ||w||_1` added, where
:math:`\alpha` is a constant and :math:`||w||_1` is the L1-norm of the
parameter vector.


This formulation is useful in some contexts due to its tendency to
prefer solutions with fewer parameter values, effectively reducing the
number of variables upon which the given solution is dependent. For
this reason, the Lasso and its variants are fundamental to the field
of compressed sensing.

This implementation uses coordinate descent as the algorithm to fit
the coefficients. See :ref:`least_angle_regression` for another implementation.

    >>> clf = linear_model.Lasso(alpha = 0.1)
    >>> clf.fit ([[0, 0], [1, 1]], [0, 1])
    Lasso(alpha=0.1, fit_intercept=True)
    >>> clf.predict ([[1, 1]])
    array([ 0.8])

The function :func:`lasso_path` computes the coefficients along the full path
of possible values.

.. topic:: Examples:

  * :ref:`example_linear_model_lasso_and_elasticnet.py`,
  * :ref:`example_linear_model_lasso_path_with_crossvalidation.py`


Elastic Net
===========
:class:`ElasticNet` is a linear model trained with L1 and L2 prior as
regularizer.


The objective function to minimize is in this case

.. math::        0.5 * ||y - X w||_2 ^ 2 + \alpha * \rho * ||w||_1 + \alpha * (1-\rho) * 0.5 * ||w||_2 ^ 2


.. topic:: Examples:

  * :ref:`example_linear_model_lasso_and_elasticnet.py`
  * :ref:`example_linear_model_plot_lasso_coordinate_descent_path.py`


.. _least_angle_regression:

Least Angle Regression
======================

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

  - It is effective in contexts where p >> n (i.e., when the number of
    dimensions is significantly greater than the number of points)

The disadvantages of the LARS method include:

  - Because LARS is based upon an iterative refitting of the
    residuals, it would appear to be especially sensitive to the
    effects of noise. This problem is discussed in detail by Weisberg
    in the discussion section of the Efron et al. (2004) Annals of
    Statistics article.

The LARS model can be used using estimator :class:`LARS`, or its
low-level implementation :func:`lars_path`.


LARS Lasso
==========

:class:`LassoLARS` is a lasso model implemented using the LARS
algorithm, and unlike the implementation based on coordinate_descent,
this yields the exact solution, which is piecewise linear as a
function of the norm of its coefficients.

.. figure:: ../auto_examples/linear_model/images/plot_lasso_lars_1.png
   :target: ../auto_examples/linear_model/plot_lasso_lars.html
   :align: center
   :scale: 50%

::

   >>> from scikits.learn import linear_model
   >>> clf = linear_model.LassoLARS(alpha=.1)
   >>> clf.fit ([[0, 0], [1, 1]], [0, 1])
   LassoLARS(alpha=0.1, verbose=False, fit_intercept=True)
   >>> clf.coef_
   array([ 0.30710678,  0.        ])


.. topic:: Examples:

 * :ref:`example_linear_model_plot_lasso_lars.py`


The LARS algorithm provides the full path of the coefficients along
the regularization parameter almost for free, thus a common operation
consist of retrieving the path with function :func:`lars_path`

 

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

.. topic:: References:

 * Original Algorithm is detailed in the paper `Least Angle Regression
   <http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf>`_
   by Hastie et al.



Bayesian Regression
===================

Bayesian regression techniques can be used to include regularization
parameters in the estimation procedure. This can be done by
introducing some prior knowledge over the parameters.  For example,
penalization by weighted :math:`\ell_{2}` norm is equivalent to
setting Gaussian priors on the weights.

The advantages of *Bayesian Regression* are:

    - It adapts to the data at hand.

    - It can be used to include regularization parameters in the
      estimation procedure.

The disadvantages of *Bayesian Regression* include:

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
There is also a Gamma prior for :math:`\lambda` and :math:`\alpha`:

.. math:: g(\alpha|\alpha_1,\alpha_2) = \frac{\alpha_2^{\alpha_1}}
    {\Gamma(\alpha_1)} \alpha^{\alpha_1-1} e^{-\alpha_2 {\alpha}}


.. math:: g(\lambda|\lambda_1,\lambda_2) = \frac{\lambda_2^{\lambda_1}}
    {\Gamma(\lambda_1)} \lambda^{\lambda_1-1} e^{-\lambda_2 {\lambda}}
    
By default :math:`\alpha_1 = \alpha_2 =  \lambda_1 = \lambda_2 = 1.e^{-6}`, *i.e.*
 very slightly informative priors.



.. figure:: ../auto_examples/linear_model/images/plot_bayesian_ridge_1.png
   :target: ../auto_examples/linear_model/plot_bayesian_ridge.html
   :align: center


*Bayesian Ridge Regression* is used for regression:

    >>> from scikits.learn import linear_model
    >>> X = [[0., 0.], [1., 1.], [2., 2.], [3., 3.]]
    >>> Y = [0., 1., 2., 3.]
    >>> clf = linear_model.BayesianRidge()
    >>> clf.fit (X, Y)
    BayesianRidge(n_iter=300, verbose=False, lambda_1=1e-06, lambda_2=1e-06,
           fit_intercept=True, eps=0.001, alpha_2=1e-06, alpha_1=1e-06,
           compute_score=False)

After being fitted, the model can then be used to predict new values::

    >>> clf.predict ([[1, 0.]])
    array([ 0.50000013])


The weights :math:`\beta` of the model can be access:

    >>> clf.coef_
    array([ 0.49999993,  0.49999993])

Due to the Bayesian framework, the weights found are slightly different to the
ones found by :ref:`ordinary_least_squares`. However, *Bayesian Ridge
Regression* is more robust to ill-posed problem.

.. topic:: Examples:

 * :ref:`example_linear_model_plot_bayesian_ridge.py`

.. topic:: References

  * More details can be found in the article `Bayesian Interpolation <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.27.9072&rep=rep1&type=pdf>`_ 
    by MacKay, David J. C.



Automatic Relevance Determination - ARD
=======================================

:class:`ARDRegression` adds a more sophisticated prior :math:`\beta`,
where we assume that each weight :math:`\beta_{i}` is drawn in a
Gaussian distribution, centered on zero and with a precision
:math:`\lambda_{i}`:

.. math:: p(\beta|\lambda) = \mathcal{N}(\beta|0,A^{-1})

with :math:`diag \; (A) = \lambda = \{\lambda_{1},...,\lambda_{p}\}`.
There is also a Gamma prior for :math:`\lambda` and :math:`\alpha`:

.. math:: g(\alpha|\alpha_1,\alpha_2) = \frac{\alpha_2^{\alpha_1}}
    {\Gamma(\alpha_1)} \alpha^{\alpha_1-1} e^{-\alpha_2 {\alpha}}


.. math:: g(\lambda|\lambda_1,\lambda_2) = \frac{\lambda_2^{\lambda_1}}
    {\Gamma(\lambda_1)} \lambda^{\lambda_1-1} e^{-\lambda_2 {\lambda}}
    
By default :math:`\alpha_1 = \alpha_2 =  \lambda_1 = \lambda_2 = 1.e-6`, *i.e.*
 very slightly informative priors.


.. figure:: ../auto_examples/linear_model/images/plot_ard_1.png
   :target: ../auto_examples/linear_model/plot_ard.html
   :align: center


.. topic:: Examples:

  * :ref:`example_linear_model_plot_ard.py`

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


.. topic:: References

 * Original Algorithm is detailed in the  book *Bayesian learning for neural
   networks* by Radford M. Neal

.. _Logistic_regression:

Logisitic regression
======================

If the task at hand is to do choose which class a sample belongs to given
a finite (hopefuly small) set of choices, the learning problem is a
classification, rather than regression. Linear models can be used for
such a decision, but it is best to use what is called a 
`logistic regression <http://en.wikipedia.org/wiki/Logistic_regression>`__, 
that doesn't try to minimize the sum of square residuals, as in regression, 
but rather a "hit or miss" cost.

The :class:`LogisticRegression` class can be used to do L1 or L2 penalized
logistic regression, in order to have sparse predicting weights.

.. topic:: Examples:

  * :ref:`example_logistic_l1_l2_coef.py`

  * :ref:`example_linear_model_plot_logistic_path.py`

Stochastic Gradient Descent - SGD
=================================

Stochastic gradient descent is a simple yet very efficient approach 
to fit linear models. It is particulary useful when the number of samples 
(and the number of features) is very large.


The classes :class:`SGDClassifier` and :class:`SGDRegressor` provide 
functionality to fit linear models for classification and regression 
using different (convex) loss functions and different penalties. 

.. topic:: References

 * :ref:`sgd`
