=========================
Generalized Linear Models
=========================

.. currentmodule:: scikits.learn.glm

The following are a set of methods intended for regression in which
the target value is expected to be a linear combination of the input
variables. In mathematical notion, if :math:`\hat{y}` is the predicted
value.

.. math::    \hat{y}(x, w) = w_0 + w_1 x_1 + ... + w_D x_D

Across the module, we designate the vector :math:`w = (w1, ..., w_D)` as
``coef_`` and :math:`w_0` as ``intercept_``.


Ordinary Least Squares
======================

This method minimizes the sum of squared distances between the
observed responses in the dataset, and the responses predicted by the
linear approximation.

Class :class:`LinearRegression` fits a model using the algorithm of
Ordinary Least Squares. It will take in its `fit` method arrays X, y
and will store the coefficients of the linear model in its `coef\_`
method.

Examples
--------
:ref:`example_glm_plot_ols.py`


Regularized Least Squares
=========================

In the following methods, We add a regularization term to an error
function in order to control over-fitting.


Ridge Regression
================

Coefficient estimates for multiple linear regression models rely on
the independence of the model terms. When terms are correlated and the
columns of the design matrix :math:`X` have an approximate linear
dependence, the matrix :math:`X(X^T X)^{-1}` becomes close to
singular. As a result, the least-squares estimate:

.. math::    \hat{\beta} = (X^T X)^{-1} X^T y

becomes highly sensitive to random errors in the observed response
:math:`y`, producing a large variance. This situation of
*multicollinearity* can arise, for example, when data are collected
without an experimental design.

Ridge regression adresses the problem by estimating regression
coefficients using:

.. math::    \hat{\beta} = (X^T X + \alpha I)^{-1} X^T y

Reference
---------
.. autoclass:: scikits.learn.glm.Ridge
   :members:


Lasso
=====

The Lasso is a linear model trained with L1 prior as regularizer. The
objective function to minimize is:

.. math::  0.5 * ||y - X w||_2 ^ 2 + \alpha * ||w||_1

The lasso estimate solves thus solves the minization of the
least-squares penalty with :math:`\alpha * ||w||_1` added, where
:math:`\alpha` is a constant and :math:`||w||_1` is the L1-norm of the
parameter vector.


This formulation is useful in some context due to its tendency to
prefer solutions with fewer parameter values, effectively reducing
the number of variables upon which the given solution is
dependent. For this reason, the LASSO and its variants are
fundamental to the field of compressed sensing.

This implementation uses coordinate descent as the algorithm to fit
the coeffcients. 

Reference
---------

.. autoclass:: scikits.learn.glm.Lasso
   :members:
   :inherited-members:


The function lasso_path computes the coefficients along the full path
of possible values.

.. autofunction:: scikits.learn.glm.lasso_path


Elastic Net
===========
Elastic Net is a linear model trained with L1 and L2 prior as
regularizer.


The objective function to minize is in this case

.. math::        0.5 * ||y - X w||_2 ^ 2 + \alpha * \rho * ||w||_1 + \alpha * (1-\rho) * 0.5 * ||w||_2 ^ 2

Reference
---------

.. autoclass:: scikits.learn.glm.ElasticNet
   :members:


Examples
--------

:ref:`example_plot_lasso_coordinate_descent_path.py`



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


Getting the full path
---------------------
See function :function:`scikits.learn.glm.lars_path`.


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
Original Algorithm is detailed in the 
`paper <http://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf>`_ by Hastie et al.

