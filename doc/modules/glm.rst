=================
Linear Regression
=================


In this model, the target value is expected to be a linear combination
of the input variables.

.. math::    y(x, w) = w_0 + w_1 x_1 + ... + w_D x_D


Ordinary Least Squares
======================

Parameter w is estimated by least squares.

.. autoclass:: scikits.learn.glm.regression.LinearRegression
   :members:


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

.. autoclass:: scikits.learn.glm.ridge.Ridge
   :members:


Lasso
=====
The Lasso is a linear model trained with L1 prior as regularizer

.. autoclass:: scikits.learn.glm.coordinate_descent.Lasso
   :members:

Elastic Net
===========
Elastic Net is a linear model trained with L1 and L2 prior as
regularizer.

.. autoclass:: scikits.learn.glm.coordinate_descent.ElasticNet
   :members:


Examples
========

:ref:`example_plot_lasso_coordinate_descent_path.py`

