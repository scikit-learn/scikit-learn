.. _linear_model:

=============
Linear Models
=============

.. currentmodule:: sklearn.linear_model

The following are a set of methods intended for regression in which
the target value is expected to be a linear combination of the features.
In mathematical notation, if :math:`\hat{y}` is the predicted
value.

.. math::    \hat{y}(w, x) = w_0 + w_1 x_1 + ... + w_p x_p

Across the module, we designate the vector :math:`w = (w_1,
..., w_p)` as ``coef_`` and :math:`w_0` as ``intercept_``.

To perform classification with generalized linear models, see
:ref:`Logistic_regression`.

.. _ordinary_least_squares:

Ordinary Least Squares
=======================

:class:`LinearRegression` fits a linear model with coefficients
:math:`w = (w_1, ..., w_p)` to minimize the residual sum
of squares between the observed targets in the dataset, and the
targets predicted by the linear approximation. Mathematically it
solves a problem of the form:

.. math:: \min_{w} || X w - y||_2^2

.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_ols_001.png
   :target: ../auto_examples/linear_model/plot_ols.html
   :align: center
   :scale: 50%

:class:`LinearRegression` will take in its ``fit`` method arrays ``X``, ``y``
and will store the coefficients :math:`w` of the linear model in its
``coef_`` member::

    >>> from sklearn import linear_model
    >>> reg = linear_model.LinearRegression()
    >>> reg.fit([[0, 0], [1, 1], [2, 2]], [0, 1, 2])
    LinearRegression()
    >>> reg.coef_
    array([0.5, 0.5])

The coefficient estimates for Ordinary Least Squares rely on the
independence of the features. When features are correlated and the
columns of the design matrix :math:`X` have an approximately linear
dependence, the design matrix becomes close to singular
and as a result, the least-squares estimate becomes highly sensitive
to random errors in the observed target, producing a large
variance. This situation of *multicollinearity* can arise, for
example, when data are collected without an experimental design.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_ols.py`

Non-Negative Least Squares
--------------------------

It is possible to constrain all the coefficients to be non-negative, which may
be useful when they represent some physical or naturally non-negative
quantities (e.g., frequency counts or prices of goods).
:class:`LinearRegression` accepts a boolean ``positive``
parameter: when set to `True` `Non-Negative Least Squares
<https://en.wikipedia.org/wiki/Non-negative_least_squares>`_ are then applied.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_nnls.py`

Ordinary Least Squares Complexity
---------------------------------

The least squares solution is computed using the singular value
decomposition of X. If X is a matrix of shape `(n_samples, n_features)`
this method has a cost of
:math:`O(n_{\text{samples}} n_{\text{features}}^2)`, assuming that
:math:`n_{\text{samples}} \geq n_{\text{features}}`.

.. _ridge_regression:

Ridge regression and classification
===================================

Regression
----------

:class:`Ridge` regression addresses some of the problems of
:ref:`ordinary_least_squares` by imposing a penalty on the size of the
coefficients. The ridge coefficients minimize a penalized residual sum
of squares:


.. math::

   \min_{w} || X w - y||_2^2 + \alpha ||w||_2^2


The complexity parameter :math:`\alpha \geq 0` controls the amount
of shrinkage: the larger the value of :math:`\alpha`, the greater the amount
of shrinkage and thus the coefficients become more robust to collinearity.

.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_ridge_path_001.png
   :target: ../auto_examples/linear_model/plot_ridge_path.html
   :align: center
   :scale: 50%


As with other linear models, :class:`Ridge` will take in its ``fit`` method
arrays ``X``, ``y`` and will store the coefficients :math:`w` of the linear model in
its ``coef_`` member::

    >>> from sklearn import linear_model
    >>> reg = linear_model.Ridge(alpha=.5)
    >>> reg.fit([[0, 0], [0, 0], [1, 1]], [0, .1, 1])
    Ridge(alpha=0.5)
    >>> reg.coef_
    array([0.34545455, 0.34545455])
    >>> reg.intercept_
    0.13636...

Note that the class :class:`Ridge` allows for the user to specify that the
solver be automatically chosen by setting `solver="auto"`. When this option
is specified, :class:`Ridge` will choose between the `"lbfgs"`, `"cholesky"`,
and `"sparse_cg"` solvers. :class:`Ridge` will begin checking the conditions
shown in the following table from top to bottom. If the condition is true,
the corresponding solver is chosen.

+-------------+----------------------------------------------------+
| **Solver**  | **Condition**                                      |
+-------------+----------------------------------------------------+
| 'lbfgs'     | The ``positive=True`` option is specified.         |
+-------------+----------------------------------------------------+
| 'cholesky'  | The input array X is not sparse.                   |
+-------------+----------------------------------------------------+
| 'sparse_cg' | None of the above conditions are fulfilled.        |
+-------------+----------------------------------------------------+


Classification
--------------

The :class:`Ridge` regressor has a classifier variant:
:class:`RidgeClassifier`. This classifier first converts binary targets to
``{-1, 1}`` and then treats the problem as a regression task, optimizing the
same objective as above. The predicted class corresponds to the sign of the
regressor's prediction. For multiclass classification, the problem is
treated as multi-output regression, and the predicted class corresponds to
the output with the highest value.

It might seem questionable to use a (penalized) Least Squares loss to fit a
classification model instead of the more traditional logistic or hinge
losses. However, in practice, all those models can lead to similar
cross-validation scores in terms of accuracy or precision/recall, while the
penalized least squares loss used by the :class:`RidgeClassifier` allows for
a very different choice of the numerical solvers with distinct computational
performance profiles.

The :class:`RidgeClassifier` can be significantly faster than e.g.
:class:`LogisticRegression` with a high number of classes because it can
compute the projection matrix :math:`(X^T X)^{-1} X^T` only once.

This classifier is sometimes referred to as a `Least Squares Support Vector
Machines
<https://en.wikipedia.org/wiki/Least-squares_support-vector_machine>`_ with
a linear kernel.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_ridge_path.py`
* :ref:`sphx_glr_auto_examples_text_plot_document_classification_20newsgroups.py`
* :ref:`sphx_glr_auto_examples_inspection_plot_linear_model_coefficient_interpretation.py`

Ridge Complexity
----------------

This method has the same order of complexity as
:ref:`ordinary_least_squares`.

.. FIXME:
.. Not completely true: OLS is solved by an SVD, while Ridge is solved by
.. the method of normal equations (Cholesky), there is a big flop difference
.. between these


Setting the regularization parameter: leave-one-out Cross-Validation
--------------------------------------------------------------------

:class:`RidgeCV` and :class:`RidgeClassifierCV` implement ridge
regression/classification with built-in cross-validation of the alpha parameter.
They work in the same way as :class:`~sklearn.model_selection.GridSearchCV` except
that it defaults to efficient Leave-One-Out :term:`cross-validation`.
When using the default :term:`cross-validation`, alpha cannot be 0 due to the
formulation used to calculate Leave-One-Out error. See [RL2007]_ for details.

Usage example::

    >>> import numpy as np
    >>> from sklearn import linear_model
    >>> reg = linear_model.RidgeCV(alphas=np.logspace(-6, 6, 13))
    >>> reg.fit([[0, 0], [0, 0], [1, 1]], [0, .1, 1])
    RidgeCV(alphas=array([1.e-06, 1.e-05, 1.e-04, 1.e-03, 1.e-02, 1.e-01, 1.e+00, 1.e+01,
          1.e+02, 1.e+03, 1.e+04, 1.e+05, 1.e+06]))
    >>> reg.alpha_
    0.01

Specifying the value of the :term:`cv` attribute will trigger the use of
cross-validation with :class:`~sklearn.model_selection.GridSearchCV`, for
example `cv=10` for 10-fold cross-validation, rather than Leave-One-Out
Cross-Validation.

.. dropdown:: References

  .. [RL2007] "Notes on Regularized Least Squares", Rifkin & Lippert (`technical report
    <http://cbcl.mit.edu/publications/ps/MIT-CSAIL-TR-2007-025.pdf>`_,
    `course slides <https://www.mit.edu/~9.520/spring07/Classes/rlsslides.pdf>`_).

.. _lasso:

Lasso
=====

The :class:`Lasso` is a linear model that estimates sparse coefficients.
It is useful in some contexts due to its tendency to prefer solutions
with fewer non-zero coefficients, effectively reducing the number of
features upon which the given solution is dependent. For this reason,
Lasso and its variants are fundamental to the field of compressed sensing.
Under certain conditions, it can recover the exact set of non-zero
coefficients (see
:ref:`sphx_glr_auto_examples_applications_plot_tomography_l1_reconstruction.py`).

Mathematically, it consists of a linear model with an added regularization term.
The objective function to minimize is:

.. math::  \min_{w} { \frac{1}{2n_{\text{samples}}} ||X w - y||_2 ^ 2 + \alpha ||w||_1}

The lasso estimate thus solves the minimization of the
least-squares penalty with :math:`\alpha ||w||_1` added, where
:math:`\alpha` is a constant and :math:`||w||_1` is the :math:`\ell_1`-norm of
the coefficient vector.

The implementation in the class :class:`Lasso` uses coordinate descent as
the algorithm to fit the coefficients. See :ref:`least_angle_regression`
for another implementation::

    >>> from sklearn import linear_model
    >>> reg = linear_model.Lasso(alpha=0.1)
    >>> reg.fit([[0, 0], [1, 1]], [0, 1])
    Lasso(alpha=0.1)
    >>> reg.predict([[1, 1]])
    array([0.8])

The function :func:`lasso_path` is useful for lower-level tasks, as it
computes the coefficients along the full path of possible values.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_and_elasticnet.py`
* :ref:`sphx_glr_auto_examples_applications_plot_tomography_l1_reconstruction.py`
* :ref:`sphx_glr_auto_examples_inspection_plot_linear_model_coefficient_interpretation.py`


.. note:: **Feature selection with Lasso**

      As the Lasso regression yields sparse models, it can
      thus be used to perform feature selection, as detailed in
      :ref:`l1_feature_selection`.

.. dropdown:: References

  The following two references explain the iterations
  used in the coordinate descent solver of scikit-learn, as well as
  the duality gap computation used for convergence control.

  * "Regularization Path For Generalized linear Models by Coordinate Descent",
    Friedman, Hastie & Tibshirani, J Stat Softw, 2010 (`Paper
    <https://www.jstatsoft.org/article/view/v033i01/v33i01.pdf>`__).
  * "An Interior-Point Method for Large-Scale L1-Regularized Least Squares,"
    S. J. Kim, K. Koh, M. Lustig, S. Boyd and D. Gorinevsky,
    in IEEE Journal of Selected Topics in Signal Processing, 2007
    (`Paper <https://web.stanford.edu/~boyd/papers/pdf/l1_ls.pdf>`__)

Setting regularization parameter
--------------------------------

The ``alpha`` parameter controls the degree of sparsity of the estimated
coefficients.

Using cross-validation
^^^^^^^^^^^^^^^^^^^^^^^

scikit-learn exposes objects that set the Lasso ``alpha`` parameter by
cross-validation: :class:`LassoCV` and :class:`LassoLarsCV`.
:class:`LassoLarsCV` is based on the :ref:`least_angle_regression` algorithm
explained below.

For high-dimensional datasets with many collinear features,
:class:`LassoCV` is most often preferable. However, :class:`LassoLarsCV` has
the advantage of exploring more relevant values of `alpha` parameter, and
if the number of samples is very small compared to the number of
features, it is often faster than :class:`LassoCV`.

.. |lasso_cv_1| image:: ../auto_examples/linear_model/images/sphx_glr_plot_lasso_model_selection_002.png
    :target: ../auto_examples/linear_model/plot_lasso_model_selection.html
    :scale: 48%

.. |lasso_cv_2| image:: ../auto_examples/linear_model/images/sphx_glr_plot_lasso_model_selection_003.png
    :target: ../auto_examples/linear_model/plot_lasso_model_selection.html
    :scale: 48%

.. centered:: |lasso_cv_1| |lasso_cv_2|

.. _lasso_lars_ic:

Information-criteria based model selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, the estimator :class:`LassoLarsIC` proposes to use the
Akaike information criterion (AIC) and the Bayes Information criterion (BIC).
It is a computationally cheaper alternative to find the optimal value of alpha
as the regularization path is computed only once instead of k+1 times
when using k-fold cross-validation.

Indeed, these criteria are computed on the in-sample training set. In short,
they penalize the over-optimistic scores of the different Lasso models by
their flexibility (cf. to "Mathematical details" section below).

However, such criteria need a proper estimation of the degrees of freedom of
the solution, are derived for large samples (asymptotic results) and assume the
correct model is candidates under investigation. They also tend to break when
the problem is badly conditioned (e.g. more features than samples).

.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_lasso_lars_ic_001.png
    :target: ../auto_examples/linear_model/plot_lasso_lars_ic.html
    :align: center
    :scale: 50%

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_model_selection.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_lars_ic.py`

.. _aic_bic:

AIC and BIC criteria
^^^^^^^^^^^^^^^^^^^^

The definition of AIC (and thus BIC) might differ in the literature. In this
section, we give more information regarding the criterion computed in
scikit-learn.

.. dropdown:: Mathematical details

  The AIC criterion is defined as:

  .. math::
      AIC = -2 \log(\hat{L}) + 2 d

  where :math:`\hat{L}` is the maximum likelihood of the model and
  :math:`d` is the number of parameters (as well referred to as degrees of
  freedom in the previous section).

  The definition of BIC replace the constant :math:`2` by :math:`\log(N)`:

  .. math::
      BIC = -2 \log(\hat{L}) + \log(N) d

  where :math:`N` is the number of samples.

  For a linear Gaussian model, the maximum log-likelihood is defined as:

  .. math::
      \log(\hat{L}) = - \frac{n}{2} \log(2 \pi) - \frac{n}{2} \ln(\sigma^2) - \frac{\sum_{i=1}^{n} (y_i - \hat{y}_i)^2}{2\sigma^2}

  where :math:`\sigma^2` is an estimate of the noise variance,
  :math:`y_i` and :math:`\hat{y}_i` are respectively the true and predicted
  targets, and :math:`n` is the number of samples.

  Plugging the maximum log-likelihood in the AIC formula yields:

  .. math::
      AIC = n \log(2 \pi \sigma^2) + \frac{\sum_{i=1}^{n} (y_i - \hat{y}_i)^2}{\sigma^2} + 2 d

  The first term of the above expression is sometimes discarded since it is a
  constant when :math:`\sigma^2` is provided. In addition,
  it is sometimes stated that the AIC is equivalent to the :math:`C_p` statistic
  [12]_. In a strict sense, however, it is equivalent only up to some constant
  and a multiplicative factor.

  At last, we mentioned above that :math:`\sigma^2` is an estimate of the
  noise variance. In :class:`LassoLarsIC` when the parameter `noise_variance` is
  not provided (default), the noise variance is estimated via the unbiased
  estimator [13]_ defined as:

  .. math::
      \sigma^2 = \frac{\sum_{i=1}^{n} (y_i - \hat{y}_i)^2}{n - p}

  where :math:`p` is the number of features and :math:`\hat{y}_i` is the
  predicted target using an ordinary least squares regression. Note, that this
  formula is valid only when `n_samples > n_features`.

  .. rubric:: References

  .. [12] :arxiv:`Zou, Hui, Trevor Hastie, and Robert Tibshirani.
          "On the degrees of freedom of the lasso."
          The Annals of Statistics 35.5 (2007): 2173-2192.
          <0712.0881.pdf>`

  .. [13] :doi:`Cherkassky, Vladimir, and Yunqian Ma.
          "Comparison of model selection for regression."
          Neural computation 15.7 (2003): 1691-1714.
          <10.1162/089976603321891864>`

Comparison with the regularization parameter of SVM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The equivalence between ``alpha`` and the regularization parameter of SVM,
``C`` is given by ``alpha = 1 / C`` or ``alpha = 1 / (n_samples * C)``,
depending on the estimator and the exact objective function optimized by the
model.

.. _multi_task_lasso:

Multi-task Lasso
================

The :class:`MultiTaskLasso` is a linear model that estimates sparse
coefficients for multiple regression problems jointly: ``y`` is a 2D array,
of shape ``(n_samples, n_tasks)``. The constraint is that the selected
features are the same for all the regression problems, also called tasks.

The following figure compares the location of the non-zero entries in the
coefficient matrix W obtained with a simple Lasso or a MultiTaskLasso.
The Lasso estimates yield scattered non-zeros while the non-zeros of
the MultiTaskLasso are full columns.

.. |multi_task_lasso_1| image:: ../auto_examples/linear_model/images/sphx_glr_plot_multi_task_lasso_support_001.png
    :target: ../auto_examples/linear_model/plot_multi_task_lasso_support.html
    :scale: 48%

.. |multi_task_lasso_2| image:: ../auto_examples/linear_model/images/sphx_glr_plot_multi_task_lasso_support_002.png
    :target: ../auto_examples/linear_model/plot_multi_task_lasso_support.html
    :scale: 48%

.. centered:: |multi_task_lasso_1| |multi_task_lasso_2|

.. centered:: Fitting a time-series model, imposing that any active feature be active at all times.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_multi_task_lasso_support.py`


.. dropdown:: Mathematical details

  Mathematically, it consists of a linear model trained with a mixed
  :math:`\ell_1` :math:`\ell_2`-norm for regularization.
  The objective function to minimize is:

  .. math::  \min_{W} { \frac{1}{2n_{\text{samples}}} ||X W - Y||_{\text{Fro}} ^ 2 + \alpha ||W||_{21}}

  where :math:`\text{Fro}` indicates the Frobenius norm

  .. math:: ||A||_{\text{Fro}} = \sqrt{\sum_{ij} a_{ij}^2}

  and :math:`\ell_1` :math:`\ell_2` reads

  .. math:: ||A||_{2 1} = \sum_i \sqrt{\sum_j a_{ij}^2}.

  The implementation in the class :class:`MultiTaskLasso` uses
  coordinate descent as the algorithm to fit the coefficients.

.. _elastic_net:

Elastic-Net
===========
:class:`ElasticNet` is a linear regression model trained with both
:math:`\ell_1` and :math:`\ell_2`-norm regularization of the coefficients.
This combination  allows for learning a sparse model where few of
the weights are non-zero like :class:`Lasso`, while still maintaining
the regularization properties of :class:`Ridge`. We control the convex
combination of :math:`\ell_1` and :math:`\ell_2` using the ``l1_ratio``
parameter.

Elastic-net is useful when there are multiple features that are
correlated with one another. Lasso is likely to pick one of these
at random, while elastic-net is likely to pick both.

A practical advantage of trading-off between Lasso and Ridge is that it
allows Elastic-Net to inherit some of Ridge's stability under rotation.

The objective function to minimize is in this case

.. math::

    \min_{w} { \frac{1}{2n_{\text{samples}}} ||X w - y||_2 ^ 2 + \alpha \rho ||w||_1 +
    \frac{\alpha(1-\rho)}{2} ||w||_2 ^ 2}


.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_lasso_lasso_lars_elasticnet_path_002.png
   :target: ../auto_examples/linear_model/plot_lasso_lasso_lars_elasticnet_path.html
   :align: center
   :scale: 50%

The class :class:`ElasticNetCV` can be used to set the parameters
``alpha`` (:math:`\alpha`) and ``l1_ratio`` (:math:`\rho`) by cross-validation.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_and_elasticnet.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_lasso_lars_elasticnet_path.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_elastic_net_precomputed_gram_matrix_with_weighted_samples.py`

.. dropdown:: References

  The following two references explain the iterations
  used in the coordinate descent solver of scikit-learn, as well as
  the duality gap computation used for convergence control.

  * "Regularization Path For Generalized linear Models by Coordinate Descent",
    Friedman, Hastie & Tibshirani, J Stat Softw, 2010 (`Paper
    <https://www.jstatsoft.org/article/view/v033i01/v33i01.pdf>`__).
  * "An Interior-Point Method for Large-Scale L1-Regularized Least Squares,"
    S. J. Kim, K. Koh, M. Lustig, S. Boyd and D. Gorinevsky,
    in IEEE Journal of Selected Topics in Signal Processing, 2007
    (`Paper <https://web.stanford.edu/~boyd/papers/pdf/l1_ls.pdf>`__)

.. _multi_task_elastic_net:

Multi-task Elastic-Net
======================

The :class:`MultiTaskElasticNet` is an elastic-net model that estimates sparse
coefficients for multiple regression problems jointly: ``Y`` is a 2D array
of shape ``(n_samples, n_tasks)``. The constraint is that the selected
features are the same for all the regression problems, also called tasks.

Mathematically, it consists of a linear model trained with a mixed
:math:`\ell_1` :math:`\ell_2`-norm and :math:`\ell_2`-norm for regularization.
The objective function to minimize is:

.. math::

    \min_{W} { \frac{1}{2n_{\text{samples}}} ||X W - Y||_{\text{Fro}}^2 + \alpha \rho ||W||_{2 1} +
    \frac{\alpha(1-\rho)}{2} ||W||_{\text{Fro}}^2}

The implementation in the class :class:`MultiTaskElasticNet` uses coordinate descent as
the algorithm to fit the coefficients.

The class :class:`MultiTaskElasticNetCV` can be used to set the parameters
``alpha`` (:math:`\alpha`) and ``l1_ratio`` (:math:`\rho`) by cross-validation.

.. _least_angle_regression:

Least Angle Regression
======================

Least-angle regression (LARS) is a regression algorithm for
high-dimensional data, developed by Bradley Efron, Trevor Hastie, Iain
Johnstone and Robert Tibshirani. LARS is similar to forward stepwise
regression. At each step, it finds the feature most correlated with the
target. When there are multiple features having equal correlation, instead
of continuing along the same feature, it proceeds in a direction equiangular
between the features.

The advantages of LARS are:

- It is numerically efficient in contexts where the number of features
  is significantly greater than the number of samples.

- It is computationally just as fast as forward selection and has
  the same order of complexity as ordinary least squares.

- It produces a full piecewise linear solution path, which is
  useful in cross-validation or similar attempts to tune the model.

- If two features are almost equally correlated with the target,
  then their coefficients should increase at approximately the same
  rate. The algorithm thus behaves as intuition would expect, and
  also is more stable.

- It is easily modified to produce solutions for other estimators,
  like the Lasso.

The disadvantages of the LARS method include:

- Because LARS is based upon an iterative refitting of the
  residuals, it would appear to be especially sensitive to the
  effects of noise. This problem is discussed in detail by Weisberg
  in the discussion section of the Efron et al. (2004) Annals of
  Statistics article.

The LARS model can be used via the estimator :class:`Lars`, or its
low-level implementation :func:`lars_path` or :func:`lars_path_gram`.


LARS Lasso
==========

:class:`LassoLars` is a lasso model implemented using the LARS
algorithm, and unlike the implementation based on coordinate descent,
this yields the exact solution, which is piecewise linear as a
function of the norm of its coefficients.

.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_lasso_lasso_lars_elasticnet_path_001.png
   :target: ../auto_examples/linear_model/plot_lasso_lasso_lars_elasticnet_path.html
   :align: center
   :scale: 50%

::

   >>> from sklearn import linear_model
   >>> reg = linear_model.LassoLars(alpha=.1)
   >>> reg.fit([[0, 0], [1, 1]], [0, 1])
   LassoLars(alpha=0.1)
   >>> reg.coef_
   array([0.6..., 0.        ])

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_lasso_lars_elasticnet_path.py`

The Lars algorithm provides the full path of the coefficients along
the regularization parameter almost for free, thus a common operation
is to retrieve the path with one of the functions :func:`lars_path`
or :func:`lars_path_gram`.

.. dropdown:: Mathematical formulation

  The algorithm is similar to forward stepwise regression, but instead
  of including features at each step, the estimated coefficients are
  increased in a direction equiangular to each one's correlations with
  the residual.

  Instead of giving a vector result, the LARS solution consists of a
  curve denoting the solution for each value of the :math:`\ell_1` norm of the
  parameter vector. The full coefficients path is stored in the array
  ``coef_path_`` of shape `(n_features, max_features + 1)`. The first
  column is always zero.

  .. rubric:: References

  * Original Algorithm is detailed in the paper `Least Angle Regression
    <https://www-stat.stanford.edu/~hastie/Papers/LARS/LeastAngle_2002.pdf>`_
    by Hastie et al.

.. _omp:

Orthogonal Matching Pursuit (OMP)
=================================
:class:`OrthogonalMatchingPursuit` and :func:`orthogonal_mp` implement the OMP
algorithm for approximating the fit of a linear model with constraints imposed
on the number of non-zero coefficients (ie. the :math:`\ell_0` pseudo-norm).

Being a forward feature selection method like :ref:`least_angle_regression`,
orthogonal matching pursuit can approximate the optimum solution vector with a
fixed number of non-zero elements:

.. math::
    \underset{w}{\operatorname{arg\,min\,}}  ||y - Xw||_2^2 \text{ subject to } ||w||_0 \leq n_{\text{nonzero_coefs}}

Alternatively, orthogonal matching pursuit can target a specific error instead
of a specific number of non-zero coefficients. This can be expressed as:

.. math::
    \underset{w}{\operatorname{arg\,min\,}} ||w||_0 \text{ subject to } ||y-Xw||_2^2 \leq \text{tol}


OMP is based on a greedy algorithm that includes at each step the atom most
highly correlated with the current residual. It is similar to the simpler
matching pursuit (MP) method, but better in that at each iteration, the
residual is recomputed using an orthogonal projection on the space of the
previously chosen dictionary elements.


.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_omp.py`

.. dropdown:: References

  * https://www.cs.technion.ac.il/~ronrubin/Publications/KSVD-OMP-v2.pdf

  * `Matching pursuits with time-frequency dictionaries
    <https://www.di.ens.fr/~mallat/papiers/MallatPursuit93.pdf>`_,
    S. G. Mallat, Z. Zhang,

.. _bayesian_regression:

Bayesian Regression
===================

Bayesian regression techniques can be used to include regularization
parameters in the estimation procedure: the regularization parameter is
not set in a hard sense but tuned to the data at hand.

This can be done by introducing `uninformative priors
<https://en.wikipedia.org/wiki/Non-informative_prior#Uninformative_priors>`__
over the hyper parameters of the model.
The :math:`\ell_{2}` regularization used in :ref:`ridge_regression` is
equivalent to finding a maximum a posteriori estimation under a Gaussian prior
over the coefficients :math:`w` with precision :math:`\lambda^{-1}`.
Instead of setting `\lambda` manually, it is possible to treat it as a random
variable to be estimated from the data.

To obtain a fully probabilistic model, the output :math:`y` is assumed
to be Gaussian distributed around :math:`X w`:

.. math::  p(y|X,w,\alpha) = \mathcal{N}(y|X w,\alpha^{-1})

where :math:`\alpha` is again treated as a random variable that is to be
estimated from the data.

The advantages of Bayesian Regression are:

- It adapts to the data at hand.

- It can be used to include regularization parameters in the
  estimation procedure.

The disadvantages of Bayesian regression include:

- Inference of the model can be time consuming.

.. dropdown:: References

  * A good introduction to Bayesian methods is given in C. Bishop: Pattern
    Recognition and Machine learning

  * Original Algorithm is detailed in the  book `Bayesian learning for neural
    networks` by Radford M. Neal

.. _bayesian_ridge_regression:

Bayesian Ridge Regression
-------------------------

:class:`BayesianRidge` estimates a probabilistic model of the
regression problem as described above.
The prior for the coefficient :math:`w` is given by a spherical Gaussian:

.. math:: p(w|\lambda) =
    \mathcal{N}(w|0,\lambda^{-1}\mathbf{I}_{p})

The priors over :math:`\alpha` and :math:`\lambda` are chosen to be `gamma
distributions <https://en.wikipedia.org/wiki/Gamma_distribution>`__, the
conjugate prior for the precision of the Gaussian. The resulting model is
called *Bayesian Ridge Regression*, and is similar to the classical
:class:`Ridge`.

The parameters :math:`w`, :math:`\alpha` and :math:`\lambda` are estimated
jointly during the fit of the model, the regularization parameters
:math:`\alpha` and :math:`\lambda` being estimated by maximizing the
*log marginal likelihood*. The scikit-learn implementation
is based on the algorithm described in Appendix A of (Tipping, 2001)
where the update of the parameters :math:`\alpha` and :math:`\lambda` is done
as suggested in (MacKay, 1992). The initial value of the maximization procedure
can be set with the hyperparameters ``alpha_init`` and ``lambda_init``.

There are four more hyperparameters, :math:`\alpha_1`, :math:`\alpha_2`,
:math:`\lambda_1` and :math:`\lambda_2` of the gamma prior distributions over
:math:`\alpha` and :math:`\lambda`. These are usually chosen to be
*non-informative*. By default :math:`\alpha_1 = \alpha_2 =  \lambda_1 = \lambda_2 = 10^{-6}`.

Bayesian Ridge Regression is used for regression::

    >>> from sklearn import linear_model
    >>> X = [[0., 0.], [1., 1.], [2., 2.], [3., 3.]]
    >>> Y = [0., 1., 2., 3.]
    >>> reg = linear_model.BayesianRidge()
    >>> reg.fit(X, Y)
    BayesianRidge()

After being fitted, the model can then be used to predict new values::

    >>> reg.predict([[1, 0.]])
    array([0.50000013])

The coefficients :math:`w` of the model can be accessed::

    >>> reg.coef_
    array([0.49999993, 0.49999993])

Due to the Bayesian framework, the weights found are slightly different to the
ones found by :ref:`ordinary_least_squares`. However, Bayesian Ridge Regression
is more robust to ill-posed problems.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_bayesian_ridge_curvefit.py`

.. dropdown:: References

  * Section 3.3 in Christopher M. Bishop: Pattern Recognition and Machine Learning, 2006

  * David J. C. MacKay, `Bayesian Interpolation <https://citeseerx.ist.psu.edu/doc_view/pid/b14c7cc3686e82ba40653c6dff178356a33e5e2c>`_, 1992.

  * Michael E. Tipping, `Sparse Bayesian Learning and the Relevance Vector Machine <https://www.jmlr.org/papers/volume1/tipping01a/tipping01a.pdf>`_, 2001.

.. _automatic_relevance_determination:

Automatic Relevance Determination - ARD
---------------------------------------

The Automatic Relevance Determination (as being implemented in
:class:`ARDRegression`) is a kind of linear model which is very similar to the
`Bayesian Ridge Regression`_, but that leads to sparser coefficients :math:`w`
[1]_ [2]_.

:class:`ARDRegression` poses a different prior over :math:`w`: it drops
the spherical Gaussian distribution for a centered elliptic Gaussian
distribution. This means each coefficient :math:`w_{i}` can itself be drawn from
a Gaussian distribution, centered on zero and with a precision
:math:`\lambda_{i}`:

.. math:: p(w|\lambda) = \mathcal{N}(w|0,A^{-1})

with :math:`A` being a positive definite diagonal matrix and
:math:`\text{diag}(A) = \lambda = \{\lambda_{1},...,\lambda_{p}\}`.

In contrast to the `Bayesian Ridge Regression`_, each coordinate of
:math:`w_{i}` has its own standard deviation :math:`\frac{1}{\lambda_i}`. The
prior over all :math:`\lambda_i` is chosen to be the same gamma distribution
given by the hyperparameters :math:`\lambda_1` and :math:`\lambda_2`.

ARD is also known in the literature as *Sparse Bayesian Learning* and *Relevance
Vector Machine* [3]_ [4]_. For a worked-out comparison between ARD and `Bayesian
Ridge Regression`_, see the example below.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_ard.py`


.. rubric:: References

.. [1] Christopher M. Bishop: Pattern Recognition and Machine Learning, Chapter 7.2.1

.. [2] David Wipf and Srikantan Nagarajan: `A New View of Automatic Relevance Determination <https://papers.nips.cc/paper/3372-a-new-view-of-automatic-relevance-determination.pdf>`_

.. [3] Michael E. Tipping: `Sparse Bayesian Learning and the Relevance Vector Machine <https://www.jmlr.org/papers/volume1/tipping01a/tipping01a.pdf>`_

.. [4] Tristan Fletcher: `Relevance Vector Machines Explained <https://citeseerx.ist.psu.edu/doc_view/pid/3dc9d625404fdfef6eaccc3babddefe4c176abd4>`_

.. _Logistic_regression:

Logistic regression
===================

The logistic regression is implemented in :class:`LogisticRegression`. Despite
its name, it is implemented as a linear model for classification rather than
regression in terms of the scikit-learn/ML nomenclature. The logistic
regression is also known in the literature as logit regression,
maximum-entropy classification (MaxEnt) or the log-linear classifier. In this
model, the probabilities describing the possible outcomes of a single trial
are modeled using a `logistic function
<https://en.wikipedia.org/wiki/Logistic_function>`_.

This implementation can fit binary, One-vs-Rest, or multinomial logistic
regression with optional :math:`\ell_1`, :math:`\ell_2` or Elastic-Net
regularization.

.. note:: **Regularization**

    Regularization is applied by default, which is common in machine
    learning but not in statistics. Another advantage of regularization is
    that it improves numerical stability. No regularization amounts to
    setting C to a very high value.

.. note:: **Logistic Regression as a special case of the Generalized Linear Models (GLM)**

    Logistic regression is a special case of
    :ref:`generalized_linear_models` with a Binomial / Bernoulli conditional
    distribution and a Logit link. The numerical output of the logistic
    regression, which is the predicted probability, can be used as a classifier
    by applying a threshold (by default 0.5) to it. This is how it is
    implemented in scikit-learn, so it expects a categorical target, making
    the Logistic Regression a classifier.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_logistic_l1_l2_sparsity.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_logistic_path.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_logistic_multinomial.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_sparse_logistic_regression_20newsgroups.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_sparse_logistic_regression_mnist.py`
* :ref:`sphx_glr_auto_examples_classification_plot_classification_probability.py`

Binary Case
-----------

For notational ease, we assume that the target :math:`y_i` takes values in the
set :math:`\{0, 1\}` for data point :math:`i`.
Once fitted, the :meth:`~sklearn.linear_model.LogisticRegression.predict_proba`
method of :class:`~sklearn.linear_model.LogisticRegression` predicts
the probability of the positive class :math:`P(y_i=1|X_i)` as

.. math:: \hat{p}(X_i) = \operatorname{expit}(X_i w + w_0) = \frac{1}{1 + \exp(-X_i w - w_0)}.


As an optimization problem, binary
class logistic regression with regularization term :math:`r(w)` minimizes the
following cost function:

.. math::
    :name: regularized-logistic-loss

    \min_{w} \frac{1}{S}\sum_{i=1}^n s_i
    \left(-y_i \log(\hat{p}(X_i)) - (1 - y_i) \log(1 - \hat{p}(X_i))\right)
    + \frac{r(w)}{S C}\,,

where :math:`{s_i}` corresponds to the weights assigned by the user to a
specific training sample (the vector :math:`s` is formed by element-wise
multiplication of the class weights and sample weights),
and the sum :math:`S = \sum_{i=1}^n s_i`.

We currently provide four choices for the regularization term  :math:`r(w)` via
the `penalty` argument:

+----------------+-------------------------------------------------+
| penalty        | :math:`r(w)`                                    |
+================+=================================================+
| `None`         | :math:`0`                                       |
+----------------+-------------------------------------------------+
| :math:`\ell_1` | :math:`\|w\|_1`                                 |
+----------------+-------------------------------------------------+
| :math:`\ell_2` | :math:`\frac{1}{2}\|w\|_2^2 = \frac{1}{2}w^T w` |
+----------------+-------------------------------------------------+
| `ElasticNet`   | :math:`\frac{1 - \rho}{2}w^T w + \rho \|w\|_1`  |
+----------------+-------------------------------------------------+

For ElasticNet, :math:`\rho` (which corresponds to the `l1_ratio` parameter)
controls the strength of :math:`\ell_1` regularization vs. :math:`\ell_2`
regularization. Elastic-Net is equivalent to :math:`\ell_1` when
:math:`\rho = 1` and equivalent to :math:`\ell_2` when :math:`\rho=0`.

Note that the scale of the class weights and the sample weights will influence
the optimization problem. For instance, multiplying the sample weights by a
constant :math:`b>0` is equivalent to multiplying the (inverse) regularization
strength `C` by :math:`b`.

Multinomial Case
----------------

The binary case can be extended to :math:`K` classes leading to the multinomial
logistic regression, see also `log-linear model
<https://en.wikipedia.org/wiki/Multinomial_logistic_regression#As_a_log-linear_model>`_.

.. note::
   It is possible to parameterize a :math:`K`-class classification model
   using only :math:`K-1` weight vectors, leaving one class probability fully
   determined by the other class probabilities by leveraging the fact that all
   class probabilities must sum to one. We deliberately choose to overparameterize the model
   using :math:`K` weight vectors for ease of implementation and to preserve the
   symmetrical inductive bias regarding ordering of classes, see [16]_. This effect becomes
   especially important when using regularization. The choice of overparameterization can be
   detrimental for unpenalized models since then the solution may not be unique, as shown in [16]_.

.. dropdown:: Mathematical details

  Let :math:`y_i \in {1, \ldots, K}` be the label (ordinal) encoded target variable for observation :math:`i`.
  Instead of a single coefficient vector, we now have
  a matrix of coefficients :math:`W` where each row vector :math:`W_k` corresponds to class
  :math:`k`. We aim at predicting the class probabilities :math:`P(y_i=k|X_i)` via
  :meth:`~sklearn.linear_model.LogisticRegression.predict_proba` as:

  .. math:: \hat{p}_k(X_i) = \frac{\exp(X_i W_k + W_{0, k})}{\sum_{l=0}^{K-1} \exp(X_i W_l + W_{0, l})}.

  The objective for the optimization becomes

  .. math::
    \min_W -\frac{1}{S}\sum_{i=1}^n \sum_{k=0}^{K-1} s_{ik} [y_i = k] \log(\hat{p}_k(X_i))
    + \frac{r(W)}{S C}\,,

  where :math:`[P]` represents the Iverson bracket which evaluates to :math:`0`
  if :math:`P` is false, otherwise it evaluates to :math:`1`.

  Again, :math:`s_{ik}` are the weights assigned by the user (multiplication of sample
  weights and class weights) with their sum :math:`S = \sum_{i=1}^n \sum_{k=0}^{K-1} s_{ik}`.

  We currently provide four choices
  for the regularization term :math:`r(W)` via the `penalty` argument, where :math:`m`
  is the number of features:

  +----------------+----------------------------------------------------------------------------------+
  | penalty        | :math:`r(W)`                                                                     |
  +================+==================================================================================+
  | `None`         | :math:`0`                                                                        |
  +----------------+----------------------------------------------------------------------------------+
  | :math:`\ell_1` | :math:`\|W\|_{1,1} = \sum_{i=1}^m\sum_{j=1}^{K}|W_{i,j}|`                        |
  +----------------+----------------------------------------------------------------------------------+
  | :math:`\ell_2` | :math:`\frac{1}{2}\|W\|_F^2 = \frac{1}{2}\sum_{i=1}^m\sum_{j=1}^{K} W_{i,j}^2`   |
  +----------------+----------------------------------------------------------------------------------+
  | `ElasticNet`   | :math:`\frac{1 - \rho}{2}\|W\|_F^2 + \rho \|W\|_{1,1}`                           |
  +----------------+----------------------------------------------------------------------------------+

.. _logistic_regression_solvers:

Solvers
-------

The solvers implemented in the class :class:`LogisticRegression`
are "lbfgs", "liblinear", "newton-cg", "newton-cholesky", "sag" and "saga":

The following table summarizes the penalties and multinomial multiclass supported by each solver:

+------------------------------+-----------------+-------------+-----------------+-----------------------+-----------+------------+
|                              |                       **Solvers**                                                                |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| **Penalties**                | **'lbfgs'** | **'liblinear'** | **'newton-cg'** | **'newton-cholesky'** | **'sag'** | **'saga'** |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| L2 penalty                   |     yes     |       no        |       yes       |     no                |    yes    |    yes     |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| L1 penalty                   |     no      |       yes       |       no        |     no                |    no     |    yes     |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| Elastic-Net (L1 + L2)        |     no      |       no        |       no        |     no                |    no     |    yes     |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| No penalty ('none')          |     yes     |       no        |       yes       |     yes               |    yes    |    yes     |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| **Multiclass support**       |                                                                                                  |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| multinomial multiclass       |     yes     |       no        |       yes       |     no                |    yes    |    yes     |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| **Behaviors**                |                                                                                                  |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| Penalize the intercept (bad) |     no      |       yes       |       no        |     no                |    no     |    no      |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| Faster for large datasets    |     no      |       no        |       no        |     no                |    yes    |    yes     |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+
| Robust to unscaled datasets  |     yes     |       yes       |       yes       |     yes               |    no     |    no      |
+------------------------------+-------------+-----------------+-----------------+-----------------------+-----------+------------+

The "lbfgs" solver is used by default for its robustness. For large datasets
the "saga" solver is usually faster.
For large dataset, you may also consider using :class:`SGDClassifier`
with `loss="log_loss"`, which might be even faster but requires more tuning.

.. _liblinear_differences:

Differences between solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There might be a difference in the scores obtained between
:class:`LogisticRegression` with ``solver=liblinear`` or
:class:`~sklearn.svm.LinearSVC` and the external liblinear library directly,
when ``fit_intercept=False`` and the fit ``coef_`` (or) the data to be predicted
are zeroes. This is because for the sample(s) with ``decision_function`` zero,
:class:`LogisticRegression` and :class:`~sklearn.svm.LinearSVC` predict the
negative class, while liblinear predicts the positive class. Note that a model
with ``fit_intercept=False`` and having many samples with ``decision_function``
zero, is likely to be a underfit, bad model and you are advised to set
``fit_intercept=True`` and increase the ``intercept_scaling``.

.. dropdown:: Solvers' details

  * The solver "liblinear" uses a coordinate descent (CD) algorithm, and relies
    on the excellent C++ `LIBLINEAR library
    <https://www.csie.ntu.edu.tw/~cjlin/liblinear/>`_, which is shipped with
    scikit-learn. However, the CD algorithm implemented in liblinear cannot learn
    a true multinomial (multiclass) model; instead, the optimization problem is
    decomposed in a "one-vs-rest" fashion so separate binary classifiers are
    trained for all classes. This happens under the hood, so
    :class:`LogisticRegression` instances using this solver behave as multiclass
    classifiers. For :math:`\ell_1` regularization :func:`sklearn.svm.l1_min_c` allows to
    calculate the lower bound for C in order to get a non "null" (all feature
    weights to zero) model.

  * The "lbfgs", "newton-cg" and "sag" solvers only support :math:`\ell_2`
    regularization or no regularization, and are found to converge faster for some
    high-dimensional data. Setting `multi_class` to "multinomial" with these solvers
    learns a true multinomial logistic regression model [5]_, which means that its
    probability estimates should be better calibrated than the default "one-vs-rest"
    setting.

  * The "sag" solver uses Stochastic Average Gradient descent [6]_. It is faster
    than other solvers for large datasets, when both the number of samples and the
    number of features are large.

  * The "saga" solver [7]_ is a variant of "sag" that also supports the
    non-smooth `penalty="l1"`. This is therefore the solver of choice for sparse
    multinomial logistic regression. It is also the only solver that supports
    `penalty="elasticnet"`.

  * The "lbfgs" is an optimization algorithm that approximates the
    Broyden–Fletcher–Goldfarb–Shanno algorithm [8]_, which belongs to
    quasi-Newton methods. As such, it can deal with a wide range of different training
    data and is therefore the default solver. Its performance, however, suffers on poorly
    scaled datasets and on datasets with one-hot encoded categorical features with rare
    categories.

  * The "newton-cholesky" solver is an exact Newton solver that calculates the hessian
    matrix and solves the resulting linear system. It is a very good choice for
    `n_samples` >> `n_features`, but has a few shortcomings: Only :math:`\ell_2`
    regularization is supported. Furthermore, because the hessian matrix is explicitly
    computed, the memory usage has a quadratic dependency on `n_features` as well as on
    `n_classes`. As a consequence, only the one-vs-rest scheme is implemented for the
    multiclass case.

  For a comparison of some of these solvers, see [9]_.

  .. rubric:: References

  .. [5] Christopher M. Bishop: Pattern Recognition and Machine Learning, Chapter 4.3.4

  .. [6] Mark Schmidt, Nicolas Le Roux, and Francis Bach: `Minimizing Finite Sums with the Stochastic Average Gradient. <https://hal.inria.fr/hal-00860051/document>`_

  .. [7] Aaron Defazio, Francis Bach, Simon Lacoste-Julien:
      :arxiv:`SAGA: A Fast Incremental Gradient Method With Support for
      Non-Strongly Convex Composite Objectives. <1407.0202>`

  .. [8] https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm

  .. [9] Thomas P. Minka `"A comparison of numerical optimizers for logistic regression"
          <https://tminka.github.io/papers/logreg/minka-logreg.pdf>`_

  .. [16] :arxiv:`Simon, Noah, J. Friedman and T. Hastie.
      "A Blockwise Descent Algorithm for Group-penalized Multiresponse and
      Multinomial Regression." <1311.6529>`


.. note:: **Feature selection with sparse logistic regression**

   A logistic regression with :math:`\ell_1` penalty yields sparse models, and can
   thus be used to perform feature selection, as detailed in
   :ref:`l1_feature_selection`.

.. note:: **P-value estimation**

    It is possible to obtain the p-values and confidence intervals for
    coefficients in cases of regression without penalization. The `statsmodels
    package <https://pypi.org/project/statsmodels/>`_ natively supports this.
    Within sklearn, one could use bootstrapping instead as well.


:class:`LogisticRegressionCV` implements Logistic Regression with built-in
cross-validation support, to find the optimal `C` and `l1_ratio` parameters
according to the ``scoring`` attribute. The "newton-cg", "sag", "saga" and
"lbfgs" solvers are found to be faster for high-dimensional dense data, due
to warm-starting (see :term:`Glossary <warm_start>`).

.. _Generalized_linear_regression:

.. _Generalized_linear_models:

Generalized Linear Models
=========================

Generalized Linear Models (GLM) extend linear models in two ways
[10]_. First, the predicted values :math:`\hat{y}` are linked to a linear
combination of the input variables :math:`X` via an inverse link function
:math:`h` as

.. math::    \hat{y}(w, X) = h(Xw).

Secondly, the squared loss function is replaced by the unit deviance
:math:`d` of a distribution in the exponential family (or more precisely, a
reproductive exponential dispersion model (EDM) [11]_).

The minimization problem becomes:

.. math::    \min_{w} \frac{1}{2 n_{\text{samples}}} \sum_i d(y_i, \hat{y}_i) + \frac{\alpha}{2} ||w||_2^2,

where :math:`\alpha` is the L2 regularization penalty. When sample weights are
provided, the average becomes a weighted average.

The following table lists some specific EDMs and their unit deviance :

================= ================================  ============================================
Distribution       Target Domain                    Unit Deviance :math:`d(y, \hat{y})`
================= ================================  ============================================
Normal            :math:`y \in (-\infty, \infty)`   :math:`(y-\hat{y})^2`
Bernoulli         :math:`y \in \{0, 1\}`            :math:`2({y}\log\frac{y}{\hat{y}}+({1}-{y})\log\frac{{1}-{y}}{{1}-\hat{y}})`
Categorical       :math:`y \in \{0, 1, ..., k\}`    :math:`2\sum_{i \in \{0, 1, ..., k\}} I(y = i) y_\text{i}\log\frac{I(y = i)}{\hat{I(y = i)}}`
Poisson           :math:`y \in [0, \infty)`         :math:`2(y\log\frac{y}{\hat{y}}-y+\hat{y})`
Gamma             :math:`y \in (0, \infty)`         :math:`2(\log\frac{\hat{y}}{y}+\frac{y}{\hat{y}}-1)`
Inverse Gaussian  :math:`y \in (0, \infty)`         :math:`\frac{(y-\hat{y})^2}{y\hat{y}^2}`
================= ================================  ============================================

The Probability Density Functions (PDF) of these distributions are illustrated
in the following figure,

.. figure:: ./glm_data/poisson_gamma_tweedie_distributions.png
   :align: center
   :scale: 100%

   PDF of a random variable Y following Poisson, Tweedie (power=1.5) and Gamma
   distributions with different mean values (:math:`\mu`). Observe the point
   mass at :math:`Y=0` for the Poisson distribution and the Tweedie (power=1.5)
   distribution, but not for the Gamma distribution which has a strictly
   positive target domain.

The Bernoulli distribution is a discrete probability distribution modelling a
Bernoulli trial - an event that has only two mutually exclusive outcomes.
The Categorical distribution is a generalization of the Bernoulli distribution
for a categorical random variable. While a random variable in a Bernoulli
distribution has two possible outcomes, a Categorical random variable can take
on one of K possible categories, with the probability of each category
specified separately.

The choice of the distribution depends on the problem at hand:

* If the target values :math:`y` are counts (non-negative integer valued) or
  relative frequencies (non-negative), you might use a Poisson distribution
  with a log-link.
* If the target values are positive valued and skewed, you might try a Gamma
  distribution with a log-link.
* If the target values seem to be heavier tailed than a Gamma distribution, you
  might try an Inverse Gaussian distribution (or even higher variance powers of
  the Tweedie family).
* If the target values :math:`y` are probabilities, you can use the Bernoulli
  distribution. The Bernoulli distribution with a logit link can be used for
  binary classification. The Categorical distribution with a softmax link can be
  used for multiclass classification.


.. dropdown:: Examples of use cases

  * Agriculture / weather modeling:  number of rain events per year (Poisson),
    amount of rainfall per event (Gamma), total rainfall per year (Tweedie /
    Compound Poisson Gamma).
  * Risk modeling / insurance policy pricing:  number of claim events /
    policyholder per year (Poisson), cost per event (Gamma), total cost per
    policyholder per year (Tweedie / Compound Poisson Gamma).
  * Credit Default: probability that a loan can't be paid back (Bernoulli).
  * Fraud Detection: probability that a financial transaction like a cash transfer
    is a fraudulent transaction (Bernoulli).
  * Predictive maintenance: number of production interruption events per year
    (Poisson), duration of interruption (Gamma), total interruption time per year
    (Tweedie / Compound Poisson Gamma).
  * Medical Drug Testing: probability of curing a patient in a set of trials or
    probability that a patient will experience side effects (Bernoulli).
  * News Classification: classification of news articles into three categories
    namely Business News, Politics and Entertainment news (Categorical).

.. rubric:: References

.. [10] McCullagh, Peter; Nelder, John (1989). Generalized Linear Models,
    Second Edition. Boca Raton: Chapman and Hall/CRC. ISBN 0-412-31760-5.

.. [11] Jørgensen, B. (1992). The theory of exponential dispersion models
    and analysis of deviance. Monografias de matemática, no. 51.  See also
    `Exponential dispersion model.
    <https://en.wikipedia.org/wiki/Exponential_dispersion_model>`_

Usage
-----

:class:`TweedieRegressor` implements a generalized linear model for the
Tweedie distribution, that allows to model any of the above mentioned
distributions using the appropriate ``power`` parameter. In particular:

- ``power = 0``: Normal distribution. Specific estimators such as
  :class:`Ridge`, :class:`ElasticNet` are generally more appropriate in
  this case.
- ``power = 1``: Poisson distribution. :class:`PoissonRegressor` is exposed
  for convenience. However, it is strictly equivalent to
  `TweedieRegressor(power=1, link='log')`.
- ``power = 2``: Gamma distribution. :class:`GammaRegressor` is exposed for
  convenience. However, it is strictly equivalent to
  `TweedieRegressor(power=2, link='log')`.
- ``power = 3``: Inverse Gaussian distribution.

The link function is determined by the `link` parameter.

Usage example::

    >>> from sklearn.linear_model import TweedieRegressor
    >>> reg = TweedieRegressor(power=1, alpha=0.5, link='log')
    >>> reg.fit([[0, 0], [0, 1], [2, 2]], [0, 1, 2])
    TweedieRegressor(alpha=0.5, link='log', power=1)
    >>> reg.coef_
    array([0.2463..., 0.4337...])
    >>> reg.intercept_
    -0.7638...


.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_poisson_regression_non_normal_loss.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_tweedie_regression_insurance_claims.py`

.. dropdown:: Practical considerations

  The feature matrix `X` should be standardized before fitting. This ensures
  that the penalty treats features equally.

  Since the linear predictor :math:`Xw` can be negative and Poisson,
  Gamma and Inverse Gaussian distributions don't support negative values, it
  is necessary to apply an inverse link function that guarantees the
  non-negativeness. For example with `link='log'`, the inverse link function
  becomes :math:`h(Xw)=\exp(Xw)`.

  If you want to model a relative frequency, i.e. counts per exposure (time,
  volume, ...) you can do so by using a Poisson distribution and passing
  :math:`y=\frac{\mathrm{counts}}{\mathrm{exposure}}` as target values
  together with :math:`\mathrm{exposure}` as sample weights. For a concrete
  example see e.g.
  :ref:`sphx_glr_auto_examples_linear_model_plot_tweedie_regression_insurance_claims.py`.

  When performing cross-validation for the `power` parameter of
  `TweedieRegressor`, it is advisable to specify an explicit `scoring` function,
  because the default scorer :meth:`TweedieRegressor.score` is a function of
  `power` itself.

Stochastic Gradient Descent - SGD
=================================

Stochastic gradient descent is a simple yet very efficient approach
to fit linear models. It is particularly useful when the number of samples
(and the number of features) is very large.
The ``partial_fit`` method allows online/out-of-core learning.

The classes :class:`SGDClassifier` and :class:`SGDRegressor` provide
functionality to fit linear models for classification and regression
using different (convex) loss functions and different penalties.
E.g., with ``loss="log"``, :class:`SGDClassifier`
fits a logistic regression model,
while with ``loss="hinge"`` it fits a linear support vector machine (SVM).

You can refer to the dedicated :ref:`sgd` documentation section for more details.

.. _perceptron:

Perceptron
==========

The :class:`Perceptron` is another simple classification algorithm suitable for
large scale learning. By default:

- It does not require a learning rate.

- It is not regularized (penalized).

- It updates its model only on mistakes.

The last characteristic implies that the Perceptron is slightly faster to
train than SGD with the hinge loss and that the resulting models are
sparser.

In fact, the :class:`Perceptron` is a wrapper around the :class:`SGDClassifier`
class using a perceptron loss and a constant learning rate. Refer to
:ref:`mathematical section <sgd_mathematical_formulation>` of the SGD procedure
for more details.

.. _passive_aggressive:

Passive Aggressive Algorithms
=============================

The passive-aggressive algorithms are a family of algorithms for large-scale
learning. They are similar to the Perceptron in that they do not require a
learning rate. However, contrary to the Perceptron, they include a
regularization parameter ``C``.

For classification, :class:`PassiveAggressiveClassifier` can be used with
``loss='hinge'`` (PA-I) or ``loss='squared_hinge'`` (PA-II).  For regression,
:class:`PassiveAggressiveRegressor` can be used with
``loss='epsilon_insensitive'`` (PA-I) or
``loss='squared_epsilon_insensitive'`` (PA-II).

.. dropdown:: References

  * `"Online Passive-Aggressive Algorithms"
    <http://jmlr.csail.mit.edu/papers/volume7/crammer06a/crammer06a.pdf>`_
    K. Crammer, O. Dekel, J. Keshat, S. Shalev-Shwartz, Y. Singer - JMLR 7 (2006)

Robustness regression: outliers and modeling errors
=====================================================

Robust regression aims to fit a regression model in the
presence of corrupt data: either outliers, or error in the model.

.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_theilsen_001.png
   :target: ../auto_examples/linear_model/plot_theilsen.html
   :scale: 50%
   :align: center

Different scenario and useful concepts
----------------------------------------

There are different things to keep in mind when dealing with data
corrupted by outliers:

.. |y_outliers| image:: ../auto_examples/linear_model/images/sphx_glr_plot_robust_fit_003.png
   :target: ../auto_examples/linear_model/plot_robust_fit.html
   :scale: 60%

.. |X_outliers| image:: ../auto_examples/linear_model/images/sphx_glr_plot_robust_fit_002.png
   :target: ../auto_examples/linear_model/plot_robust_fit.html
   :scale: 60%

.. |large_y_outliers| image:: ../auto_examples/linear_model/images/sphx_glr_plot_robust_fit_005.png
   :target: ../auto_examples/linear_model/plot_robust_fit.html
   :scale: 60%

* **Outliers in X or in y**?

  ==================================== ====================================
  Outliers in the y direction          Outliers in the X direction
  ==================================== ====================================
  |y_outliers|                         |X_outliers|
  ==================================== ====================================

* **Fraction of outliers versus amplitude of error**

  The number of outlying points matters, but also how much they are
  outliers.

  ==================================== ====================================
  Small outliers                       Large outliers
  ==================================== ====================================
  |y_outliers|                         |large_y_outliers|
  ==================================== ====================================

An important notion of robust fitting is that of breakdown point: the
fraction of data that can be outlying for the fit to start missing the
inlying data.

Note that in general, robust fitting in high-dimensional setting (large
`n_features`) is very hard. The robust models here will probably not work
in these settings.


.. topic:: Trade-offs: which estimator ?

  Scikit-learn provides 3 robust regression estimators:
  :ref:`RANSAC <ransac_regression>`,
  :ref:`Theil Sen <theil_sen_regression>` and
  :ref:`HuberRegressor <huber_regression>`.

  * :ref:`HuberRegressor <huber_regression>` should be faster than
    :ref:`RANSAC <ransac_regression>` and :ref:`Theil Sen <theil_sen_regression>`
    unless the number of samples are very large, i.e. ``n_samples`` >> ``n_features``.
    This is because :ref:`RANSAC <ransac_regression>` and :ref:`Theil Sen <theil_sen_regression>`
    fit on smaller subsets of the data. However, both :ref:`Theil Sen <theil_sen_regression>`
    and :ref:`RANSAC <ransac_regression>` are unlikely to be as robust as
    :ref:`HuberRegressor <huber_regression>` for the default parameters.

  * :ref:`RANSAC <ransac_regression>` is faster than :ref:`Theil Sen <theil_sen_regression>`
    and scales much better with the number of samples.

  * :ref:`RANSAC <ransac_regression>` will deal better with large
    outliers in the y direction (most common situation).

  * :ref:`Theil Sen <theil_sen_regression>` will cope better with
    medium-size outliers in the X direction, but this property will
    disappear in high-dimensional settings.

  When in doubt, use :ref:`RANSAC <ransac_regression>`.

.. _ransac_regression:

RANSAC: RANdom SAmple Consensus
--------------------------------

RANSAC (RANdom SAmple Consensus) fits a model from random subsets of
inliers from the complete data set.

RANSAC is a non-deterministic algorithm producing only a reasonable result with
a certain probability, which is dependent on the number of iterations (see
`max_trials` parameter). It is typically used for linear and non-linear
regression problems and is especially popular in the field of photogrammetric
computer vision.

The algorithm splits the complete input sample data into a set of inliers,
which may be subject to noise, and outliers, which are e.g. caused by erroneous
measurements or invalid hypotheses about the data. The resulting model is then
estimated only from the determined inliers.

.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_ransac_001.png
   :target: ../auto_examples/linear_model/plot_ransac.html
   :align: center
   :scale: 50%

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_ransac.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_robust_fit.py`

.. dropdown:: Details of the algorithm

  Each iteration performs the following steps:

  1. Select ``min_samples`` random samples from the original data and check
     whether the set of data is valid (see ``is_data_valid``).
  2. Fit a model to the random subset (``estimator.fit``) and check
     whether the estimated model is valid (see ``is_model_valid``).
  3. Classify all data as inliers or outliers by calculating the residuals
     to the estimated model (``estimator.predict(X) - y``) - all data
     samples with absolute residuals smaller than or equal to the
     ``residual_threshold`` are considered as inliers.
  4. Save fitted model as best model if number of inlier samples is
     maximal. In case the current estimated model has the same number of
     inliers, it is only considered as the best model if it has better score.

  These steps are performed either a maximum number of times (``max_trials``) or
  until one of the special stop criteria are met (see ``stop_n_inliers`` and
  ``stop_score``). The final model is estimated using all inlier samples (consensus
  set) of the previously determined best model.

  The ``is_data_valid`` and ``is_model_valid`` functions allow to identify and reject
  degenerate combinations of random sub-samples. If the estimated model is not
  needed for identifying degenerate cases, ``is_data_valid`` should be used as it
  is called prior to fitting the model and thus leading to better computational
  performance.

.. dropdown:: References

  * https://en.wikipedia.org/wiki/RANSAC
  * `"Random Sample Consensus: A Paradigm for Model Fitting with Applications to
    Image Analysis and Automated Cartography"
    <https://www.cs.ait.ac.th/~mdailey/cvreadings/Fischler-RANSAC.pdf>`_
    Martin A. Fischler and Robert C. Bolles - SRI International (1981)
  * `"Performance Evaluation of RANSAC Family"
    <http://www.bmva.org/bmvc/2009/Papers/Paper355/Paper355.pdf>`_
    Sunglok Choi, Taemin Kim and Wonpil Yu - BMVC (2009)

.. _theil_sen_regression:

Theil-Sen estimator: generalized-median-based estimator
--------------------------------------------------------

The :class:`TheilSenRegressor` estimator uses a generalization of the median in
multiple dimensions. It is thus robust to multivariate outliers. Note however
that the robustness of the estimator decreases quickly with the dimensionality
of the problem. It loses its robustness properties and becomes no
better than an ordinary least squares in high dimension.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_theilsen.py`
* :ref:`sphx_glr_auto_examples_linear_model_plot_robust_fit.py`


.. dropdown:: Theoretical considerations

  :class:`TheilSenRegressor` is comparable to the :ref:`Ordinary Least Squares
  (OLS) <ordinary_least_squares>` in terms of asymptotic efficiency and as an
  unbiased estimator. In contrast to OLS, Theil-Sen is a non-parametric
  method which means it makes no assumption about the underlying
  distribution of the data. Since Theil-Sen is a median-based estimator, it
  is more robust against corrupted data aka outliers. In univariate
  setting, Theil-Sen has a breakdown point of about 29.3% in case of a
  simple linear regression which means that it can tolerate arbitrary
  corrupted data of up to 29.3%.

  .. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_theilsen_001.png
    :target: ../auto_examples/linear_model/plot_theilsen.html
    :align: center
    :scale: 50%

  The implementation of :class:`TheilSenRegressor` in scikit-learn follows a
  generalization to a multivariate linear regression model [#f1]_ using the
  spatial median which is a generalization of the median to multiple
  dimensions [#f2]_.

  In terms of time and space complexity, Theil-Sen scales according to

  .. math::
      \binom{n_{\text{samples}}}{n_{\text{subsamples}}}

  which makes it infeasible to be applied exhaustively to problems with a
  large number of samples and features. Therefore, the magnitude of a
  subpopulation can be chosen to limit the time and space complexity by
  considering only a random subset of all possible combinations.

  .. rubric:: References

  .. [#f1] Xin Dang, Hanxiang Peng, Xueqin Wang and Heping Zhang: `Theil-Sen Estimators in a Multiple Linear Regression Model. <http://home.olemiss.edu/~xdang/papers/MTSE.pdf>`_

  .. [#f2] T. Kärkkäinen and S. Äyrämö: `On Computation of Spatial Median for Robust Data Mining. <http://users.jyu.fi/~samiayr/pdf/ayramo_eurogen05.pdf>`_

  Also see the `Wikipedia page <https://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator>`_


.. _huber_regression:

Huber Regression
----------------

The :class:`HuberRegressor` is different from :class:`Ridge` because it applies a
linear loss to samples that are defined as outliers by the `epsilon` parameter.
A sample is classified as an inlier if the absolute error of that sample is
lesser than the threshold `epsilon`. It differs from :class:`TheilSenRegressor`
and :class:`RANSACRegressor` because it does not ignore the effect of the outliers
but gives a lesser weight to them.

.. figure:: /auto_examples/linear_model/images/sphx_glr_plot_huber_vs_ridge_001.png
   :target: ../auto_examples/linear_model/plot_huber_vs_ridge.html
   :align: center
   :scale: 50%

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_huber_vs_ridge.py`

.. dropdown:: Mathematical details

  :class:`HuberRegressor` minimizes

  .. math::

    \min_{w, \sigma} {\sum_{i=1}^n\left(\sigma + H_{\epsilon}\left(\frac{X_{i}w - y_{i}}{\sigma}\right)\sigma\right) + \alpha {||w||_2}^2}

  where the loss function is given by

  .. math::

    H_{\epsilon}(z) = \begin{cases}
          z^2, & \text {if } |z| < \epsilon, \\
          2\epsilon|z| - \epsilon^2, & \text{otherwise}
    \end{cases}

  It is advised to set the parameter ``epsilon`` to 1.35 to achieve 95%
  statistical efficiency.

  .. rubric:: References

  * Peter J. Huber, Elvezio M. Ronchetti: Robust Statistics, Concomitant scale
    estimates, p. 172.

The :class:`HuberRegressor` differs from using :class:`SGDRegressor` with loss set to `huber`
in the following ways.

- :class:`HuberRegressor` is scaling invariant. Once ``epsilon`` is set, scaling ``X`` and ``y``
  down or up by different values would produce the same robustness to outliers as before.
  as compared to :class:`SGDRegressor` where ``epsilon`` has to be set again when ``X`` and ``y`` are
  scaled.

- :class:`HuberRegressor` should be more efficient to use on data with small number of
  samples while :class:`SGDRegressor` needs a number of passes on the training data to
  produce the same robustness.

Note that this estimator is different from the `R implementation of Robust
Regression <https://stats.oarc.ucla.edu/r/dae/robust-regression/>`_  because the R
implementation does a weighted least squares implementation with weights given to each
sample on the basis of how much the residual is greater than a certain threshold.

.. _quantile_regression:

Quantile Regression
===================

Quantile regression estimates the median or other quantiles of :math:`y`
conditional on :math:`X`, while ordinary least squares (OLS) estimates the
conditional mean.

Quantile regression may be useful if one is interested in predicting an
interval instead of point prediction. Sometimes, prediction intervals are
calculated based on the assumption that prediction error is distributed
normally with zero mean and constant variance. Quantile regression provides
sensible prediction intervals even for errors with non-constant (but
predictable) variance or non-normal distribution.

.. figure:: /auto_examples/linear_model/images/sphx_glr_plot_quantile_regression_002.png
   :target: ../auto_examples/linear_model/plot_quantile_regression.html
   :align: center
   :scale: 50%

Based on minimizing the pinball loss, conditional quantiles can also be
estimated by models other than linear models. For example,
:class:`~sklearn.ensemble.GradientBoostingRegressor` can predict conditional
quantiles if its parameter ``loss`` is set to ``"quantile"`` and parameter
``alpha`` is set to the quantile that should be predicted. See the example in
:ref:`sphx_glr_auto_examples_ensemble_plot_gradient_boosting_quantile.py`.

Most implementations of quantile regression are based on linear programming
problem. The current implementation is based on
:func:`scipy.optimize.linprog`.

.. rubric:: Examples

* :ref:`sphx_glr_auto_examples_linear_model_plot_quantile_regression.py`

.. dropdown:: Mathematical details

  As a linear model, the :class:`QuantileRegressor` gives linear predictions
  :math:`\hat{y}(w, X) = Xw` for the :math:`q`-th quantile, :math:`q \in (0, 1)`.
  The weights or coefficients :math:`w` are then found by the following
  minimization problem:

  .. math::
      \min_{w} {\frac{1}{n_{\text{samples}}}
      \sum_i PB_q(y_i - X_i w) + \alpha ||w||_1}.

  This consists of the pinball loss (also known as linear loss),
  see also :class:`~sklearn.metrics.mean_pinball_loss`,

  .. math::
      PB_q(t) = q \max(t, 0) + (1 - q) \max(-t, 0) =
      \begin{cases}
          q t, & t > 0, \\
          0,    & t = 0, \\
          (q-1) t, & t < 0
      \end{cases}

  and the L1 penalty controlled by parameter ``alpha``, similar to
  :class:`Lasso`.

  As the pinball loss is only linear in the residuals, quantile regression is
  much more robust to outliers than squared error based estimation of the mean.
  Somewhat in between is the :class:`HuberRegressor`.

.. dropdown:: References

  * Koenker, R., & Bassett Jr, G. (1978). `Regression quantiles.
    <https://gib.people.uic.edu/RQ.pdf>`_
    Econometrica: journal of the Econometric Society, 33-50.

  * Portnoy, S., & Koenker, R. (1997). :doi:`The Gaussian hare and the Laplacian
    tortoise: computability of squared-error versus absolute-error estimators.
    Statistical Science, 12, 279-300 <10.1214/ss/1030037960>`.

  * Koenker, R. (2005). :doi:`Quantile Regression <10.1017/CBO9780511754098>`.
    Cambridge University Press.


.. _polynomial_regression:

Polynomial regression: extending linear models with basis functions
===================================================================

.. currentmodule:: sklearn.preprocessing

One common pattern within machine learning is to use linear models trained
on nonlinear functions of the data.  This approach maintains the generally
fast performance of linear methods, while allowing them to fit a much wider
range of data.

.. dropdown:: Mathematical details

  For example, a simple linear regression can be extended by constructing
  **polynomial features** from the coefficients.  In the standard linear
  regression case, you might have a model that looks like this for
  two-dimensional data:

  .. math::    \hat{y}(w, x) = w_0 + w_1 x_1 + w_2 x_2

  If we want to fit a paraboloid to the data instead of a plane, we can combine
  the features in second-order polynomials, so that the model looks like this:

  .. math::    \hat{y}(w, x) = w_0 + w_1 x_1 + w_2 x_2 + w_3 x_1 x_2 + w_4 x_1^2 + w_5 x_2^2

  The (sometimes surprising) observation is that this is *still a linear model*:
  to see this, imagine creating a new set of features

  .. math::  z = [x_1, x_2, x_1 x_2, x_1^2, x_2^2]

  With this re-labeling of the data, our problem can be written

  .. math::    \hat{y}(w, z) = w_0 + w_1 z_1 + w_2 z_2 + w_3 z_3 + w_4 z_4 + w_5 z_5

  We see that the resulting *polynomial regression* is in the same class of
  linear models we considered above (i.e. the model is linear in :math:`w`)
  and can be solved by the same techniques.  By considering linear fits within
  a higher-dimensional space built with these basis functions, the model has the
  flexibility to fit a much broader range of data.

Here is an example of applying this idea to one-dimensional data, using
polynomial features of varying degrees:

.. figure:: ../auto_examples/linear_model/images/sphx_glr_plot_polynomial_interpolation_001.png
   :target: ../auto_examples/linear_model/plot_polynomial_interpolation.html
   :align: center
   :scale: 50%

This figure is created using the :class:`PolynomialFeatures` transformer, which
transforms an input data matrix into a new data matrix of a given degree.
It can be used as follows::

    >>> from sklearn.preprocessing import PolynomialFeatures
    >>> import numpy as np
    >>> X = np.arange(6).reshape(3, 2)
    >>> X
    array([[0, 1],
           [2, 3],
           [4, 5]])
    >>> poly = PolynomialFeatures(degree=2)
    >>> poly.fit_transform(X)
    array([[ 1.,  0.,  1.,  0.,  0.,  1.],
           [ 1.,  2.,  3.,  4.,  6.,  9.],
           [ 1.,  4.,  5., 16., 20., 25.]])

The features of ``X`` have been transformed from :math:`[x_1, x_2]` to
:math:`[1, x_1, x_2, x_1^2, x_1 x_2, x_2^2]`, and can now be used within
any linear model.

This sort of preprocessing can be streamlined with the
:ref:`Pipeline <pipeline>` tools. A single object representing a simple
polynomial regression can be created and used as follows::

    >>> from sklearn.preprocessing import PolynomialFeatures
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.pipeline import Pipeline
    >>> import numpy as np
    >>> model = Pipeline([('poly', PolynomialFeatures(degree=3)),
    ...                   ('linear', LinearRegression(fit_intercept=False))])
    >>> # fit to an order-3 polynomial data
    >>> x = np.arange(5)
    >>> y = 3 - 2 * x + x ** 2 - x ** 3
    >>> model = model.fit(x[:, np.newaxis], y)
    >>> model.named_steps['linear'].coef_
    array([ 3., -2.,  1., -1.])

The linear model trained on polynomial features is able to exactly recover
the input polynomial coefficients.

In some cases it's not necessary to include higher powers of any single feature,
but only the so-called *interaction features*
that multiply together at most :math:`d` distinct features.
These can be gotten from :class:`PolynomialFeatures` with the setting
``interaction_only=True``.

For example, when dealing with boolean features,
:math:`x_i^n = x_i` for all :math:`n` and is therefore useless;
but :math:`x_i x_j` represents the conjunction of two booleans.
This way, we can solve the XOR problem with a linear classifier::

    >>> from sklearn.linear_model import Perceptron
    >>> from sklearn.preprocessing import PolynomialFeatures
    >>> import numpy as np
    >>> X = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    >>> y = X[:, 0] ^ X[:, 1]
    >>> y
    array([0, 1, 1, 0])
    >>> X = PolynomialFeatures(interaction_only=True).fit_transform(X).astype(int)
    >>> X
    array([[1, 0, 0, 0],
           [1, 0, 1, 0],
           [1, 1, 0, 0],
           [1, 1, 1, 1]])
    >>> clf = Perceptron(fit_intercept=False, max_iter=10, tol=None,
    ...                  shuffle=False).fit(X, y)

And the classifier "predictions" are perfect::

    >>> clf.predict(X)
    array([0, 1, 1, 0])
    >>> clf.score(X, y)
    1.0
