===========================
Stochastic Gradient Descent
===========================

.. currentmodule:: scikits.learn.sgd

**Stochastic Gradient Descent (SGD)** is ...
.

The advantages of Stochastic Gradient Descent are:

    - Efficiency

    - Ease of implementation (lots of opportunities for code tuning). 

The disadvantages of Stochastic Gradient Descent include:

    - SGD requires a number of hyperparameters including the regularization
      parameter and the number of iterations. 

    - SGD is sensitive to feature scaling. If your features do not have 
      an intrinsic scale (e.g. after applying PCA), make sure you standardize
      them to zero mean and unit variance. 

    - SVMs do not directly provide probability estimates, so these
      must be calculated using indirect techniques. In our case, these
      techniques imply conducting five-fold cross-validation, so
      performance can suffer.  See method predict_proba for more
      information.


.. _sgd_classification:

Classification
==============



.. currentmodule:: scikits.learn.svm.sparse


Stochastic Gradient Descent for sparse data
===========================================

There is support for sparse data given in any matrix in a format
supported by scipy.sparse. Classes have the same name, just prefixed
by the `sparse` namespace, and take the same arguments, with the
exception of training and test data, which is expected to be in a
matrix format defined in scipy.sparse.

For maximum efficiency, use the CSR matrix format as defined in
`scipy.sparse.csr_matrix
<http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html>`_.

Implemented classes are :class:`SGD`.


Complexity
==========

Support Vector Machines are powerful tools, but their compute and
storage requirements increase rapidly with the number of training
vectors. The core of an SVM is a quadratic programming problem (QP),
separating support vectors from the rest of the training data. The QP
solver used by this `libsvm`_-based implementation scales between
:math:O(n_{features} \times n_{samples}^2) and
:math:O(n_{features} \times n_{samples}^3) depending on how efficiently
the `libsvm`_ cache is used in practice (dataset dependent). If the data
is very sparse :math:n_{features} should be replaced by the average number
of non-zero features in a sample vector

Also note that for the linear case, the algorithm used in
:class:`LinearSVC` by the `liblinear`_ implementation is much more
efficient than its `libsvm`_-based :class:`SVC` counterpart and can
scale almost linearly to millions of samples and/or features.


Tips on Practical Use
=====================

  * Stochastic Gradient Descent is sensitive to feature scaling, so it
    is highly recommended to scale your data. For example, scale each
    attribute on the input vector X to [0,1] or [-1,+1], or standardize
    it to have mean 0 and variance 1. Note that the *same* scaling
    must be applied to the test vector to obtain meaningful
    results. See `The CookBook
    <https://sourceforge.net/apps/trac/scikit-learn/wiki/CookBook>`_
    for some examples on scaling. If your attributes have an inherent
    scale (e.g. word frequencies, indicator features) scaling is
    not needed. 


.. _svm_mathematical_formulation:

Mathematical formulation
========================

A support vector machine constructs a hyper-plane or set of hyper-planes
in a high or infinite dimensional space, which can be used for
classification, regression or other tasks. Intuitively, a good
separation is achieved by the hyper-plane that has the largest distance
to the nearest training data points of any class (so-called functional
margin), since in general the larger the margin the lower the
generalization error of the classifier.


.. figure:: ../auto_examples/svm/images/plot_separating_hyperplane.png
   :align: center
   :scale: 50

SVC
---

Given training vectors :math:`x_i \in R^n`, i=1,..., l, in two
classes, and a vector :math:`y \in R^l` such that :math:`y_i \in {1,
-1}`, SVC solves the following primal problem:


.. math::

    \min_ {w, b, \zeta} \frac{1}{2} w^T w + C \sum_{i=1, l} \zeta_i



    \textrm {subject to } & y_i (w^T \phi (x_i) + b) \geq 1 - \zeta_i,\\
    & \zeta_i \geq 0, i=1, ..., l

Its dual is

.. math::

   \min_{\alpha} \frac{1}{2} \alpha^T Q \alpha - e^T \alpha


   \textrm {subject to } & y^T \alpha = 0\\
   & 0 \leq \alpha_i \leq C, i=1, ..., l

where :math:`e` is the vector of all ones, C > 0 is the upper bound, Q
is an l by l positive semidefinite matrix, :math:`Q_ij \equiv K(x_i,
x_j)` and :math:`\phi (x_i)^T \phi (x)` is the kernel. Here training
vectors are mapped into a higher (maybe infinite) dimensional space by
the function :math:`\phi`


The decision function is:

.. math:: sgn(\sum_{i=1}^l y_i \alpha_i K(x_i, x) + \rho)


.. TODO multiclass case ?/

This parameters can be accessed through the members support\_ and intercept\_:

     - Member support\_ holds the product :math:`y^T \alpha`

     - Member intercept\_ of the classifier holds :math:`-\rho`

.. topic:: References:

 * `"Pegasos: Primal estimated sub-gradient solver for svm" 
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.74.8513>`_
   S. Shalev-Shwartz, Y. Singer, N. Srebro - In Proceedings of ICML '07.

 * `"Stochastic gradient descent training for l1-regularized log-linear models with cumulative penalty"
   <www.aclweb.org/anthology/P/P09/P09-1054.pdf>`_
   Y. Tsuruoka, J. Tsujii, S. Ananiadou -  In Proceedings of the AFNLP/ACL '09.

 * `"Solving large scale linear prediction problems using stochastic gradient descent algorithms"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.58.7377>`_
   T. Zhang - In Proceedings of ICML '04.
   
  * `"Regularization and variable selection via the elastic net"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.124.4696>`_
   H. Zou, T. Hastie - Journal of the Royal Statistical Society Series B, 67 (2), 301-320.


Frequently Asked Questions
==========================

     * Q: Can I get the indices of the support vectors instead of the
       support vectors ?

       A: The underlying C implementation does not provide this
       information.


Implementation details
======================

The implementation of SGD is based on the Stochastic Gradient SVM of Léon Bottou. 
Similar to SvmSGD in `sgd`_ the weight vector is represented as the product of a scalar 
and a vector which allows an efficient weight update in the case of L2 regularization. 
The code is written in Cython.

.. _`sgd`: http://leon.bottou.org/projects/sgd

.. topic:: References:

