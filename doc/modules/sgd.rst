===========================
Stochastic Gradient Descent
===========================

.. currentmodule:: scikits.learn.sgd

**Stochastic Gradient Descent (SGD)** is a simple yet very efficient approach 
to discriminative learning of linear classifiers under convex loss functions 
such as Support Vector Machines and Logistic Regression. 
Even though SGD has been around in the ML community for a long time, 
it has received a considerable amount of attention just recently in the 
context of large-scale learning. 

SGD has been successfully applied to large-scale and sparse machine learning 
problems often encountered in text classification and natural language processing. 
Given that the data is sparse, the classifiers in this module easily scale 
to problems with more than 10^5 training examples and more than 10^4 features. 

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

The major advantage of SGD is its efficiency, which is basically 
linear in the number of training examples. Recent theoretical 
results, however, show that the runtime to get some desired optimization 
accuracy does not increase as the training set size increases. 
In fact, for PEGASOS training indeed decreases as the size of the training set 
increases.

Tips on Practical Use
=====================

  * Stochastic Gradient Descent is sensitive to feature scaling, so it
    is highly recommended to scale your data. For example, scale each
    attribute on the input vector X to [0,1] or [-1,+1], or standardize
    it to have mean 0 and variance 1. Note that the *same* scaling
    must be applied to the test vector to obtain meaningful
    results. See `The CookBook
    <https://sourceforge.net/apps/trac/scikit-learn/wiki/CookBook>`_
    for some examples on scaling. If your attributes have an intrinsic
    scale (e.g. word frequencies, indicator features) scaling is 
    not needed. 

  * Finding a reasonable regularization term :math:`\alpha` is 
    best done using grid search. 

  * Empirically, I found that SGD converges after observing 
    approx. 10^6 training examples. Thus, a reasonable first guess 
    for the number of iterations is `n_iter = 10**6 / n`.

.. _sgd_mathematical_formulation:

Mathematical formulation
========================

Given a set of training examples :math:`(x_1, y_1), \ldots, (x_n, y_n)` where 
:math:`x_i \in \mathbf{R}^n` and :math:`y_i \in \{-1,1\}`, our goal is to 
learn a linear scoring function :math:`f(x) = w^T x + b` with model parameters
:math:`w \in \mathbf{R}^n` and intercept :math:`b \in \mathbf{R}`. In order
to make predictions, we simply look at the sign of :math:`f(x)`.
A common choice to find the model parameters is by minimizing the regularized 
training error given by

.. math::

    E(w,b) = \sum_{i=1, l} L(y_i, f(x_i)) + \alpha R(w)

where :math:`L` is a loss function that measures model (mis)fit and :math:`R` is a
regularization term (aka penalty) that penalizes model complexity; :math:`\alpha > 0`
is a non-negative hyperparameter. 




.. figure:: ../auto_examples/sgd/images/plot_loss_functions.png
   :align: center
   :scale: 50

.. figure:: ../auto_examples/sgd/images/plot_penalties.png
   :align: center
   :scale: 50

SGD
---

Stochastic gradient descent is an optimization method for unconstrained 
optimization problems. In contrast to (batch) gradient descent, SGD
approximates the true gradient of :math:`E(w,b)` by considering a 
single training example at a time. 

The class SGD implements a simple first-order SGD learning routine, which 
considers the first derivative of the loss function. The algorithm iterates
over the training examples and for each example updates the model parameters 
according to the update rule given by

.. math::

    w <- w + ...




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

The implementation of SGD is based on the `Stochastic Gradient SVM 
<http://leon.bottou.org/projects/sgd>`_  of Léon Bottou. Similar to SvmSGD, 
the weight vector is represented as the product of a scalar and a vector 
which allows an efficient weight update in the case of L2 regularization. 
The code is written in Cython.

.. topic:: References:

