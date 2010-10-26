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

    - Efficiency.

    - Ease of implementation (lots of opportunities for code tuning). 

The disadvantages of Stochastic Gradient Descent include:

    - SGD requires a number of hyperparameters including the regularization
      parameter and the number of iterations. 

    - SGD is sensitive to feature scaling. 

Classification
==============

.. topic:: Examples:

 * :ref:`example_sgd_plot_separating_hyperplane.py`,

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
For PEGASOS, training time indeed decreases as the size of the training set 
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
    best done using grid search `for alpha in 10.0**-np.arange(1,7)`.

  * Empirically, we found that SGD converges after observing 
    approx. 10^6 training examples. Thus, a reasonable first guess 
    for the number of iterations is `n_iter = np.ceil(10**6 / n)`, 
    where `n` is the size of the training set.

.. topic:: References:

 * `"Efficient BackProp" <yann.lecun.com/exdb/publis/pdf/lecun-98b.pdf>`_ 
   Y. LeCun, L. Bottou, G. Orr, K. Müller - In Neural Networks: Tricks of the Trade 1998.

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

    E(w,b) = \sum_{i=1}^{n} L(y_i, f(x_i)) + \alpha R(w)

where :math:`L` is a loss function that measures model (mis)fit and :math:`R` is a
regularization term (aka penalty) that penalizes model complexity; :math:`\alpha > 0`
is a non-negative hyperparameter. 

Different choices for :math:`L` entail different classifiers such as 

   - Hinge: (soft-margin) Support Vector Machines.
   - Log:   Logistic Regression.
   - Least-Squares: Ridge Regression. 

All of the above loss functions can be regarded as an upper bound on the 
misclassification error (0-1 loss) as shown on the Figure below. 

.. figure:: ../auto_examples/sgd/images/plot_loss_functions.png
   :align: center
   :scale: 50

Popular choices for :math:`R` include:

   - L2 norm: :math:`R(w) := \frac{1}{2} \sum_{i=1}^{n} w_i^2`, 
   - L1 norm: :math:`R(w) := \sum_{i=1}^{n} |w_i|`, which leadsin sparse solutions.
   - Elastic Net: :math:`R(w) := \rho \frac{1}{2} \sum_{i=1}^{n} w_i^2 + (1-\rho) \sum_{i=1}^{n} |w_i|`, a convex combination of L2 and L1. 

The Figure below shows the contours of the different regularization terms 
in the parameter space when :math:`R(w) = 1`. 

.. figure:: ../auto_examples/sgd/images/plot_penalties.png
   :align: center
   :scale: 50

SGD
---

Stochastic gradient descent is an optimization method for unconstrained 
optimization problems. In contrast to (batch) gradient descent, SGD
approximates the true gradient of :math:`E(w,b)` by considering a 
single training example at a time. 

The class SGD implements a first-order SGD learning routine. The algorithm 
iterates over the training examples and for each example updates the model 
parameters according to the update rule given by

.. math::

    w \leftarrow w - \eta (\alpha \frac{\partial R(w)}{\partial w} 
    + \frac{\partial L(w^T x_i + b, y_i)}{\partial w})

where :math:`\eta` is the learning rate which controls the step-size 
in the parameter space. 
The intercept :math:`b` is updated similarly but without regularization.

.. TODO multiclass case ?/

The model parameters can be accessed through the members coef\_ and intercept\_:

     - Member coef\_ holds the weights :math:`w`

     - Member intercept\_ holds :math:`b`

.. topic:: References:

 * `"Solving large scale linear prediction problems using stochastic gradient descent algorithms"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.58.7377>`_ 
   T. Zhang - In Proceedings of ICML '04.

 * `"Pegasos: Primal estimated sub-gradient solver for svm" 
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.74.8513>`_
   S. Shalev-Shwartz, Y. Singer, N. Srebro - In Proceedings of ICML '07.
   
 * `"Regularization and variable selection via the elastic net"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.124.4696>`_ 
   H. Zou, T. Hastie - Journal of the Royal Statistical Society Series B, 67 (2), 301-320.


Frequently Asked Questions
==========================




Implementation details
======================

The implementation of SGD is influenced by the `Stochastic Gradient SVM 
<http://leon.bottou.org/projects/sgd>`_  of Léon Bottou. Similar to SvmSGD, 
the weight vector is represented as the product of a scalar and a vector 
which allows an efficient weight update in the case of L2 regularization
and the intercept is updated with a smaller learning rate (multiplied by 0.01) 
to account for the fact that it is updated more frequently. 

For L1 regularization (and the Elastic Net) we use the truncated gradient
algorithm proposed by Tsuruoka et al. 2009. 

The code is written in Cython.

.. topic:: References:

 * `"Stochastic Gradient Descent" <leon.bottou.org/projects/sgd>`_
    L. Bottou - Website, 2010

 * `"Stochastic gradient descent training for l1-regularized log-linear models with cumulative penalty"
   <www.aclweb.org/anthology/P/P09/P09-1054.pdf>`_
   Y. Tsuruoka, J. Tsujii, S. Ananiadou -  In Proceedings of the AFNLP/ACL '09.
