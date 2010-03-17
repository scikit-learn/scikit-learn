=======================
Support Vector Machines
=======================

**Support vector machines (SVMs)** are a set of supervised learning
methods used for classification and regression. In simple words, given
a set of training examples, witheach sample marked as belonging to one
of the multiple categories, an SVM training algorithm builds a model
that predicts whether a new example falls into one category or the
other.

More formally, a support vector machine constructs a hyperplane or set
of hyperplanes in a high or infinite dimensional space, which can be
used for classification, regression or other tasks. Intuitively, a
good separation is achieved by the hyperplane that has the largest
distance to the nearest training datapoints of any class (so-called
functional margin), since in general the larger the margin the lower
the generalization error of the classifier.

Suppose some given data points each belong to one of two classes, and
the goal is to decide which class a new data point will be in. This
classification will be performed by creating a hyperplane that
maximizes the distance between any two classes.

.. literalinclude:: ../../examples/plot_svm_hyperplane.py
.. image:: svm_data/separating_hyperplane_2D.png

The original optimal hyperplane algorithm was a linear
classifier. However, in 1992, Bernhard Boser, Isabelle Guyon and
Vapnik suggested a way to replace every dot product by a non-linear
kernel function. This allows the algorithm to fit the maximum-margin
hyperplane in a transformed feature space. The transformation may be
non-linear and the transformed space high dimensional; thus though the
classifier is a hyperplane in the high-dimensional feature space, it
may be non-linear in the original input space.

If the kernel used is a Gaussian radial basis function, the
corresponding feature space is a Hilbert space of infinite
dimension. Maximum margin classifiers are well regularized, so the
infinite dimension does not spoil the results. Available kernels are,

  * linear :math:`(1 + <x, x'>)`
  * polynomial :math:`(1 + <x, x'>)^d`
  * radial :math:`exp(-\gamma |x-x'|^2)`
  * sigmoid :math:`tanh(x_i x_j + c)`


The exclusive-OR is the simplest problem that cannot be solved using a
linear kernel. In this problem, point (x, y) belongs has target 1 if
and only if x > 0 XOR y > 0. In the following example, we create a
training set of random points X with target Y = XOR(X). We see that
the SVM correctly draws the decision function.

.. literalinclude:: ../../examples/plot_svm_nonlinear.py
.. image:: svm_data/separating_nonlinear.png

Complete class reference:

.. autoclass:: scikits.learn.svm.SVC
   :members:

Regression
----------
The method of Support Vector Classification can be extended to solve
the regression problem. This method is called Support Vector
Regression.

The model produced by support vector classification (as described
above) depends only on a subset of the training data, because the cost
function for building the model does not care about training points
that lie beyond the margin. Analogously, the model produced by Support
Vector Regression depends only on a subset of the training data,
because the cost function for building the model ignores any training
data close to the model prediction.


.. autoclass:: scikits.learn.svm.SVR
   :members:


Distribution estimation
=======================
One-class SVM was proposed by Scholkopf et al. (2001) for estimating
the support of a high-dimensional distribution. Given training vectors
:math:`x_i \in \mathbb{R}^n, i=1, .., l` without any class
information, the primal form is:

.. math::    \min_{w, b, \xi} {1 \over 2} w^T w - \rho + {1 \over \nu l} \sum_{i=1}^l \xi_i

             \textrm{subject to} w^T \phi(x_i) \geq \rho - \xi_i

             \xi_i \geq 0, i=1,...,l

Scaling
=======

TODO

.. Mathematical formulation (Model selection)
.. ========================


.. C-support vector classification (C-SVC)
.. ---------------------------------------
.. Given training vectors :math:`x_i \in \mathbb{R}^n , i=1, ..., l` in
.. two classes, and a vector :math:`y \in \mathbb{R}^l` such that
.. :math:`y_i \in {1, -1}`, C-SVC solves the following primal problem:

.. .. math:: \min_{w, b, \xi} {1 \over 2} w^T w + C \sum_{i=1}^l \xi_i
.. .. math:: \textrm{subject to}\ y_i (w^T \phi(x_i) + b) \geq 1 - \xi_i
.. .. math:: \xi_i >= 0, i=1, .., l

.. Here training vectors :math:`x_i` are mapped into a higher (maybe
.. infinite) dimensional space by the function :math:`phi`. The decision
.. function is

.. .. math::    sgn(\sum_{i=0}^l y_i \alpha_i K(x_i, x) + b)

.. This is implemented in class SVC


.. Nu-Support Vector Classification
.. --------------------------------
.. The nu-Support Vector Classification uses a new parameter :math:`\nu`
.. which controls the number of support vectors and trainign errors. The
.. parameter :math:`nu \in (0, 1]` is an upper bound on the fraction of
.. training errors and a lower bound of the fraction of support vectors.

.. Given training vectors :math:`x_i \in \mathbb{R}^n , i=1, ..., l` in
.. two classes, and a vector :math:`y \in \mathbb{R}^l` such that
.. :math:`y_i \in {1, -1}`, C-SVC solves the following primal problem:

.. .. math:: \min_{w, b, \xi} {1 \over 2} w^T w - \nu \rho + {1 \over 2} \sum_{i=1}^l \xi_i

..           \textrm{subject to}\ y_i (w^T \phi(x_i) + b) \geq \rho - \xi_i

..           \xi_i \geq 0, i=1, .., l, \rho \geq 0

.. The decision function is:

.. .. math::    sgn(\sum_{i=1}^l y_i \alpha_i K(x_i, x) + b

.. This is implemented in SVC(impl='nu-svc')



Low-level implementation
========================

Internally, we use libsvm[1] to handle all computations. Libsvm is binded
through some wrappers written in C and Cython.

.. [1] http://www.csie.ntu.edu.tw/~cjlin/libsvm/


References
==========

For a description of the implementation and details of the algorithms
used, please refer to
http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf

http://en.wikipedia.org/wiki/Support_vector_machine
