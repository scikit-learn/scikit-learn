=======================
Support Vector Machines
=======================

Support vector machines (SVMs) are a set of supervised learning
methods used for classification and regression. In simple words, given
a set of training examples, witheach sample marked as belonging to one
of the multiple categories, an SVM training algorithm builds a model
that predicts whether a new example falls into one category or the
other.



Classification
==============

Classification is implemented in class SVC. There are two variants of the algorithm, C-SVC and Nu-SVC.


.. autoclass:: scikits.learn.svm.SVC
   :members:


C-support vector classification (C-SVC)
---------------------------------------
Given training vectors :math:`x_i \in \mathbb{R}^n , i=1, ..., l` in two classes, and a vector :math:`y \in \mathbb{R}^l` such that :math:`y_i \in {1, -1}`, C-SVC solves the following primal problem:

.. math::    \min_{w, b, \xi} {1 \over 2} w^T w + C \sum_{i=1}^l \xi_i

              \textrm{subject to}\ y_i (w^T \phi(x_i) + b) \geq 1 - \xi_i

              \xi_i >= 0, i=1, .., l

Here training vectors :math:`x_i` are mapped into a higher (maybe infinite) dimensional space by the function :math:`phi`. The decision function is

.. math::    sgn(\sum_{i=0}^l y_i \alpha_i K(x_i, x) + b)


Nu-Support Vector Classification
--------------------------------
The nu-Support Vector Classification uses a new parameter :math:`\nu`
which controls the number of support vectors and trainign errors. The
parameter :math:`nu \in (0, 1]` is an upper bound on the fraction of
training errors and a lower bound of the fraction of support vectors.

Given training vectors :math:`x_i \in \mathbb{R}^n , i=1, ..., l` in two classes, and a vector :math:`y \in \mathbb{R}^l` such that :math:`y_i \in {1, -1}`, C-SVC solves the following primal problem:

.. math::    \min_{w, b, \xi} {1 \over 2} w^T w - \nu \rho + {1 \over 2} \sum_{i=1}^l \xi_i

              \textrm{subject to}\ y_i (w^T \phi(x_i) + b) \geq \rho - \xi_i

              \xi_i \geq 0, i=1, .., l, \rho \geq 0

The decision function is:

.. math::    sgn(\sum_{i=1}^l y_i \alpha_i K(x_i, x) + b

Implementation
--------------

Both problems are implemented in class scikits.learn.svm.SVC . This class follows the pattern of an estimator. See section Parameters for more details about available parameters.

Examples
--------
.. literalinclude:: ../../examples/plot_svm.py



Distribution estimation
=======================
One-class

Regression
==========



epsilon-support vector regression (epsilon-SVR), and ν-support vector regression
(ν-SVR)


Parameters
==========

Kernel Type
-----------
The kernel of a Support Vector Machine is the function that computes
inner products in the transformed space.

Choices are one of 

  * linear :math:`(1 + <x, x'>)`
  * poly :math:`(1 + <x, x'>)^d`
  * rbf :math:`exp(-\gamma |x-x'|^2)`
  * 'sigmoid'
  * 'precomputed'.


Support Vectors
===============

Support vector machines ultimately depend on a subset of the training
set (the support vectors) to represent the classification boundary.

You can access the support vectors of the object in all classes using the member ``support_``

>>> from scikits.learn import svm
>>> X = [[-1, 1], [0, 0], [1, 1]]
>>> Y = [1, 2, 2]
>>> clf = svm.SVC(kernel='linear')
>>> clf.fit(X, Y) #doctest: +ELLIPSIS
<scikits.learn.svm.SVC object at ...>
>>> print clf.support_ #doctest: +NORMALIZE_WHITESPACE
          [[-1.  1.]
           [ 0.  0.]]

Coefficient for support vectors

TODO: include image

Low-level implementation
========================

Internally, we use libsvm[1] to handle all computations. Libsvm is binded
through some wrappers written in C and Cython.

.. [1] http://www.csie.ntu.edu.tw/~cjlin/libsvm/
