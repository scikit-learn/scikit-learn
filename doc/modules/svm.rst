=======================
Support Vector Machines
=======================

Support vector machines (SVMs) are a set of supervised learning
methods used for classification and regression. In simple words, given
a set of training examples, witheach sample marked as belonging to one
of the multiple categories, an SVM training algorithm builds a model
that predicts whether a new example falls into one category or the
other.


.. .. automodule:: scikits.learn.machine.svm
   :members:


Formulations
============

C-support vector classification (C-SVC), ν-support vector
classification (ν-SVC), distribution estimation (one- class SVM),
epsilon-support vector regression (epsilon-SVR), and ν-support vector regression
(ν-SVR)


Parameters
==========

Kernel Type
-----------
The kernel of a Support Vector Machine is the function that computes
inner products in the transformed space.

Choices are one of 

  * 'linear' :math:`(1 + <x, x'>)`
  * 'poly' :math:`(1 + <x, x'>)^d`
  * 'rbf' :math:`exp(-\gamma |x-x'|^2)`
  * 'sigmoid'
  * 'precomputed'.


Support Vectors
===============

Support vector machines ultimately depend on a subset of the training
set (the support vectors) to represent the classification boundary
Visualize support vectors

Support vectors can be retrived through the member support_ of an
instance of SVM. Example:

>>> SV.


Coefficient for support vectors

Examples
========

.. literalinclude:: ../../examples/plot_svm.py

TODO: include image

Low-level implementation
========================

Internally, we use libsvm[1] to handle all computations. Libsvm is binded
through some wrappers written in C and Cython.

.. [1] http://www.csie.ntu.edu.tw/~cjlin/libsvm/
