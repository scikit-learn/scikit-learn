=======================
Support Vector Machines
=======================

**Support vector machines (SVMs)** are a set of supervised learning
methods used for classification, regression and outlayer detection.

The advantages of Support Vector Machines are:

    - Effective high dimensional spaces.

    - Still effective in cases where number of dimensions is greater
      than the number of samples.

    - Uses a subset of training points in the decission function (called
      support vectors), so it is also memory efficient.

    - Versatile: different :ref:`svm_kernels <kernels>` can be
      specified for the decission function. Common kernels are
      provided, but it is also possibly to specify custom kernels.

The dissadvantes of Support Vector Machines include:

    - If the number of features is much greater than the number of
      samples, the method is likely to give poor performance.

    - SVMs do not directly provide probability estimates, so these
      must be calculated using indirect techniques. In our case, these
      techniques imply conducting five-fold cross-validation, so
      performance can suffer.  See method predict_proba for more
      information.


Classification
==============

Suppose some given data points each belong to one of N classes, and
the goal is to decide which class a new data point will be in. The
classes that permform this task are SVC, NuSVC and LinearSVC.

SVC and NuSVC are similar methods, but accept slightly different set
of parameters and have different mathematical formulations (see
section :ref:`svm_mathematical_formulation`). On the other hand,
LinearSVC is another implemntation of SVC optimized in the case of a
linear kernel. Note that other classes this one does not accept
keyword 'kernel', as this is assumed to be linear. It also lacks some
of the memebrs of SVC and NuSVC, like support\_.


.. figure:: ../auto_examples/svm/images/plot_iris.png
   :align: center
   :scale: 50


**Complete class reference:**

.. autoclass:: scikits.learn.svm.SVC
   :members:
   :inherited-members:

.. autoclass:: scikits.learn.svm.NuSVC
   :members:
   :inherited-members:


The following class 
.. autoclass:: scikits.learn.svm.LinearSVC
   :members:
   :inherited-members:


Using Custom Kernels
--------------------

.. TODO: this is not restricted to classification

You can also use your own defined kernels by passing a function to the
keyword `kernel` in the constructor.

Your kernel must take as arguments two matrices and return a third matrix.

The following code defines a linear kernel and creates a classifier
instance that will use that kernel::

    >>> import numpy as np
    >>> from scikits.learn import svm
    >>> def my_kernel(x, y):
    ...     return np.dot(x, y.T)
    ... 
    >>> clf = svm.SVC(kernel=my_kernel)


Classifiers with custom kernels behave the same way as any other
classifiers, except that:

    * Support vectors do no longer represent the vectors, but rather are
      indices of the support vectors for the training vectors.

    * A reference (and not a copy) of the first argument in the fit()
      method is stored for future reference. If that array changes
      between the use of fit() and predict() you will have
      unexpected results.

.. note::

Examples
--------

For a complete example of custom kernels see
:ref:`example_svm_plot_custom_kernel.py`. 

.. TODO: precomputed kernels.

Regression
==========

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

____

.. autoclass:: scikits.learn.svm.SVR
   :members:
   :inherited-members:

.. autoclass:: scikits.learn.svm.NuSVR
   :members:
   :inherited-members:

Density estimation
=======================

One-class SVM is used for outliers detection, that is, given a set of
samples, it will detect the soft boundary of that set so as to classify
new points as belonging to that set or not.

____

.. autoclass:: scikits.learn.svm.OneClassSVM
   :members:
   :inherited-members:

.. note::

    For a complete example on one class SVM see 
    :ref:`example_svm_plot_svm_oneclass.py` example.

.. figure:: ../auto_examples/svm/images/plot_svm_oneclass.png
   :align: center
   :scale: 50


Examples
========

See :ref:`svm_examples` for a complete list of examples.


Support Vector machines for sparse data
=======================================

There is support for sparse data given in any matrix in a format
supported by scipy.sparse. See module scikits.learn.sparse.svm.


Tips on Practical Use
=====================

  * Support Vector Machine algorithms are not scale-invariant, so it
  is highly recommended to scale your data. For example, scale each
  attribute on the input vector X to [0,1] or [-1,+1], or standarize
  it to have mean 0 and variance 1. Note that the *same* scaling must
  be applied to the test vector to obtain meaningful results. See `The
  CookBook
  <https://sourceforge.net/apps/trac/scikit-learn/wiki/CookBook>`_ for
  some examples on scaling.

  * nu in nu-SVC/one-class-SVM/nu-SVR approximates the fraction of
  training errors and support vectors.

  * If data for classification are unbalanced (e.g. many positive and
    few negative), try different penalty parameters C.

  * Specify larger cache size (keyworkd cache) for huge problems.


.. _svm_mathematical_formulation:

Mathematical formulation
========================


A support vector machine constructs a hyperplane or set of hyperplanes
in a high or infinite dimensional space, which can be used for
classification, regression or other tasks. Intuitively, a good
separation is achieved by the hyperplane that has the largest distance
to the nearest training datapoints of any class (so-called functional
margin), since in general the larger the margin the lower the
generalization error of the classifier.


.. figure:: ../auto_examples/svm/images/plot_separating_hyperplane.png
   :align: center
   :scale: 50


In SVC The decision function in this case will be:

.. math:: sgn(\sum_{i=1}^l \alpha_i K(x_i, x) + \rho)
where :math:`\alpha, \rho` can be accessed through fields `support_` and
`intercept_` of the classifier instance, respectevely.

    - *kernel function*.
      Can be any of the following: linear :math:`(<x_i, x_j'>)`,
      polynomial :math:`(\gamma <x, x'> + r)^d, \gamma > 0`
      radial basis :math:`exp(-\gamma |x-x'|^2), \gamma > 0`
      sigmoid :math:`tanh(<x_i, x_j> + r)`,
      where :math:`\gamma, r`, and :math:`d` are kernel parameters.

    - *penalty*. C > 0 is the penalty parameter of the error term.


Implementation details
======================
Internally, we use libsvm[1] to handle all computations. Libsvm is wrapped
using C and Cython.

.. [1] http://www.csie.ntu.edu.tw/~cjlin/libsvm/


References
==========
For a description of the implementation and details of the algorithms
used, please refer to

    - http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf
    - http://en.wikipedia.org/wiki/Support_vector_machine
