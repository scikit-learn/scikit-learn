=======================
Support Vector Machines
=======================

.. currentmodule:: scikits.learn.svm

**Support vector machines (SVMs)** are a set of supervised learning
methods used for :ref:`classification <svm_classification>`,
:ref:`regression <svm_regression>` and :ref:`outliers detection
<svm_outlier_detection>`.

The advantages of Support Vector Machines are:

    - Effective in high dimensional spaces.

    - Still effective in cases where number of dimensions is greater
      than the number of samples.

    - Uses a subset of training points in the decision function (called
      support vectors), so it is also memory efficient.

    - Versatile: different :ref:`svm_kernels` can be
      specified for the decision function. Common kernels are
      provided, but it is also possible to specify custom kernels.

The disadvantages of Support Vector Machines include:

    - If the number of features is much greater than the number of
      samples, the method is likely to give poor performances.

    - SVMs do not directly provide probability estimates, so these
      must be calculated using indirect techniques. In our case, these
      techniques imply conducting five-fold cross-validation, so
      performance can suffer.  See method predict_proba for more
      information.


.. _svm_classification:

Classification
==============

Suppose some given data points each belonging to one of N classes, and
the goal is to decide which class a new data point will be in. The
classes that perform this task are :class:`SVC`, :class:`NuSVC` and
:class:`LinearSVC`.

:class:`SVC` and :class:`NuSVC` are similar methods, but accept
slightly different sets of parameters and have different mathematical
formulations (see section :ref:`svm_mathematical_formulation`). On the
other hand, :class:`LinearSVC` is another implementation of SVC
optimized in the case of a linear kernel. Note that :class:`LinearSVC`
does not accept keyword 'kernel', as this is assumed to be linear. It
also lacks some of the members of SVC and NuSVC, like support\_.


.. figure:: ../auto_examples/svm/images/plot_iris.png
   :target: ../auto_examples/svm/plot_iris.html
   :align: center


As other classifiers, SVC and NuSVC have to be fitted with two arrays:
an array X of size [m_samples, n_features] holding the training
samples, and an array Y of size [n_samples] holding the target values
(class labels) for the training samples::


    >>> from scikits.learn import svm
    >>> X = [[0., 0.], [1., 1.]]
    >>> Y = [0, 1]
    >>> clf = svm.SVC()
    >>> clf.fit (X, Y)
    SVC(kernel='rbf', C=1.0, probability=False, degree=3, coef0=0.0, eps=0.001,
      cache_size=100.0, shrinking=True, gamma=0.5)

After being fitted, the model can then be used to predict new values::

    >>> clf.predict ([[2., 2.]])
    array([ 1.])

SVMs perform classification as a function of some subset of the
training data, called the support vectors. These vectors can be
accessed in member `support_`:

    >>> clf.support_
    array([[ 0.,  0.],
           [ 1.,  1.]])

Member `n_support_` holds the number of support vectors for each class:

    >>> clf.n_support_
    array([1, 1], dtype=int32)


.. topic:: Examples:

 * :ref:`example_svm_plot_iris.py`,
 * :ref:`example_svm_plot_separating_hyperplane.py`,
 * :ref:`example_svm_plot_svm_anova.py`,
 * :ref:`example_svm_plot_svm_nonlinear.py`

.. _svm_regression:

Regression
==========

The method of Support Vector Classification can be extended to solve
regression problems. This method is called Support Vector Regression.

The model produced by support vector classification (as described
above) depends only on a subset of the training data, because the cost
function for building the model does not care about training points
that lie beyond the margin. Analogously, the model produced by Support
Vector Regression depends only on a subset of the training data,
because the cost function for building the model ignores any training
data close to the model prediction.

There are two flavors of Support Vector Regression: :class:`SVR` and
:class:`NuSVR`.

As with classification classes, the fit method will take as
argument vectors X, y, only that in this case y is expected to have
floating point values instead of integer values.


.. topic:: Examples:

 * :ref:`example_svm_plot_svm_regression.py`

.. _svm_outlier_detection:

Density estimation, outliers detection
=======================================

One-class SVM is used for outliers detection, that is, given a set of
samples, it will detect the soft boundary of that set so as to
classify new points as belonging to that set or not. The class that
implements this is called :class:`OneClassSVM`


In this case, as it is a type of unsupervised learning, the fit method
will only take as input an array X, as there are no class labels.


.. figure:: ../auto_examples/svm/images/plot_oneclass.png
   :target: ../auto_examples/svm/plot_oneclass.html
   :align: center
   :scale: 50


.. topic:: Examples:

 * :ref:`example_svm_plot_oneclass.py`


.. currentmodule:: scikits.learn.svm.sparse


Support Vector machines for sparse data
=======================================

There is support for sparse data given in any matrix in a format
supported by scipy.sparse. Classes have the same name, just prefixed
by the `sparse` namespace, and take the same arguments, with the
exception of training and test data, which is expected to be in a
matrix format defined in scipy.sparse.

For maximum efficiency, use the CSR matrix format as defined in
`scipy.sparse.csr_matrix
<http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html>`_.

Implemented classes are :class:`SVC`, :class:`NuSVC`,
:class:`SVR`, :class:`NuSVR`, :class:`OneClassSVM`,
:class:`LinearSVC`.


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

  * Support Vector Machine algorithms are not scale invariant, so it
    is highly recommended to scale your data. For example, scale each
    attribute on the input vector X to [0,1] or [-1,+1], or standardize
    it to have mean 0 and variance 1. Note that the *same* scaling
    must be applied to the test vector to obtain meaningful
    results. See `The CookBook
    <https://sourceforge.net/apps/trac/scikit-learn/wiki/CookBook>`_
    for some examples on scaling.

  * Parameter nu in NuSVC/OneClassSVM/NuSVR approximates the fraction
    of training errors and support vectors.

  * If data for classification are unbalanced (e.g. many positive and
    few negative), try different penalty parameters C.

  * Specify larger cache size (keyword cache) for huge problems.


.. _svm_kernels:

Kernel functions
================

The *kernel function* can be any of the following:

  * linear: :math:`<x_i, x_j'>`.

  * polynomial: :math:`(\gamma <x, x'> + r)^d`. d is specified by
    keyword `degree`.

  * rbf (:math:`exp(-\gamma |x-x'|^2), \gamma > 0`). :math:`\gamma` is
    specified by keyword gamma.

  * sigmoid (:math:`tanh(<x_i,x_j> + r)`).

Different kernels are specified by keyword kernel at initialization::

    >>> linear_svc = svm.SVC(kernel='linear')
    >>> linear_svc.kernel
    'linear'
    >>> rbf_svc = svm.SVC (kernel='rbf')
    >>> rbf_svc.kernel
    'rbf'


Custom Kernels
--------------

You can define your own kernels by either giving the kernel as a
python function or by precomputing the Gram matrix.

Classifiers with custom kernels behave the same way as any other
classifiers, except that:

    * Support vectors do no longer represent the vectors, but rather are
      indices of the support vectors for the training vectors.

    * A reference (and not a copy) of the first argument in the fit()
      method is stored for future reference. If that array changes
      between the use of fit() and predict() you will have
      unexpected results.


Using python functions as kernels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Passing the gram matrix
~~~~~~~~~~~~~~~~~~~~~~~

Set kernel='precomputed' and pass the gram matrix instead of X in the
fit method.


.. topic:: Examples:

 * :ref:`example_svm_plot_custom_kernel.py`.


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

 * `"Automatic Capacity Tuning of Very Large VC-dimension Classifiers"
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.17.7215>`_
   I Guyon, B Boser, V Vapnik - Advances in neural information
   processing 1993,


 * `"Support-vector networks"
   <http://www.springerlink.com/content/k238jx04hm87j80g/>`_
   C. Cortes, V. Vapnik, Machine Leaming, 20, 273-297 (1995)



NuSVC
-----

We introduce a new parameter :math:`\nu` which controls the number of
support vectors and training errors. The parameter :math:`\nu \in (0,
1]` is an upper bound on the fraction of training errors and a lower
bound of the fraction of support vectors.


Frequently Asked Questions
==========================

     * Q: Can I get the indices of the support vectors instead of the
       support vectors ?

       A: The underlying C implementation does not provide this
       information.


Implementation details
======================

Internally, we use `libsvm`_ and `liblinear`_ to handle all
computations. These libraries are wrapped using C and Cython.

.. _`libsvm`: http://www.csie.ntu.edu.tw/~cjlin/libsvm/
.. _`liblinear`: http://www.csie.ntu.edu.tw/~cjlin/liblinear/

.. topic:: References:

  For a description of the implementation and details of the algorithms
  used, please refer to

    - `LIBSVM: a library for Support Vector Machines
      <http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf>`_

    - `LIBLINEAR -- A Library for Large Linear Classification
      <http://www.csie.ntu.edu.tw/~cjlin/liblinear/>`_


