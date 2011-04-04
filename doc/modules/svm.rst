
.. _svm:

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

    - SVMs do not directly provide probability estimates, these are
      calculated using five-fold cross-validation, and thus
      performance can suffer.

.. TODO: add reference to probability estimates

.. _svm_classification:

Classification
==============

:class:`SVC`, :class:`NuSVC` and :class:`LinearSVC` are classes
capable of performing multi-class classification on a dataset.


.. figure:: ../auto_examples/svm/images/plot_iris_1.png
   :target: ../auto_examples/svm/plot_iris.html
   :align: center


:class:`SVC` and :class:`NuSVC` are similar methods, but accept
slightly different sets of parameters and have different mathematical
formulations (see section :ref:`svm_mathematical_formulation`). On the
other hand, :class:`LinearSVC` is another implementation of Support
Vector Classification for the case of a linear kernel. Note that
:class:`LinearSVC` does not accept keyword 'kernel', as this is
assumed to be linear. It also lacks some of the members of
:class:`SVC` and :class:`NuSVC`, like support\_.

As other classifiers, :class:`SVC`, :class:`NuSVC` and
:class:`LinearSVC` take as input two arrays: an array X of size
[n_samples, n_features] holding the training samples, and an array Y
of integer values, size [n_samples], holding the class labels for the
training samples::


    >>> from scikits.learn import svm
    >>> X = [[0, 0], [1, 1]]
    >>> Y = [0, 1]
    >>> clf = svm.SVC()
    >>> clf.fit(X, Y)
    SVC(kernel='rbf', C=1.0, probability=False, degree=3, coef0=0.0, tol=0.001,
      cache_size=100.0, shrinking=True, gamma=0.5)

After being fitted, the model can then be used to predict new values::

    >>> clf.predict([[2., 2.]])
    array([ 1.])

SVMs decision function depends on some subset of the training data,
called the support vectors. Some properties of these support vectors
can be found in members `support_vectors_`, `support_` and
`n_support`::

    >>> # get support vectors
    >>> clf.support_vectors_
    array([[ 0.,  0.],
           [ 1.,  1.]])
    >>> # get indices of support vectors
    >>> clf.support_ # doctest: +ELLIPSIS
    array([0, 1]...)
    >>> # get number of support vectors for each class
    >>> clf.n_support_ # doctest: +ELLIPSIS
    array([1, 1]...)


Multi-class classification
--------------------------

:class:`SVC` and :class:`NuSVC` implement the "one-against-one"
approach (Knerr et al., 1990) for multi- class classification. If
n_class is the number of classes, then n_class * (n_class - 1)/2
classifiers are constructed and each one trains data from two classes.


    >>> X = [[0], [1], [2], [3]]
    >>> Y = [0, 1, 2, 3]
    >>> clf = svm.SVC()
    >>> clf.fit(X, Y)
    SVC(kernel='rbf', C=1.0, probability=False, degree=3, coef0=0.0, tol=0.001,
      cache_size=100.0, shrinking=True, gamma=0.25)
    >>> dec = clf.decision_function([[1]])
    >>> dec.shape[1] # 4 classes: 4*3/2 = 6
    6


On the other hand, :class:`LinearSVC` implements "one-vs-the-rest"
multi-class strategy, thus training n_class models. If there are only
two classes, only one model is trained.


    >>> lin_clf = svm.LinearSVC()
    >>> lin_clf.fit(X, Y)
    LinearSVC(loss='l2', C=1.0, dual=True, fit_intercept=True, penalty='l2',
         multi_class=False, tol=0.0001, intercept_scaling=1)
    >>> dec = lin_clf.decision_function([[1]])
    >>> dec.shape[1]
    4


See :ref:`svm_mathematical_formulation` for a complete description of
the decision function.


Unbalanced problems
--------------------

In problems where it is desired to give more importance to certain
classes or certain individual samples keywords ``class_weight`` and
``sample_weight`` can be used.

:class:`SVC` (but not :class:`NuSVC`) implement a keyword
``class_weight`` in the fit method. It's a dictionary of the form
``{class_label : value}``, where value is a floating point number > 0
that sets the parameter C of class ``class_label`` to C * value.

.. figure:: ../auto_examples/svm/images/plot_separating_hyperplane_unbalanced_1.png
   :target: ../auto_examples/svm/plot_separating_hyperplane_unbalanced.html
   :align: center
   :scale: 75


:class:`SVC`, :class:`NuSVC`, :class:`SVR`, :class:`NuSVR` and
:class:`OneClassSVM` implement also weights for individual samples in method
``fit`` through keyword sample_weight.


.. figure:: ../auto_examples/svm/images/plot_weighted_samples_1.png
   :target: ../auto_examples/svm/plot_weighted_samples.html
   :align: center
   :scale: 75


.. topic:: Examples:

 * :ref:`example_svm_plot_iris.py`,
 * :ref:`example_svm_plot_separating_hyperplane.py`,
 * :ref:`example_svm_plot_separating_hyperplane_unbalanced.py`
 * :ref:`example_svm_plot_svm_anova.py`,
 * :ref:`example_svm_plot_svm_nonlinear.py`
 * :ref:`example_svm_plot_weighted_samples.py`,


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

    >>> from scikits.learn import svm
    >>> X = [[0, 0], [2, 2]]
    >>> y = [0.5, 2.5]
    >>> clf = svm.SVR()
    >>> clf.fit(X, y)
    SVR(kernel='rbf', C=1.0, probability=False, degree=3, shrinking=True, p=0.1,
      tol=0.001, cache_size=100.0, coef0=0.0, nu=0.5, gamma=0.5)
    >>> clf.predict([[1, 1]])
    array([ 1.5])


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


.. figure:: ../auto_examples/svm/images/plot_oneclass_1.png
   :target: ../auto_examples/svm/plot_oneclass.html
   :align: center
   :scale: 75


.. topic:: Examples:

 * :ref:`example_svm_plot_oneclass.py`
 * :ref:`example_applications_plot_species_distribution_modeling.py`

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
:math:`O(n_{features} \times n_{samples}^2)` and
:math:`O(n_{features} \times n_{samples}^3)` depending on how efficiently
the `libsvm`_ cache is used in practice (dataset dependent). If the data
is very sparse :math:`n_{features}` should be replaced by the average number
of non-zero features in a sample vector.

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

  * In SVC, if data for classification are unbalanced (e.g. many
    positive and few negative), set class_weight='auto' and/or try
    different penalty parameters C.

  * Specify larger cache size (keyword cache) for huge problems.

  * The underlying :class:`LinearSVC` implementation uses a random
    number generator to select features when fitting the model. It is
    thus not uncommon, to have slightly different results for the same
    input data. If that happens, try with a smaller tol parameter.


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
    >>> rbf_svc = svm.SVC(kernel='rbf')
    >>> rbf_svc.kernel
    'rbf'


Custom Kernels
--------------

You can define your own kernels by either giving the kernel as a
python function or by precomputing the Gram matrix.

Classifiers with custom kernels behave the same way as any other
classifiers, except that:

    * Field `support_vectors\_` is now empty, only indices of support
      vectors are stored in `support_`

    * A reference (and not a copy) of the first argument in the fit()
      method is stored for future reference. If that array changes
      between the use of fit() and predict() you will have unexpected
      results.


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

Using the Gram matrix
~~~~~~~~~~~~~~~~~~~~~

Set kernel='precomputed' and pass the Gram matrix instead of X in the
fit method.

.. TODO: inline example

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


.. figure:: ../auto_examples/svm/images/plot_separating_hyperplane_1.png
   :align: center
   :scale: 75

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

This parameters can be accessed through the members `dual_coef\_`
which holds the product :math:`y_i \alpha_i`, `support_vectors\_` which
holds the support vectors, and `intercept\_` which holds the independent
term :math:`-\rho` :

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


