.. _mkl:

========================
Multiple Kernel Learning
========================

.. currentmodule:: sklearn.mkl

**Multiple Kernel Learning (MKL)** is a framework for learning an optimal combination of
several kernel functions, each capturing different patterns or similarities in the data.
This implementation provides several MKL algorithms that optimize the kernel weights
according to different strategies.
Once the weights are optimized, an internal :ref:`SVM <svm>` (Support Vector Machine)
is trained on the resulting combined kernel.

MKL classes can be used in two main ways:

- As an estimator, using the :ref:`internal SVM <mkl_svm>` to perform classification,
  regression or outlier detection.
- As a transformer, returning the combined kernel matrix, which can then be used with
  any kernel-based method.


.. _mkl_kernels:

Kernels Configuration
=====================

The MKLC classifier supports a flexible configuration of kernels through three
parameters:

- **kernels**: list of kernel functions or strings (or `None`). Each element can be a
  callable (function taking `X`, `Y`) or a string among the
  `sklearn.metrics.pairwise.PAIRWISE_KERNEL_FUNCTIONS` kernel names
  (e.g., `"linear"`, `"rbf"`, etc.). If `None`, the input `X` must be a list of
  precomputed kernel matrices.

- **kernels_scopes**: list of `"single"` or `"all"`, one per kernel. If `"all"`, the
  kernel is applied on the entire feature space. If `"single"`, the kernel is applied
  separately to each feature, producing `n_features` distinct kernels.

- **kernels_param_grids**: list of dictionaries defining grid parameters for each
  kernel. The grid is expanded via a Cartesian product: for example, the element
  `{"gamma": [0.1, 1.0], "degree": [2, 3]}` yields 4 kernels with all combinations
  of `gamma` and `degree`.

Usage example with precomputed kernels::

    >>> from sklearn.mkl import MKLC
    >>> from sklearn.metrics.pairwise import rbf_kernel
    >>> X = [[0, 0], [1, 1]]
    >>> K1 = rbf_kernel(X, gamma=0.1)
    >>> K2 = rbf_kernel(X, gamma=1.0)
    >>> K3 = rbf_kernel(X, gamma=10.0)
    >>> clf = MKLC()
    >>> clf.fit([K1, K2, K3], y=[0, 1])
    MKLC()

Usage example with kernel functions::

    >>> from sklearn.mkl import MKLC
    >>> from sklearn.metrics.pairwise import polynomial_kernel
    >>> X = [[0, 0], [1, 1]]
    >>> clf = MKLC(
    ...     kernels=["rbf", polynomial_kernel],  # Same as ["rbf", "poly"]
    ...     kernels_scopes=["all", "single"],
    ...     kernels_param_grids=[
    ...         {"gamma": [0.1, 1.0]},
    ...         {"gamma": [0.01, 0.1], "degree": [2, 3, 4]}
    ...     ],
    ... )
    >>> clf.fit(X, y=[0, 1])
    MKLC(kernels=['rbf', 'poly'],
         kernels_param_grids=[{'gamma': [0.1, 1.0]},
                              {'degree': [2, 3, 4], 'gamma': [0.01, 0.1]}],
         kernels_scopes=['all', 'single'])

In this example:

- The RBF kernel is applied on the full feature set with two gamma values → :math:`2`
  kernels.
- The polynomial kernel is applied per feature (:math:`n_{features} = 2`), with 2
  gammas and 3 degrees → :math:`6 \times n_{features} = 12` kernels.

Total number of kernels is :math:`2 + 12 = 14`.


.. _mkl_usage:

Using a Fitted MKL Instance
===========================

Once the MKL model is fitted, there are different ways to use the learned kernel
combination. Below are the different usage methods:

1. **Using MKL as a Transformer**

   After fitting the model, MKL can be used as a transformer to apply the learned
   kernel combination to new data. This is useful if you want to extract the learned
   kernel weights and apply them to a new dataset.

   Example:
    >>> from sklearn.mkl import MKLC
    >>> mkl = MKLC(kernels=["linear", "rbf"])
    >>> mkl.fit(X=[[0, 1], [1, 2]], y=[0, 1])
    >>> mkl.transform([[2, 1], [2, 3]])
    array([[1., 4.],
           [3., 8.]])

2. **Using MKL as a Callable Kernel**

   MKL can also be used in a kernel-based estimator by passing the fitted
   instance as the kernel function parameter. This allows the integration of MKL
   into the kernel-based estimator training process.

   Example:
    >>> from sklearn.mkl import MKLR
    >>> from sklearn.svm import SVR
    >>> X = [[-1, -1], [-2, -1], [1, 1], [2, 1]]
    >>> y = [-2.0, -3.0, 3.0, 4.0]
    >>> mkl = MKLR(kernels=["linear", "rbf"])
    >>> mkl.fit(X, y)
    >>> reg = SVR(kernel=mkl)  # The MKL instance is used as a kernel
    >>> reg.fit(X, y)
    >>> reg.predict(X)
    array([-1.9, -3.1,  2.9,  4.1])

3. **Using MKL as an Estimator**

   The MKL model itself can be used as an estimator by directly using the fitted model
   for prediction, scoring, etc. In this case, MKL uses its
   :ref:`internal SVM <mkl_svm>` with the optimal combination of kernels learned
   during the fitting process.

   Example:
    >>> from sklearn.mkl import OneClassMKL
    >>> mkl = OneClassMKL(kernels=["linear", "rbf"])
    >>> X = [[0], [0.44], [0.45], [0.46], [1]]
    >>> mkl.fit(X)
    >>> mkl.predict(X)
    array([-1,  1,  1,  1, -1])

4. **Using MKL in a Pipeline**

   MKL can also be integrated into machine learning pipelines to automate the process
   of fitting and transforming data. Below are examples of using MKL in a pipeline both
   as a transformer and as an estimator.

   - **MKL as a Transformer in a Pipeline**

     In this case, MKL is used as a part of a data preprocessing pipeline, where
     it transforms the data before being passed to a classifier.

     Example:
      >>> from sklearn.mkl import MKLC
      >>> from sklearn.pipeline import Pipeline
      >>> from sklearn.svm import SVC
      >>> mkl = MKLC(kernels=["linear", "rbf"])
      >>> clf = SVC()
      >>> pipeline = Pipeline([('mkl', mkl), ('svc', clf)])
      >>> pipeline.fit(X=[[-1, -1], [-2, -1], [1, 1], [2, 1]], y=[1, 1, 2, 2])
      >>> pipeline.predict([[-0.8, -1], [0.8, 1.2]])
      array([1, 2])

     Here, MKL first transforms the data using the fitted kernel combination, and
     then the transformed data is passed to the `SVC` classifier.

   - **MKL as an Estimator in a Pipeline**

     MKL can also be used as the final step in a pipeline, where it directly makes
     predictions.

     Example:
      >>> from sklearn.mkl import MKLR
      >>> from sklearn.pipeline import Pipeline
      >>> from sklearn.preprocessing import StandardScaler
      >>> X = [[-1, -1], [-2, -1], [1, 1], [2, 1]]
      >>> y = [-2.0, -3.0, 3.0, 4.0]
      >>> mkl = MKLR(kernels=["linear", "rbf"])
      >>> pipeline = Pipeline([('sc', StandardScaler()), ('mkl', mkl)])
      >>> pipeline.fit(X, y)
      >>> pipeline.predict(X)
      array([-1.90000006, -2.90017257,  2.90000006,  3.90017257])

     In this example, MKL is used as the final estimator, making predictions with its
     internal SVM.


.. _mkl_svm:

Internal SVM
============

In the MKL framework, the underlying SVM differs based on the model type: for
:class:`MKLC`, it uses a :class:`sklearn.svm.SVC` (Support Vector Classifier); for
:class:`MKLR`, it's a :class:`sklearn.svm.SVR` (Support Vector Regressor); and for
:class:`OneClassMKL`, a :class:`sklearn.svm.OneClassSVM` (One-Class Support Vector
Machine).

Setting SVM Parameters
----------------------

You can customize the internal SVM used by the MKL estimator by providing
additional parameters through the `svm_params` argument. This allows you to modify
settings for the support vector machine that are not directly handled by the MKL class,
such as `kernel`, `C`, `epsilon`, `nu`, and `random_state`, which are managed by
the MKL estimator itself.

Example of setting SVM parameters::

    >>> from sklearn.mkl import MKLC
    >>> svm_params = {
    ...     'cache_size': 2000,
    ...     'max_iter': 100_000_000,
    ...     'decision_function_shape': 'ovo',
    ... }
    >>> mkl = MKLC(
    ...     kernels=["linear", "rbf"],
    ...     kernels_scopes=["single", "all"],
    ...     svm_params=svm_params
    ... )
    >>> mkl.fit(X=[[0, 0], [1, 1], [2, 2], [3, 3]], y=[0, 1, 2, 3])
    >>> mkl.decision_function([[0, 0], [1, 1]]).shape
    (2, 6)

Accessing the internal SVM
--------------------------

If you need to access the internal SVM for further inspection or to retrieve specific
attributes, you can do so by accessing the `mkl._svm` attribute. Keep in mind that
any modifications to this SVM instance will impact the behavior of the MKL class.

For example, after training the model, you can access the support vectors or the
fitted coefficients like this::

    >>> from sklearn.mkl import MKLR
    >>> X = [[-1, -1], [-2, -1], [1, 1], [2, 1]]
    >>> y = [-2.0, -3.0, 3.0, 4.0]
    >>> mkl = MKLR(kernels=["rbf"], kernels_scopes=["single"], epsilon=0.01)
    >>> mkl.fit(X, y)
    >>> # Access the internal SVM after fitting
    >>> internal_svm = mkl._svm
    >>> internal_svm.dual_coef_
    array([[-1., -1.,  1.,  1.]])
