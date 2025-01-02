.. _array_api:

================================
Array API support (experimental)
================================

.. currentmodule:: sklearn

The `Array API <https://data-apis.org/array-api/latest/>`_ specification defines
a standard API for all array manipulation libraries with a NumPy-like API.
Scikit-learn's Array API support requires
`array-api-compat <https://github.com/data-apis/array-api-compat>`__ to be installed,
and the environment variable `SCIPY_ARRAY_API` must be set to `1` before importing
`scipy` and `scikit-learn`:

.. prompt:: bash $

   export SCIPY_ARRAY_API=1

Please note that this environment variable is intended for temporary use.
For more details, refer to SciPy's `Array API documentation
<https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html#using-array-api-standard-support>`_.


Some scikit-learn estimators that primarily rely on NumPy (as opposed to using
Cython) to implement the algorithmic logic of their `fit`, `predict` or
`transform` methods can be configured to accept any Array API compatible input
datastructures and automatically dispatch operations to the underlying namespace
instead of relying on NumPy.

At this stage, this support is **considered experimental** and must be enabled
explicitly as explained in the following.

.. note::
    Currently, only `array-api-strict`, `cupy`, and `PyTorch` are known to work
    with scikit-learn's estimators.

The following video provides an overview of the standard's design principles
and how it facilitates interoperability between array libraries:

- `Scikit-learn on GPUs with Array API <https://www.youtube.com/watch?v=c_s8tr1AizA>`_
  by :user:`Thomas Fan <thomasjpfan>` at PyData NYC 2023.

Example usage
=============

Here is an example code snippet to demonstrate how to use `CuPy
<https://cupy.dev/>`_ to run
:class:`~discriminant_analysis.LinearDiscriminantAnalysis` on a GPU::

    >>> from sklearn.datasets import make_classification
    >>> from sklearn import config_context
    >>> from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    >>> import cupy

    >>> X_np, y_np = make_classification(random_state=0)
    >>> X_cu = cupy.asarray(X_np)
    >>> y_cu = cupy.asarray(y_np)
    >>> X_cu.device
    <CUDA Device 0>

    >>> with config_context(array_api_dispatch=True):
    ...     lda = LinearDiscriminantAnalysis()
    ...     X_trans = lda.fit_transform(X_cu, y_cu)
    >>> X_trans.device
    <CUDA Device 0>

After the model is trained, fitted attributes that are arrays will also be
from the same Array API namespace as the training data. For example, if CuPy's
Array API namespace was used for training, then fitted attributes will be on the
GPU. We provide a experimental `_estimator_with_converted_arrays` utility that
transfers an estimator attributes from Array API to a ndarray::

    >>> from sklearn.utils._array_api import _estimator_with_converted_arrays
    >>> cupy_to_ndarray = lambda array : array.get()
    >>> lda_np = _estimator_with_converted_arrays(lda, cupy_to_ndarray)
    >>> X_trans = lda_np.transform(X_np)
    >>> type(X_trans)
    <class 'numpy.ndarray'>

PyTorch Support
---------------

PyTorch Tensors are supported by setting `array_api_dispatch=True` and passing in
the tensors directly::

    >>> import torch
    >>> X_torch = torch.asarray(X_np, device="cuda", dtype=torch.float32)
    >>> y_torch = torch.asarray(y_np, device="cuda", dtype=torch.float32)

    >>> with config_context(array_api_dispatch=True):
    ...     lda = LinearDiscriminantAnalysis()
    ...     X_trans = lda.fit_transform(X_torch, y_torch)
    >>> type(X_trans)
    <class 'torch.Tensor'>
    >>> X_trans.device.type
    'cuda'

.. _array_api_supported:

Support for `Array API`-compatible inputs
=========================================

Estimators and other tools in scikit-learn that support Array API compatible inputs.

Estimators
----------

- :class:`decomposition.PCA` (with `svd_solver="full"`,
  `svd_solver="randomized"` and `power_iteration_normalizer="QR"`)
- :class:`linear_model.Ridge` (with `solver="svd"`)
- :class:`discriminant_analysis.LinearDiscriminantAnalysis` (with `solver="svd"`)
- :class:`preprocessing.KernelCenterer`
- :class:`preprocessing.LabelEncoder`
- :class:`preprocessing.MaxAbsScaler`
- :class:`preprocessing.MinMaxScaler`
- :class:`preprocessing.Normalizer`

Meta-estimators
---------------

Meta-estimators that accept Array API inputs conditioned on the fact that the
base estimator also does:

- :class:`model_selection.GridSearchCV`
- :class:`model_selection.RandomizedSearchCV`
- :class:`model_selection.HalvingGridSearchCV`
- :class:`model_selection.HalvingRandomSearchCV`

Metrics
-------

- :func:`sklearn.metrics.cluster.entropy`
- :func:`sklearn.metrics.accuracy_score`
- :func:`sklearn.metrics.d2_tweedie_score`
- :func:`sklearn.metrics.explained_variance_score`
- :func:`sklearn.metrics.f1_score`
- :func:`sklearn.metrics.fbeta_score`
- :func:`sklearn.metrics.max_error`
- :func:`sklearn.metrics.mean_absolute_error`
- :func:`sklearn.metrics.mean_absolute_percentage_error`
- :func:`sklearn.metrics.mean_gamma_deviance`
- :func:`sklearn.metrics.mean_pinball_loss`
- :func:`sklearn.metrics.mean_poisson_deviance` (requires `enabling array API support for SciPy <https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html#using-array-api-standard-support>`_)
- :func:`sklearn.metrics.mean_squared_error`
- :func:`sklearn.metrics.mean_squared_log_error`
- :func:`sklearn.metrics.mean_tweedie_deviance`
- :func:`sklearn.metrics.multilabel_confusion_matrix`
- :func:`sklearn.metrics.pairwise.additive_chi2_kernel`
- :func:`sklearn.metrics.pairwise.chi2_kernel`
- :func:`sklearn.metrics.pairwise.cosine_similarity`
- :func:`sklearn.metrics.pairwise.cosine_distances`
- :func:`sklearn.metrics.pairwise.euclidean_distances` (see :ref:`device_support_for_float64`)
- :func:`sklearn.metrics.pairwise.linear_kernel`
- :func:`sklearn.metrics.pairwise.paired_cosine_distances`
- :func:`sklearn.metrics.pairwise.paired_euclidean_distances`
- :func:`sklearn.metrics.pairwise.polynomial_kernel`
- :func:`sklearn.metrics.pairwise.rbf_kernel` (see :ref:`device_support_for_float64`)
- :func:`sklearn.metrics.pairwise.sigmoid_kernel`
- :func:`sklearn.metrics.precision_score`
- :func:`sklearn.metrics.precision_recall_fscore_support`
- :func:`sklearn.metrics.r2_score`
- :func:`sklearn.metrics.recall_score`
- :func:`sklearn.metrics.root_mean_squared_error`
- :func:`sklearn.metrics.root_mean_squared_log_error`
- :func:`sklearn.metrics.zero_one_loss`

Tools
-----

- :func:`model_selection.train_test_split`
- :func:`utils.check_consistent_length`

Coverage is expected to grow over time. Please follow the dedicated `meta-issue on GitHub
<https://github.com/scikit-learn/scikit-learn/issues/22352>`_ to track progress.

Type of return values and fitted attributes
-------------------------------------------

When calling functions or methods with Array API compatible inputs, the
convention is to return array values of the same array container type and
device as the input data.

Similarly, when an estimator is fitted with Array API compatible inputs, the
fitted attributes will be arrays from the same library as the input and stored
on the same device. The `predict` and `transform` method subsequently expect
inputs from the same array library and device as the data passed to the `fit`
method.

Note however that scoring functions that return scalar values return Python
scalars (typically a `float` instance) instead of an array scalar value.

Common estimator checks
=======================

Add the `array_api_support` tag to an estimator's set of tags to indicate that
it supports the Array API. This will enable dedicated checks as part of the
common tests to verify that the estimators result's are the same when using
vanilla NumPy and Array API inputs.

To run these checks you need to install
`array_api_compat <https://github.com/data-apis/array-api-compat>`_ in your
test environment. To run the full set of checks you need to install both
`PyTorch <https://pytorch.org/>`_ and `CuPy <https://cupy.dev/>`_ and have
a GPU. Checks that can not be executed or have missing dependencies will be
automatically skipped. Therefore it's important to run the tests with the
`-v` flag to see which checks are skipped:

.. prompt:: bash $

    pip install array-api-compat  # and other libraries as needed
    pytest -k "array_api" -v

.. _mps_support:

Note on MPS device support
--------------------------

On macOS, PyTorch can use the Metal Performance Shaders (MPS) to access
hardware accelerators (e.g. the internal GPU component of the M1 or M2 chips).
However, the MPS device support for PyTorch is incomplete at the time of
writing. See the following github issue for more details:

- https://github.com/pytorch/pytorch/issues/77764

To enable the MPS support in PyTorch, set the environment variable
`PYTORCH_ENABLE_MPS_FALLBACK=1` before running the tests:

.. prompt:: bash $

    PYTORCH_ENABLE_MPS_FALLBACK=1 pytest -k "array_api" -v

At the time of writing all scikit-learn tests should pass, however, the
computational speed is not necessarily better than with the CPU device.

.. _device_support_for_float64:

Note on device support for ``float64``
--------------------------------------

Certain operations within scikit-learn will automatically perform operations
on floating-point values with `float64` precision to prevent overflows and ensure
correctness (e.g., :func:`metrics.pairwise.euclidean_distances`). However,
certain combinations of array namespaces and devices, such as `PyTorch on MPS`
(see :ref:`mps_support`) do not support the `float64` data type. In these cases,
scikit-learn will revert to using the `float32` data type instead. This can result in
different behavior (typically numerically unstable results) compared to not using array
API dispatching or using a device with `float64` support.
