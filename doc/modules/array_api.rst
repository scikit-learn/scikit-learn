.. _array_api:

================================
Array API support (experimental)
================================

.. currentmodule:: sklearn

The `Array API <https://data-apis.org/array-api/latest/>`_ specification defines
a standard API for all array manipulation libraries with a NumPy-like API.
Scikit-learn vendors pinned copies of
`array-api-compat <https://github.com/data-apis/array-api-compat>`__
and `array-api-extra <https://github.com/data-apis/array-api-extra>`__.

Scikit-learn's support for the array API standard requires the environment variable
`SCIPY_ARRAY_API` to be set to `1` before importing `scipy` and `scikit-learn`:

.. prompt:: bash $

   export SCIPY_ARRAY_API=1

Please note that this environment variable is intended for temporary use.
For more details, refer to SciPy's `Array API documentation
<https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html#using-array-api-standard-support>`_.

Some scikit-learn estimators that primarily rely on NumPy (as opposed to using
Cython) to implement the algorithmic logic of their `fit`, `predict` or
`transform` methods can be configured to accept any Array API compatible input
data structures and automatically dispatch operations to the underlying namespace
instead of relying on NumPy.

At this stage, this support is **considered experimental** and must be enabled
explicitly by the `array_api_dispatch` configuration. See below for details.

.. note::
    Currently, only `array-api-strict`, `cupy`, and `PyTorch` are known to work
    with scikit-learn's estimators.

The following video provides an overview of the standard's design principles
and how it facilitates interoperability between array libraries:

- `Scikit-learn on GPUs with Array API <https://www.youtube.com/watch?v=c_s8tr1AizA>`_
  by :user:`Thomas Fan <thomasjpfan>` at PyData NYC 2023.

Example usage
=============

The configuration `array_api_dispatch=True` needs to be set to `True` to enable array
API support. We recommend setting this configuration globally to ensure consistent
behaviour and prevent accidental mixing of array namespaces.
Note that we set it with :func:`config_context` below to avoid having to call
:func:`set_config(array_api_dispatch=False)` at the end of every code snippet
that uses the array API.
The example code snippet below demonstrates how to use `CuPy
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
GPU. We provide an experimental `_estimator_with_converted_arrays` utility that
transfers an estimator attributes from Array API to a ndarray::

    >>> from sklearn.utils._array_api import _estimator_with_converted_arrays
    >>> cupy_to_ndarray = lambda array : array.get()
    >>> lda_np = _estimator_with_converted_arrays(lda, cupy_to_ndarray)
    >>> X_trans = lda_np.transform(X_np)
    >>> type(X_trans)
    <class 'numpy.ndarray'>

PyTorch Support
---------------

PyTorch Tensors can also be passed directly::

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
- :class:`preprocessing.Binarizer`
- :class:`preprocessing.KernelCenterer`
- :class:`preprocessing.LabelEncoder`
- :class:`preprocessing.MaxAbsScaler`
- :class:`preprocessing.MinMaxScaler`
- :class:`preprocessing.Normalizer`
- :class:`preprocessing.PolynomialFeatures`
- :class:`mixture.GaussianMixture` (with `init_params="random"` or
  `init_params="random_from_data"` and `warm_start=False`)

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

- :func:`sklearn.metrics.accuracy_score`
- :func:`sklearn.metrics.confusion_matrix`
- :func:`sklearn.metrics.d2_tweedie_score`
- :func:`sklearn.metrics.explained_variance_score`
- :func:`sklearn.metrics.f1_score`
- :func:`sklearn.metrics.fbeta_score`
- :func:`sklearn.metrics.hamming_loss`
- :func:`sklearn.metrics.jaccard_score`
- :func:`sklearn.metrics.max_error`
- :func:`sklearn.metrics.mean_absolute_error`
- :func:`sklearn.metrics.mean_absolute_percentage_error`
- :func:`sklearn.metrics.mean_gamma_deviance`
- :func:`sklearn.metrics.mean_pinball_loss`
- :func:`sklearn.metrics.mean_poisson_deviance` (requires `enabling array API support for SciPy <https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html#using-array-api-standard-support>`_)
- :func:`sklearn.metrics.mean_squared_error`
- :func:`sklearn.metrics.mean_squared_log_error`
- :func:`sklearn.metrics.mean_tweedie_deviance`
- :func:`sklearn.metrics.median_absolute_error`
- :func:`sklearn.metrics.multilabel_confusion_matrix`
- :func:`sklearn.metrics.pairwise.additive_chi2_kernel`
- :func:`sklearn.metrics.pairwise.chi2_kernel`
- :func:`sklearn.metrics.pairwise.cosine_similarity`
- :func:`sklearn.metrics.pairwise.cosine_distances`
- :func:`sklearn.metrics.pairwise.pairwise_distances` (only supports "cosine", "euclidean" and "l2" metrics)
- :func:`sklearn.metrics.pairwise.euclidean_distances` (see :ref:`device_support_for_float64`)
- :func:`sklearn.metrics.pairwise.linear_kernel`
- :func:`sklearn.metrics.pairwise.paired_cosine_distances`
- :func:`sklearn.metrics.pairwise.paired_euclidean_distances`
- :func:`sklearn.metrics.pairwise.pairwise_kernels` (supports all `sklearn.pairwise.PAIRWISE_KERNEL_FUNCTIONS` except :func:`sklearn.metrics.pairwise.laplacian_kernel`)
- :func:`sklearn.metrics.pairwise.polynomial_kernel`
- :func:`sklearn.metrics.pairwise.rbf_kernel` (see :ref:`device_support_for_float64`)
- :func:`sklearn.metrics.pairwise.sigmoid_kernel`
- :func:`sklearn.metrics.precision_score`
- :func:`sklearn.metrics.precision_recall_fscore_support`
- :func:`sklearn.metrics.r2_score`
- :func:`sklearn.metrics.recall_score`
- :func:`sklearn.metrics.roc_curve`
- :func:`sklearn.metrics.root_mean_squared_error`
- :func:`sklearn.metrics.root_mean_squared_log_error`
- :func:`sklearn.metrics.zero_one_loss`

Tools
-----

- :func:`model_selection.train_test_split`
- :func:`utils.check_consistent_length`

Coverage is expected to grow over time. Please follow the dedicated `meta-issue on GitHub
<https://github.com/scikit-learn/scikit-learn/issues/22352>`_ to track progress.

Input and output array type handling
====================================

Estimators and scoring functions are able to accept input arrays
from different array libraries and/or devices. When a mixed set of input arrays is
passed, scikit-learn converts arrays as needed to make them all consistent.

For estimators, the rule is **"everything follows `X`"** - mixed array inputs are
converted so that they all match the array library and device of `X`.
For scoring functions the rule is **"everything follows `y_pred`"** - mixed array
inputs are converted so that they all match the array library and device of `y_pred`.

When a function or method has been called with array API compatible inputs, the
convention is to return arrays from the same array library and on the same
device as the input data.

Estimators
----------

When an estimator is fitted with an array API compatible `X`, all other
array inputs, including constructor arguments, (e.g., `y`, `sample_weight`)
will be converted to match the array library and device of `X`, if they do not already.
This behaviour enables switching from processing on the CPU to processing
on the GPU at any point within a pipeline.

This allows estimators to accept mixed input types, enabling `X` to be moved
to a different device within a pipeline, without explicitly moving `y`.
Note that scikit-learn pipelines do not allow transformation of `y` (to avoid
:ref:`leakage <data_leakage>`).

Take for example a pipeline where `X` and `y` both start on CPU, and go through
the following three steps:

* :class:`~sklearn.preprocessing.TargetEncoder`, which will transform categorial
  `X` but also requires `y`, meaning both `X` and `y` need to be on CPU.
* :class:`FunctionTransformer(func=partial(torch.asarray, device="cuda")) <sklearn.preprocessing.FunctionTransformer>`,
  which moves `X` to GPU, to improve performance in the next step.
* :class:`~sklearn.linear_model.Ridge`, whose performance can be improved when
  passed arrays on a GPU, as they can handle large matrix operations very efficiently.

`X` initially contains categorical string data (thus needs to be on CPU), which is
target encoded to numerical values in :class:`~sklearn.preprocessing.TargetEncoder`.
`X` is then explicitly moved to GPU to improve the performance of
:class:`~sklearn.linear_model.Ridge`. `y` cannot be transformed by the pipeline
(recall scikit-learn pipelines do not allow transformation of `y`) but as
:class:`~sklearn.linear_model.Ridge` is able to accept mixed input types,
this is not a problem and the pipeline is able to be run.

The fitted attributes of an estimator fitted with an array API compatible `X`, will
be arrays from the same library as the input and stored on the same device.
The `predict` and `transform` method subsequently expect
inputs from the same array library and device as the data passed to the `fit`
method.

Scoring functions
-----------------

When an array API compatible `y_pred` is passed to a scoring function,
all other array inputs (e.g., `y_true`, `sample_weight`) will be converted
to match the array library and device of `y_pred`, if they do not already.
This allows scoring functions to accept mixed input types, enabling them to be
used within a :term:`meta-estimator` (or function that accepts estimators), with a
pipeline that moves input arrays between devices (e.g., CPU to GPU).

For example, to be able to use the pipeline described above within e.g.,
:func:`~sklearn.model_selection.cross_validate` or
:class:`~sklearn.model_selection.GridSearchCV`, the scoring function internally
called needs to be able to accept mixed input types.

The output type of scoring functions depends on the number of output values.
When a scoring function returns a scalar value, it will return a Python
scalar (typically a `float` instance) instead of an array scalar value.
For scoring functions that support :term:`multiclass` or :term:`multioutput`,
an array from the same array library and device as `y_pred` will be returned when
multiple values need to be output.

Common estimator checks
=======================

Add the `array_api_support` tag to an estimator's set of tags to indicate that
it supports the array API. This will enable dedicated checks as part of the
common tests to verify that the estimators' results are the same when using
vanilla NumPy and array API inputs.

To run these checks you need to install
`array-api-strict <https://data-apis.org/array-api-strict/>`_ in your
test environment. This allows you to run checks without having a
GPU. To run the full set of checks you also need to install
`PyTorch <https://pytorch.org/>`_, `CuPy <https://cupy.dev/>`_ and have
a GPU. Checks that can not be executed or have missing dependencies will be
automatically skipped. Therefore it's important to run the tests with the
`-v` flag to see which checks are skipped:

.. prompt:: bash $

    pip install array-api-strict  # and other libraries as needed
    pytest -k "array_api" -v

Running the scikit-learn tests against `array-api-strict` should help reveal
most code problems related to handling multiple device inputs via the use of
simulated non-CPU devices. This allows for fast iterative development and debugging of
array API related code.

However, to ensure full handling of PyTorch or CuPy inputs allocated on actual GPU
devices, it is necessary to run the tests against those libraries and hardware.
This can either be achieved by using
`Google Colab <https://gist.github.com/EdAbati/ff3bdc06bafeb92452b3740686cc8d7c>`_
or leveraging our CI infrastructure on pull requests (manually triggered by maintainers
for cost reasons).

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
