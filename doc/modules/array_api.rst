.. _array_api:

================================
Array API support (experimental)
================================

.. currentmodule:: sklearn

The `array API <https://data-apis.org/array-api/latest/>`__ specification defines
a standard API for all array manipulation libraries with a NumPy-like API.

Some scikit-learn estimators that primarily rely on NumPy (as opposed to using
Cython) to implement the algorithmic logic of their `fit`, `predict` or
`transform` methods can be configured to accept any Array API compatible input
data structures and automatically dispatch operations to the underlying namespace
instead of relying on NumPy.

At this stage, this support is **considered experimental**, must be enabled
explicitly by the `array_api_dispatch` configuration and assumes that the latest
versions libraries are installed. See below for details.

The following video provides an overview of the standard's design principles
and how it facilitates interoperability between array libraries:

- `Scikit-learn on GPUs with Array API <https://www.youtube.com/watch?v=c_s8tr1AizA>`__
  by :user:`Thomas Fan <thomasjpfan>` at PyData NYC 2023.

Supported array libraries
=========================
The following table lists the libraries and hardware for which we run automated
compliance tests on a regular basis. Other array API conforming libraries and
devices might also work out of the box.

.. list-table::
   :header-rows: 1
   :widths: 15 10 20 30

   * - Library
     - Install
     - Supported Devices
     - Notes
   * - `PyTorch <https://pytorch.org/>`_
     - `Install <https://pytorch.org/get-started/locally/>`_
     - CPU, NVIDIA GPU (CUDA), Apple GPU (MPS), Intel GPU (XPU)
     - See :ref:`mps_support`; see :ref:`xpu_support`;
       see :ref:`device_support_for_float64`
   * - `CuPy <https://cupy.dev/>`__
     - `Install <https://docs.cupy.dev/en/stable/install.html>`__
     - CUDA
     -
   * - `dpnp <https://intelpython.github.io/dpnp/>`__
     - `Install <https://intelpython.github.io/dpnp/quick_start_guide.html>`__
     - CPU, Intel GPU (SYCL)
     - See install link for driver setup; see :ref:`xpu_support`;
       see :ref:`device_support_for_float64`

Coverage is expected to grow over time.

Enabling array API support
==========================

The configuration parameter `array_api_dispatch` needs to be set to `True` to enable array
API support. We recommend setting this configuration globally to ensure consistent
behaviour and prevent accidental mixing of array namespaces.
Note that in the examples below, we use a context manager (:func:`config_context`)
to avoid having to reset it to `False` at the end of every code snippet, so as to
not affect the rest of the documentation.

Scikit-learn's support for the array API standard requires the environment variable
`SCIPY_ARRAY_API` to be set to `1` before importing `scipy` and `scikit-learn`:

.. code-block:: console

   $ export SCIPY_ARRAY_API=1

Please note that this environment variable is intended for temporary use.
For more details, refer to SciPy's `Array API documentation
<https://docs.scipy.org/doc/scipy/dev/api-dev/array_api.html#using-array-api-standard-support>`__.

The array API functionality assumes that the latest versions of scikit-learn's dependencies are
installed. Older versions might work, but we make no promises. While array API support is marked
as experimental, backwards compatibility is not guaranteed. In particular, when a newer version
of a dependency fixes a bug we will not introduce additional code to backport the fix or
maintain compatibility with older versions.

Scikit-learn accepts :term:`array-like` inputs for all :mod:`metrics`
and some estimators. When `array_api_dispatch=False`, these inputs are converted
into NumPy arrays using :func:`numpy.asarray` (or :func:`numpy.array`).
While this will successfully convert some array API inputs (e.g., JAX array),
we generally recommend setting `array_api_dispatch=True` when using array API inputs.
This is because NumPy conversion can often fail, e.g., torch tensor allocated on GPU.
Also note that the output array will be NumPy when `array_api_dispatch=False` whereas
when `array_api_dispatch=True` the output array will depend on the input array (see
:ref:`input_output_array_api` for details).

Example usage
=============

The example code snippet below demonstrates how to use `PyTorch
<https://pytorch.org/>`__ to run
:class:`~linear_model.LogisticRegression` on a CUDA GPU::

    >>> from sklearn.datasets import make_classification
    >>> from sklearn import config_context
    >>> from sklearn.linear_model import LogisticRegression
    >>> import torch

    >>> X_np, y_np = make_classification(random_state=0)
    >>> X_torch = torch.asarray(X_np, device="cuda", dtype=torch.float32)
    >>> y_torch = torch.asarray(y_np, device="cuda")

    >>> with config_context(array_api_dispatch=True):
    ...     clf = LogisticRegression().fit(X_torch, y_torch)
    ...     y_pred = clf.predict(X_torch)
    >>> type(y_pred)
    <class 'torch.Tensor'>
    >>> y_pred.device.type
    'cuda'

This pattern works identically with any supported array library. For example,
replace ``torch.asarray(..., device="cuda")`` with ``cupy.asarray(...)`` for CuPy
or ``dpnp.asarray(...)`` for dpnp. You can also target different devices within
PyTorch by changing the ``device=`` argument (e.g., ``"cpu"``, ``"xpu"``, ``"mps"``).

The computations of ``fit`` ran on the GPU, but the fitted attributes are stored
as NumPy arrays, whatever the array library and device of the training data::

    >>> type(clf.coef_)
    <class 'numpy.ndarray'>

``predict``, ``transform`` and other methods run on the array library and device
of the data they are called with, and return arrays of that library on that device.
The same fitted estimator can therefore be used with NumPy input::

    >>> with config_context(array_api_dispatch=True):
    ...     y_pred = clf.predict(X_np)
    >>> type(y_pred)
    <class 'numpy.ndarray'>

See :ref:`array_api_fitted_attributes` for details.

.. TODO(#34604): remove this note once `_NUMPY_FITTED_ATTRS` is removed.

.. note::
   Storing fitted attributes as NumPy arrays is being rolled out estimator by
   estimator. At the moment only :class:`~linear_model.LogisticRegression` and
   :class:`~linear_model.LogisticRegressionCV` do so. All other estimators store
   their fitted attributes in the array library and on the device of the
   training data, and can require the data passed to ``predict``, ``transform``
   and other methods to use the same array library and device.

.. _array_api_avoid_transfers:

Avoiding repeated transfers at inference time
---------------------------------------------

Each call to ``predict``, ``transform`` or another method copies the fitted
attributes it needs to the device of its input. A copy is reused within that
call, but not across calls. For a few calls on large inputs this copy is
negligible. When an estimator is called many times on small inputs, or has
large fitted attributes, the repeated copies can dominate the run time.

In that case, use ``move_estimator_to`` to store the fitted attributes on the
target array library and device once. Pass it the namespace and device of the
data you will call the estimator with::

    >>> from sklearn.utils._array_api import get_namespace_and_device
    >>> from sklearn.utils._array_api import move_estimator_to
    >>> with config_context(array_api_dispatch=True):
    ...     xp, _, device = get_namespace_and_device(X_torch)
    ...     clf_cuda = move_estimator_to(clf, xp, device)
    ...     y_pred = clf_cuda.predict(X_torch)
    >>> clf_cuda.coef_.device.type
    'cuda'

The moved estimator still accepts input from any array library and device, but
input from the library and device it was moved to no longer needs a copy.
Calling ``fit`` on the moved estimator stores its fitted attributes as NumPy
arrays again. Note that a moved estimator can only be unpickled in an
environment where its array library and device are available.
``move_estimator_to`` does not yet support estimators with string class labels
(see :issue:`35109`).

.. _array_api_supported:

Support for array API compatible inputs
=======================================

Estimators and other tools in scikit-learn that support array API compatible inputs.

Estimators
----------

- :class:`covariance.LedoitWolf` (see :ref:`device_support_for_float64`)
- :class:`covariance.OAS` (see :ref:`device_support_for_float64`)
- :class:`decomposition.PCA` (with `svd_solver="full"`, `svd_solver="covariance_eigh"`, or
  `svd_solver="randomized"` (`svd_solver="randomized"` only if `power_iteration_normalizer="QR"`))
- :class:`kernel_approximation.Nystroem`
- :class:`linear_model.LogisticRegression` (with `solver="lbfgs"` and `solver="newton-cg"`)
- :class:`linear_model.LogisticRegressionCV` (with `solver="lbfgs"` and `solver="newton-cg"`)
- :class:`linear_model.PoissonRegressor` (with `solver="lbfgs"`)
- :class:`linear_model.Ridge` (with `solver="svd"`)
- :class:`linear_model.RidgeCV` (see :ref:`device_support_for_float64`)
- :class:`linear_model.RidgeClassifier` (with `solver="svd"`)
- :class:`linear_model.RidgeClassifierCV` (see :ref:`device_support_for_float64`)
- :class:`discriminant_analysis.LinearDiscriminantAnalysis` (with `solver="svd"`)
- :class:`naive_bayes.GaussianNB`
- :class:`preprocessing.Binarizer`
- :class:`preprocessing.KernelCenterer`
- :class:`preprocessing.LabelBinarizer` (with `sparse_output=False`)
- :class:`preprocessing.LabelEncoder`
- :class:`preprocessing.MaxAbsScaler`
- :class:`preprocessing.MinMaxScaler`
- :class:`preprocessing.Normalizer`
- :class:`preprocessing.PolynomialFeatures`
- :class:`preprocessing.StandardScaler` (see :ref:`device_support_for_float64`)
- :class:`mixture.GaussianMixture` (with `init_params="random"` or
  `init_params="random_from_data"` and `warm_start=False`)

Meta-estimators
---------------

Meta-estimators that accept Array API inputs conditioned on the fact that the
base estimator also does:

- :class:`calibration.CalibratedClassifierCV` (with `method="temperature"`)
- :class:`pipeline.FeatureUnion`
- :class:`model_selection.GridSearchCV`
- :class:`model_selection.RandomizedSearchCV`
- :class:`model_selection.HalvingGridSearchCV`
- :class:`model_selection.HalvingRandomSearchCV`

Metrics
-------

- :func:`sklearn.metrics.accuracy_score`
- :func:`sklearn.metrics.average_precision_score`
- :func:`sklearn.metrics.balanced_accuracy_score`
- :func:`sklearn.metrics.brier_score_loss`
- :func:`sklearn.metrics.calinski_harabasz_score`
- :func:`sklearn.metrics.cluster.contingency_matrix`
- :func:`sklearn.metrics.cohen_kappa_score`
- :func:`sklearn.metrics.confusion_matrix`
- :func:`sklearn.metrics.d2_absolute_error_score`
- :func:`sklearn.metrics.d2_brier_score`
- :func:`sklearn.metrics.d2_log_loss_score`
- :func:`sklearn.metrics.d2_pinball_score`
- :func:`sklearn.metrics.d2_tweedie_score`
- :func:`sklearn.metrics.det_curve`
- :func:`sklearn.metrics.explained_variance_score`
- :func:`sklearn.metrics.f1_score`
- :func:`sklearn.metrics.fbeta_score`
- :func:`sklearn.metrics.hamming_loss`
- :func:`sklearn.metrics.jaccard_score`
- :func:`sklearn.metrics.log_loss`
- :func:`sklearn.metrics.matthews_corrcoef` (see :ref:`device_support_for_float64`)
- :func:`sklearn.metrics.max_error`
- :func:`sklearn.metrics.mean_absolute_error`
- :func:`sklearn.metrics.mean_absolute_percentage_error`
- :func:`sklearn.metrics.mean_gamma_deviance`
- :func:`sklearn.metrics.mean_pinball_loss`
- :func:`sklearn.metrics.mean_poisson_deviance`
- :func:`sklearn.metrics.mean_squared_error`
- :func:`sklearn.metrics.mean_squared_log_error`
- :func:`sklearn.metrics.mean_tweedie_deviance`
- :func:`sklearn.metrics.median_absolute_error`
- :func:`sklearn.metrics.multilabel_confusion_matrix`
- :func:`sklearn.metrics.pairwise.additive_chi2_kernel`
- :func:`sklearn.metrics.pairwise.chi2_kernel`
- :func:`sklearn.metrics.pairwise.cosine_similarity`
- :func:`sklearn.metrics.pairwise.cosine_distances`
- :func:`sklearn.metrics.pairwise_distances` (only supports "cosine", "euclidean", "manhattan" and "l2" metrics)
- :func:`sklearn.metrics.pairwise_distances_argmin`
- :func:`sklearn.metrics.pairwise.euclidean_distances` (see :ref:`device_support_for_float64`)
- :func:`sklearn.metrics.pairwise.laplacian_kernel`
- :func:`sklearn.metrics.pairwise.linear_kernel`
- :func:`sklearn.metrics.pairwise.manhattan_distances`
- :func:`sklearn.metrics.pairwise.paired_cosine_distances`
- :func:`sklearn.metrics.pairwise.paired_euclidean_distances`
- :func:`sklearn.metrics.pairwise.paired_manhattan_distances`
- :func:`sklearn.metrics.pairwise.pairwise_kernels`
- :func:`sklearn.metrics.pairwise.polynomial_kernel`
- :func:`sklearn.metrics.pairwise.rbf_kernel` (see :ref:`device_support_for_float64`)
- :func:`sklearn.metrics.pairwise.sigmoid_kernel`
- :func:`sklearn.metrics.precision_score`
- :func:`sklearn.metrics.precision_recall_curve`
- :func:`sklearn.metrics.precision_recall_fscore_support`
- :func:`sklearn.metrics.r2_score`
- :func:`sklearn.metrics.recall_score`
- :func:`sklearn.metrics.roc_curve`
- :func:`sklearn.metrics.root_mean_squared_error`
- :func:`sklearn.metrics.root_mean_squared_log_error`
- :func:`sklearn.metrics.zero_one_loss`

Tools
-----

- :func:`preprocessing.label_binarize` (with `sparse_output=False`)
- :func:`model_selection.cross_val_predict`
- :func:`model_selection.train_test_split`
- :func:`utils.check_consistent_length`

.. _input_output_array_api:

Input and output array type handling
====================================

Estimators and scoring functions are able to accept input arrays
from different array libraries and/or devices. When a mixed set of input arrays is
passed, scikit-learn converts arrays as needed to make them all consistent. The
following section describes how input and output array types are handled
when `array_api_dispatch=True`.

For estimators, the rule is **"everything follows** `X` **"** - mixed array inputs are
converted so that they all match the array library and device of `X`.
For scoring functions the rule is **"everything follows** `y_pred` **"** - mixed array
inputs are converted so that they all match the array library and device of `y_pred`.

Mixed array input support also covers the particular case where `y` is a NumPy
array of string values, while other inputs are numerical arrays of any
container type. As all models require numerical input, scikit-learn often converts
`y` to a numerical representation internally (for instance via one-hot encoding or
ordinal encoding). The result of this is able to be moved to the namespace and device
of the other inputs, as string-valued arrays are not covered by the array API
specifications.

Note that the `classes_` estimator attribute is a NumPy array, like all other
fitted attributes, which allows it to hold string class labels (see
:ref:`array_api_fitted_attributes`).

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

.. _array_api_fitted_attributes:

Fitted attributes
-----------------

The computations of `fit` run on the array library and device of `X`, but
fitted attributes that are arrays are stored as NumPy arrays, whatever the
array library and device of the data passed to `fit`. This includes `classes_`.

Methods such as `predict`, `predict_proba`, `decision_function` and `transform`
run on the array library and device of the `X` they are called with, and return
arrays from that library on that device. The fitted attributes they need are
copied to that library and device for the duration of the call. For
classifiers, `predict` returns numerical class labels in the array library and
on the device of `X`, and string class labels as a NumPy array.

As a result:

- fitted attributes can be inspected or plotted with NumPy based tools, without
  first converting them;
- an estimator fitted on one device can be used on another, for example fitted
  on a GPU and then used for prediction on a CPU, or the other way around;
- an estimator fitted on a GPU can be unpickled on a machine without a GPU, or
  without the array library it was fitted with.

The copies made at every call can be avoided, see
:ref:`array_api_avoid_transfers`.

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

Information for developers
==========================

See :ref:`array_api_developers` for how to add array API support to an
estimator and how to test it.

Notes
=====

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

.. code-block:: console

    $ PYTORCH_ENABLE_MPS_FALLBACK=1 pytest -k "array_api" -v

At the time of writing all scikit-learn tests should pass, however, the
computational speed is not necessarily better than with the CPU device.

.. _device_support_for_float64:

Note on device support for ``float64``
--------------------------------------

Certain operations within scikit-learn will automatically perform operations
on floating-point values with `float64` precision to prevent overflows and ensure
correctness (e.g., :func:`metrics.pairwise.euclidean_distances`,
:class:`preprocessing.StandardScaler`). However, certain combinations of array
namespaces and devices, such as PyTorch on MPS (see :ref:`mps_support`) and
some Intel GPU devices with PyTorch XPU or dpnp, do not support the `float64`
data type. In these cases, scikit-learn will revert to using the `float32` data
type instead. This can result in different behavior (typically numerically
unstable results) compared to not using array API dispatching or using a device
with `float64` support.

.. _xpu_support:

Note on Intel GPU support
-------------------------

PyTorch XPU support is only available in PyTorch >= 2.12. For detailed compatibility
information and setup instructions, see the `PyTorch XPU documentation
<https://docs.pytorch.org/docs/2.12/notes/get_start_xpu.html>`__.
