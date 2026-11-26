.. _array_api_developers:

========================================
Adding array API support to an estimator
========================================

.. currentmodule:: sklearn

This section is intended for contributors who add or maintain array API support
in scikit-learn estimators. For how to use array API support, see the
:ref:`user guide <array_api>`.

The helpers described here are private (they live in
``sklearn/utils/_array_api.py``) and can change without deprecation.

Overview
========

An estimator that supports the array API:

- runs the computations of ``fit`` in the array library (namespace) and on the
  device of ``X``;
- stores its fitted attributes as NumPy arrays;
- runs ``predict``, ``transform`` and other inference methods in the namespace
  and on the device of the ``X`` they are called with, and returns arrays from
  that namespace on that device.

See :ref:`array_api_fitted_attributes` for this contract from the point of view
of users.

.. TODO(#34604): remove this note once `_NUMPY_FITTED_ATTRS` is removed.

.. note::
   Storing fitted attributes as NumPy arrays is being rolled out estimator by
   estimator. Only estimators whose class name is listed in
   ``_NUMPY_FITTED_ATTRS`` follow the contract described here. The other
   estimators store their fitted attributes in the namespace of ``X``, and their
   inference methods can require ``X`` to use the same namespace and device.
   Some of them check this upfront with ``check_same_namespace``. Those that
   share inference code with a converted estimator, for example the subclasses
   of ``LinearModel`` such as ``Ridge``, already read their attributes through
   ``_fitted_attrs_like``: they accept ``X`` from any namespace and move their
   attributes on every call.

Fitting
=======

Fitted attributes are converted to NumPy by the ``_fit_context`` decorator
(``sklearn/base.py``) after the decorated method returns. Methods decorated with
it, usually ``fit`` and ``partial_fit``, therefore do not need any extra code,
and every ``return`` statement of such a method is covered. A method that sets
fitted attributes without being decorated, or without calling a decorated
method, does not get its attributes converted.

The conversion only looks at the estimator's own attributes and the containers
(``dict``, ``list``, ``tuple``, ...) stored in them. Nested estimators are left
alone, because each of them converts its own attributes when its own decorated
fit method returns.

An attribute that holds arrays but is neither an array, a container nor an
estimator can define ``__sklearn_array_api_convert__(self, converter)``, which
returns a copy of the object with ``converter`` applied to each of its arrays.

Inference methods
=================

In inference methods, access fitted attributes through ``_fitted_attrs_like``,
which returns a read-only view of the estimator. Reading an attribute from the
view moves it to the namespace and device of ``X``::

    from sklearn.utils._array_api import _fitted_attrs_like

    def decision_function(self, X):
        check_is_fitted(self)
        X = validate_data(self, X, reset=False)
        fitted = _fitted_attrs_like(self, X)
        return X @ fitted.coef_.T + fitted.intercept_

If the namespace and device of ``X`` are already known, pass them as well to
avoid inspecting ``X`` again::

    xp, _, device = get_namespace_and_device(X)
    fitted = _fitted_attrs_like(self, X, xp=xp, device=device)

Keep the following in mind:

- If ``X`` is not a NumPy array and has a floating point dtype, floating point
  attributes are also cast to the dtype of ``X``. Namespaces like PyTorch do not
  promote mixed float dtypes, so an estimator fitted on float64 NumPy data would
  otherwise fail on float32 GPU data. NumPy inputs are not cast: NumPy promotes
  mixed float dtypes, so their output dtype does not depend on whether array API
  dispatch is enabled.

- Read attributes from ``self`` when only their metadata is needed (for example
  ``self.coef_.ndim`` or ``self.classes_.shape[0]``). Metadata does not depend
  on the namespace, and reading it from the view would copy the array for
  nothing.
- A moved attribute is cached in the view, so reading it a second time from the
  same view does not copy it again. The cache lives as long as the view. An
  inference method that calls another inference method (for example ``predict``
  calling ``decision_function``) copies the attributes that both of them use
  twice.
- The view cannot be assigned to. Assign to the estimator instead.

Code shared between fitting and inference
-----------------------------------------

Some computations that read fitted attributes run both during fitting, when the
attributes are still in the namespace of ``X``, and in inference methods, when
they are NumPy arrays. Do not pass "the estimator or a view of it" around to
handle both cases. Instead:

- Where possible, write the computation as a function of arrays, and let each
  caller pass the arrays: ``self.precision_`` from fit or from a method without
  ``X`` such as ``get_precision``, ``fitted.precision_`` from an inference
  method.
- Otherwise, for example in a method that subclasses override, let the method
  build its own view from the ``X`` it computes with,
  ``fitted = _fitted_attrs_like(self, X)``. During fitting the attributes
  already follow ``X``, so the view passes them through unchanged.

Arrays used during fitting
==========================

Fitted state that is read back during fitting is a NumPy array, even if the
previous fit was done with a different namespace. With ``warm_start=True`` or in
``partial_fit``, move the stored state to the namespace and device used for the
computation with ``_move_state`` before using it. The same applies to other
arrays that a fit method computes with and writes back to the estimator, for
example array-valued constructor parameters used as initial values::

    from sklearn.utils._array_api import _move_state

    # Running statistics, kept in a higher precision than `X`.
    self.mean_, self.var_ = _move_state(self.mean_, self.var_, xp=xp, device=device)
    # Parameters that are combined with `X` in the dtype of `X`.
    means = _move_state(self.means_init, xp=xp, device=device, dtype=X.dtype)

``_move_state`` accepts containers of arrays and leaves values that are not
arrays unchanged. Pass ``dtype`` to cast real floating point arrays: use the
dtype that the estimator would use for this computation in a fresh fit, which is
usually the floating dtype of ``X``. Leave it out to keep the stored dtype, for
example for accumulators that the estimator deliberately keeps in a higher
precision than ``X``. ``dtype`` is ignored when the target namespace is NumPy, for
the same reason that ``_fitted_attrs_like`` does not cast for NumPy ``X``: NumPy
promotes mixed float dtypes, so NumPy inputs keep the behaviour they have without
array API dispatch. Call sites can therefore always pass the dtype of a fresh fit,
for example ``dtype=X.dtype``, without checking the namespace. Cast explicitly
where a NumPy computation needs a specific dtype.

The target is where the computation happens, which is not always the namespace
of ``X``: ``LogisticRegression`` moves its stored coefficients to NumPy for its
solvers, and ``classes_`` is combined with ``y``.

Use ``_fitted_attrs_like`` to read attributes that an inference method computes
with, and ``_move_state`` for arrays whose results are written back to the
estimator. The view casts with a fixed rule that is right for outputs that follow
``X`` but not for every kind of state, and it caches the arrays it moved, so
reading an attribute through the same view after assigning to the estimator
returns the old value.

``classes_`` is a NumPy array, which allows it to hold string class labels that
cannot be represented in every namespace. Predicted labels follow ``X``, unless
they cannot be represented in the namespace of ``X``, in which case they are a
NumPy array. Which one applies depends only on the dtype of ``classes_``, so a
fitted classifier always returns the same type. To turn scores into class
labels, compute the indices in the namespace of ``X`` and pass them to
``_take_classes``::

    from sklearn.utils._array_api import _take_classes

    indices = xp.argmax(scores, axis=1)
    return _take_classes(self.classes_, indices, xp=xp, device=device)

To compare labels with ``classes_``, for example to encode ``y`` or to check that
its labels are known, move ``classes_`` to the namespace and device of the labels
with ``move_to(self.classes_, xp=xp, device=device)``. Code outside the estimator
that combines ``classes_`` with labels, for example a metric called with
``labels=estimator.classes_``, does the same.

Labels are compared by value, so in both directions they keep the dtype of
``classes_``: they are never cast to the dtype of ``X`` or ``y``. The only
exception is a device without float64 support, to which ``move_to`` moves
float64 ``classes_`` as float32. For the same reason ``classes_`` cannot be read
through ``_fitted_attrs_like``, which casts floating point attributes: doing so
raises a ``TypeError``.

Meta-estimators and nested fits
===============================

.. TODO: update this section once the conversion of nested fits is settled,
   see the discussion in https://github.com/scikit-learn/scikit-learn/pull/35053

The attributes of a sub-estimator are converted to NumPy as soon as the
sub-estimator's own fit method returns. This also happens while the fit method
of the meta-estimator is still running. A meta-estimator that uses the fitted
attributes of its sub-estimators during ``fit`` therefore sees NumPy arrays,
even when ``X`` is from another namespace, and has to move them (or use
``_fitted_attrs_like``) before combining them with ``X``.

Testing
=======

Add the ``array_api_support`` tag to an estimator's set of tags to indicate that
it supports the array API. This enables dedicated checks as part of the common
tests. They verify that:

- the fitted attributes and the outputs of the inference methods are the same
  when using NumPy and array API inputs, including inputs where ``y`` and
  ``sample_weight`` use a different namespace or device than ``X``, and string
  ``y``;
- all fitted attributes that are arrays are NumPy arrays, and the fitted
  estimator does not refer to the array library it was fitted with, for example
  through a module or device object, so that it can be pickled and loaded where
  that library is not available;
- the outputs of the inference methods are in the namespace and on the device of
  their input;
- an estimator fitted with one namespace can be used for inference with another,
  and the other way around (``check_array_api_cross_namespace_inference``),
  also when the input has a different floating point dtype than the data used
  to fit (``check_array_api_cross_namespace_inference_mixed_dtypes``). These
  checks call ``predict``, ``predict_proba``, ``predict_log_proba``,
  ``predict_joint_log_proba``, ``decision_function``, ``transform``,
  ``inverse_transform``, ``score`` and ``score_samples``, so estimator-specific
  tests are only needed for other methods;
- ``partial_fit`` accepts a batch from one namespace after a batch from another,
  and gives the same fitted attributes as two NumPy batches
  (``check_array_api_cross_namespace_partial_fit``).

To run these checks you need to install
`array-api-strict <https://data-apis.org/array-api-strict/>`__ in your
test environment. This allows you to run checks without having a
GPU. To run checks on real GPU devices you also need to install
`PyTorch <https://pytorch.org/>`__, `CuPy <https://cupy.dev/>`__, and/or
`dpnp <https://intelpython.github.io/dpnp/>`__, and have compatible GPU
hardware. Full GPU coverage is expected to be split across machines because
CUDA, MPS, and Intel GPU backends require different hardware. Checks that can
not be executed or have missing dependencies will be automatically skipped.
Therefore it's important to run the tests with the `-v` flag to see which
checks are skipped:

.. code-block:: console

    $ pip install array-api-strict  # and other libraries as needed
    $ pytest -k "array_api" -v

On macOS, set ``PYTORCH_ENABLE_MPS_FALLBACK=1`` to run the checks on the MPS
device (see :ref:`mps_support`).

Running the scikit-learn tests against `array-api-strict` should help reveal
most code problems related to handling multiple device inputs via the use of
simulated non-CPU devices. This allows for fast iterative development and debugging of
array API related code.

However, to ensure full handling of PyTorch, CuPy, or dpnp inputs allocated on
actual GPU devices, it is necessary to run the tests against those libraries and
hardware. This can either be achieved by using
`Google Colab <https://gist.github.com/EdAbati/ff3bdc06bafeb92452b3740686cc8d7c>`__
for CUDA or leveraging our CI infrastructure on pull requests. CUDA and Intel
GPU tests are manually triggered by maintainers. Intel GPU testing for PyTorch
XPU and dpnp is run on a dedicated self-hosted runner:
`probabl-ai/scikit-learn-intel-workflow
<https://github.com/probabl-ai/scikit-learn-intel-workflow>`__.
